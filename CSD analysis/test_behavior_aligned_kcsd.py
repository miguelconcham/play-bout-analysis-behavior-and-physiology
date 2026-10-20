import csv
import os
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from behavior_aligned_kcsd import (
    average_kcsd_columns,
    average_saved_behavior_kcsd_from_root,
    discover_behavior_aligned_pairs,
    epoch_paired_with_onset_baseline,
    estimate_animal_pooled_kcsd,
    load_paired_raw,
    organize_animal_pooled_figures,
    organize_behavior_aligned_figures,
    run_behavior_aligned_kcsd,
    trial_weighted_psth,
)


def _write_leaf(
    folder,
    lfp,
    *,
    t0=5,
    fs=1000.0,
    event_ids=None,
):
    os.makedirs(folder, exist_ok=True)
    np.save(os.path.join(folder, "lfp.npy"), lfp)
    np.save(os.path.join(folder, "x_coords.npy"), np.arange(lfp.shape[1]))
    np.save(os.path.join(folder, "y_coords.npy"), np.arange(lfp.shape[1]) * 20)
    np.save(os.path.join(folder, "T0.npy"), np.array(t0))
    np.save(os.path.join(folder, "fs.npy"), np.array(fs))
    if event_ids is not None:
        np.save(os.path.join(folder, "event_ids.npy"), np.asarray(event_ids))


class PairedBaselineTests(unittest.TestCase):
    def test_offset_uses_onset_baseline(self):
        onset = np.array(
            [
                [[0, 0, 0, 1, 3, 10, 12, 0]],
                [[0, 0, 0, 4, 8, 20, 24, 0]],
            ],
            dtype=float,
        )
        # Offset pre-alignment samples have deliberately large values. They
        # must never be used as the baseline.
        offset = np.array(
            [
                [[0, 0, 0, 100, 200, 20, 22, 0]],
                [[0, 0, 0, 300, 500, 40, 44, 0]],
            ],
            dtype=float,
        )

        out = epoch_paired_with_onset_baseline(
            onset,
            offset,
            t0_onset=5,
            t0_offset=5,
            fs=1000.0,
            t_pre=2.0,
            t_post=2.0,
            step_ms=1.0,
            onset_baseline_window_ms=(-2.0, 0.0),
        )

        np.testing.assert_allclose(out["onset_baseline"][:, 0, 0], [2.0, 6.0])
        np.testing.assert_allclose(
            out["offset_epochs"][0, 0], [98.0, 198.0, 18.0, 20.0]
        )
        np.testing.assert_allclose(
            out["offset_epochs"][1, 0], [294.0, 494.0, 34.0, 38.0]
        )

        chunked = epoch_paired_with_onset_baseline(
            onset,
            offset,
            t0_onset=5,
            t0_offset=5,
            fs=1000.0,
            t_pre=2.0,
            t_post=2.0,
            step_ms=1.0,
            onset_baseline_window_ms=(-2.0, 0.0),
            retain_epochs=False,
            event_chunk_size=1,
        )
        self.assertIsNone(chunked["onset_epochs"])
        self.assertIsNone(chunked["offset_epochs"])
        np.testing.assert_allclose(chunked["onset_evoked"], out["onset_evoked"])
        np.testing.assert_allclose(chunked["offset_evoked"], out["offset_evoked"])

    def test_loader_rejects_mismatched_event_ids(self):
        with tempfile.TemporaryDirectory() as root:
            onset_folder = os.path.join(root, "onset")
            offset_folder = os.path.join(root, "offset")
            lfp = np.zeros((2, 2, 10))
            _write_leaf(onset_folder, lfp, event_ids=[1, 2])
            _write_leaf(offset_folder, lfp, event_ids=[2, 1])

            with self.assertRaisesRegex(ValueError, "event_ids"):
                load_paired_raw(onset_folder, offset_folder)

    def test_discovery_handles_alignment_in_middle_of_tree(self):
        with tempfile.TemporaryDirectory() as root:
            onset_folder = os.path.join(
                root, "behavior", "Pounce", "onset", "session1", "probe2"
            )
            offset_folder = os.path.join(
                root, "behavior", "Pounce", "offset", "session1", "probe2"
            )
            lfp = np.zeros((1, 2, 10))
            _write_leaf(onset_folder, lfp, event_ids=[7])
            _write_leaf(offset_folder, lfp, event_ids=[7])

            pairs = discover_behavior_aligned_pairs(root)

            self.assertEqual(len(pairs), 1)
            self.assertEqual(
                pairs[0]["name"], "behavior/Pounce/session1/probe2"
            )
            self.assertEqual(pairs[0]["onset_folder"], onset_folder)
            self.assertEqual(pairs[0]["offset_folder"], offset_folder)

    def test_runner_records_onset_baseline_metadata(self):
        with tempfile.TemporaryDirectory() as root:
            onset_folder = os.path.join(root, "play_bout", "onset", "session1")
            offset_folder = os.path.join(root, "play_bout", "offset", "session1")
            onset = np.zeros((2, 2, 10))
            onset[:, :, 3:5] = np.array([[[2]], [[4]]])
            offset = np.full((2, 2, 10), 20.0)
            _write_leaf(onset_folder, onset, event_ids=[1, 2])
            _write_leaf(offset_folder, offset, event_ids=[1, 2])

            result = run_behavior_aligned_kcsd(
                onset_folder,
                offset_folder,
                t_pre=2.0,
                t_post=2.0,
                onset_baseline_window_ms=(-2.0, 0.0),
                compute_kcsd=False,
                compute_fd=False,
                save=True,
                save_plots=False,
                retain_epochs=True,
            )

            self.assertEqual(result["meta"]["baseline_source"], "onset")
            self.assertEqual(
                result["meta"]["onset_baseline_window_ms"], [-2.0, 0.0]
            )
            self.assertEqual(result["meta"]["pairing"], "event_ids")
            np.testing.assert_allclose(result["offset_epochs"][0], 18.0)
            np.testing.assert_allclose(result["offset_epochs"][1], 16.0)
            for folder in (onset_folder, offset_folder):
                output = os.path.join(folder, "csd_1d_onset_baseline")
                self.assertTrue(
                    os.path.isfile(
                        os.path.join(output, "paired_baseline_metadata.json")
                    )
                )
                saved_baseline = np.load(
                    os.path.join(output, "onset_baseline.npy")
                )
                np.testing.assert_allclose(saved_baseline, result["onset_baseline"])

    def test_kcsd_progress_labels_include_behavior_dataset(self):
        with tempfile.TemporaryDirectory() as root:
            onset_folder = os.path.join(root, "behavior", "Pounce", "onset", "s1")
            offset_folder = os.path.join(root, "behavior", "Pounce", "offset", "s1")
            lfp = np.zeros((2, 2, 10), dtype=np.float32)
            _write_leaf(onset_folder, lfp, event_ids=[1, 2])
            _write_leaf(offset_folder, lfp, event_ids=[1, 2])
            np.save(os.path.join(onset_folder, "x_coords.npy"), np.zeros(2))
            np.save(os.path.join(offset_folder, "x_coords.npy"), np.zeros(2))

            fake_csd = {0.0: np.zeros((2, 4))}
            with patch(
                "behavior_aligned_kcsd.compute_kcsd_1d_columns",
                return_value=fake_csd,
            ) as compute:
                run_behavior_aligned_kcsd(
                    onset_folder,
                    offset_folder,
                    name="behavior/Pounce/s1",
                    t_pre=2.0,
                    t_post=2.0,
                    onset_baseline_window_ms=(-2.0, 0.0),
                    compute_kcsd=True,
                    compute_fd=False,
                    save=False,
                    save_plots=False,
                    plot=False,
                )

            labels = [
                call.kwargs["progress_label"] for call in compute.call_args_list
            ]
            self.assertEqual(
                labels,
                [
                    "behavior/Pounce/s1 | onset",
                    "behavior/Pounce/s1 | offset",
                ],
            )

    def test_column_average_interpolates_staggered_depths(self):
        x = np.array([0, 0, 1, 1], dtype=float)
        y = np.array([0, 20, 10, 30], dtype=float)
        evoked = np.array(
            [
                [0, 2],
                [20, 22],
                [110, 112],
                [130, 132],
            ],
            dtype=float,
        )
        csd = {
            0.0: np.array([[0, 2], [20, 22]], dtype=float),
            1.0: np.array([[110, 112], [130, 132]], dtype=float),
        }

        averaged = average_kcsd_columns(csd, x, y, evoked)

        np.testing.assert_allclose(averaged["y"], [10, 20])
        # At y=10: first column interpolates to 10, second is 110.
        # At y=20: first is 20, second interpolates to 120.
        np.testing.assert_allclose(
            averaged["csd_mean"], [[60, 62], [70, 72]]
        )
        np.testing.assert_allclose(
            averaged["voltage_mean"], [[60, 62], [70, 72]]
        )
        np.testing.assert_array_equal(
            averaged["n_columns_contributing"], [2, 2]
        )
        self.assertEqual(float(averaged["x_center"]), 0.5)

    def test_figure_organizer_creates_both_sort_orders_and_index(self):
        with tempfile.TemporaryDirectory() as root:
            onset_folder = os.path.join(
                root, "category", "play", "onset", "B1 D1", "probe2"
            )
            offset_folder = os.path.join(
                root, "category", "play", "offset", "B1 D1", "probe2"
            )
            lfp = np.zeros((1, 2, 10), dtype=np.float32)
            _write_leaf(onset_folder, lfp, event_ids=[1])
            _write_leaf(offset_folder, lfp, event_ids=[1])
            for folder in (onset_folder, offset_folder):
                figure_dir = os.path.join(
                    folder, "csd_1d_onset_baseline", "figures"
                )
                os.makedirs(figure_dir)
                with open(os.path.join(figure_dir, "kcsd_x8.0.png"), "wb") as stream:
                    stream.write(b"fake png")
                with open(os.path.join(figure_dir, "kcsd_xmean.png"), "wb") as stream:
                    stream.write(b"fake averaged png")

            result = organize_behavior_aligned_figures(root)

            expected_session_name = (
                "B1_D1__category__play__probe2__xmean.png"
            )
            expected_behavior_name = (
                "category__play__B1_D1__probe2__xmean.png"
            )
            self.assertTrue(
                os.path.isfile(
                    os.path.join(
                        result["output_root"],
                        "by_session",
                        "onset",
                        expected_session_name,
                    )
                )
            )
            self.assertTrue(
                os.path.isfile(
                    os.path.join(
                        result["output_root"],
                        "by_behavior",
                        "offset",
                        expected_behavior_name,
                    )
                )
            )
            self.assertEqual(result["n_source_figures"], 2)
            self.assertEqual(result["n_copies"], 4)
            with open(result["index_path"], newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["scope"], "category")
            self.assertEqual(rows[0]["behavior"], "play")
            self.assertEqual(rows[0]["probe"], "probe2")

    def test_saved_matrix_postprocessing_does_not_rerun_kcsd(self):
        with tempfile.TemporaryDirectory() as root:
            onset_folder = os.path.join(
                root, "behavior", "Pounce", "onset", "session1"
            )
            offset_folder = os.path.join(
                root, "behavior", "Pounce", "offset", "session1"
            )
            lfp = np.zeros((1, 4, 10), dtype=np.float32)
            _write_leaf(onset_folder, lfp, event_ids=[1])
            _write_leaf(offset_folder, lfp, event_ids=[1])

            for folder, shift in ((onset_folder, 0.0), (offset_folder, 5.0)):
                output = os.path.join(folder, "csd_1d_onset_baseline")
                figure_dir = os.path.join(output, "figures")
                os.makedirs(figure_dir)
                np.save(
                    os.path.join(output, "csd_x0.0.npy"),
                    np.array([[0, 2], [20, 22]], dtype=float) + shift,
                )
                np.save(
                    os.path.join(output, "csd_x1.0.npy"),
                    np.array([[110, 112], [130, 132]], dtype=float) + shift,
                )
                np.save(os.path.join(output, "y_x0.0.npy"), [0, 20])
                np.save(os.path.join(output, "y_x1.0.npy"), [10, 30])
                np.save(
                    os.path.join(output, "V_x0.0.npy"),
                    np.array([[0, 2], [20, 22]], dtype=float) + shift,
                )
                np.save(
                    os.path.join(output, "V_x1.0.npy"),
                    np.array([[110, 112], [130, 132]], dtype=float) + shift,
                )
                np.save(os.path.join(output, "time_voltage_x0.0.npy"), [-1, 0])
                np.save(os.path.join(output, "time_voltage_x1.0.npy"), [-1, 0])
                np.save(os.path.join(output, "time.npy"), [-1, 0])
                with open(
                    os.path.join(figure_dir, "kcsd_x0.0.png"), "wb"
                ) as stream:
                    stream.write(b"legacy")
                with open(
                    os.path.join(figure_dir, "kcsd_x1.0.png"), "wb"
                ) as stream:
                    stream.write(b"legacy")

            with patch(
                "behavior_aligned_kcsd.compute_kcsd_1d_columns"
            ) as compute:
                results = average_saved_behavior_kcsd_from_root(root)
            compute.assert_not_called()

            self.assertNotIn("error", results[0])
            for folder in (onset_folder, offset_folder):
                output = os.path.join(folder, "csd_1d_onset_baseline")
                np.testing.assert_allclose(
                    np.load(os.path.join(output, "y_mean_columns.npy")),
                    [10, 20],
                )
                self.assertEqual(
                    os.listdir(os.path.join(output, "figures")),
                    ["kcsd_xmean.png"],
                )


class AnimalPooledTests(unittest.TestCase):
    def test_trial_weighted_psth_aligns_different_depths_without_extrapolation(self):
        time_ms = np.array([0.0, 1.0])
        first = {
            "session": "day1",
            "n_events": 1,
            "x": np.array([0.0, 0.0, 0.0]),
            "y": np.array([0.0, 20.0, 40.0]),
            "onset_evoked": np.array([[0.0, 0.0], [10.0, 10.0], [20.0, 20.0]]),
            "offset_evoked": np.zeros((3, 2)),
            "time_ms": time_ms,
        }
        second = {
            "session": "day2",
            "n_events": 1,
            "x": np.array([0.0, 0.0, 0.0]),
            "y": np.array([10.0, 30.0, 50.0]),
            "onset_evoked": np.array([[10.0, 10.0], [30.0, 30.0], [50.0, 50.0]]),
            "offset_evoked": np.zeros((3, 2)),
            "time_ms": time_ms,
        }
        pooled = trial_weighted_psth([first, second])
        np.testing.assert_allclose(pooled["y"], [10.0, 20.0, 30.0, 40.0])
        np.testing.assert_allclose(
            pooled["onset_evoked"][:, 0], [7.5, 15.0, 22.5, 30.0]
        )
        self.assertEqual(pooled["aligned_by"], "shared_column_depths")
        self.assertEqual(pooled["n_events"], 2)

    def test_same_animal_trials_are_pooled_before_kcsd_for_categories_only(self):
        day1 = "B1A2 20240101 RatA"
        day2 = "B1A2 20240202 RatA"
        other = "B1B2 20240101 RatB"
        with tempfile.TemporaryDirectory() as root:
            self._write_constant_pair(
                root, "category", "play", day1, n_events=1, response=8.0
            )
            self._write_constant_pair(
                root, "category", "play", day2, n_events=3, response=0.0
            )
            self._write_constant_pair(
                root, "category", "play", other, n_events=1, response=100.0
            )
            self._write_constant_pair(
                root, "behavior", "Pounce", day1, n_events=1, response=50.0
            )
            self._write_constant_pair(
                root, "play_bout", "all", day1, n_events=2, response=4.0
            )
            self._write_constant_pair(
                root, "category", "aggression", day1, n_events=2, response=6.0
            )
            self._write_constant_pair(
                root, "category", "exploration", day1, n_events=2, response=3.0
            )
            self._write_constant_pair(
                root, "play_role", "self_play", day1, n_events=4, response=9.0
            )

            results = estimate_animal_pooled_kcsd(
                root,
                t_pre=2.0,
                t_post=2.0,
                onset_baseline_window_ms=(-2.0, 0.0),
                compute_kcsd=False,
                compute_fd=False,
                save=False,
                save_plots=False,
                plot=False,
            )

            self.assertTrue(all("error" not in result for result in results))
            play = np.load(
                os.path.join(
                    root,
                    "animal_pooled",
                    "category",
                    "play",
                    "onset",
                    "Batch1RatA2",
                    "lfp.npy",
                )
            )
            np.testing.assert_allclose(play[0, :, 2:], 2.0)
            self.assertFalse(
                os.path.exists(
                    os.path.join(root, "animal_pooled", "behavior", "Pounce")
                )
            )
            self.assertFalse(
                os.path.exists(
                    os.path.join(root, "animal_pooled", "play_role", "self_play")
                )
            )
            for behavior, response in (("aggression", 6.0), ("exploration", 3.0)):
                pooled_lfp = np.load(
                    os.path.join(
                        root,
                        "animal_pooled",
                        "category",
                        behavior,
                        "onset",
                        "Batch1RatA2",
                        "lfp.npy",
                    )
                )
                np.testing.assert_allclose(pooled_lfp[0, :, 2:], response)
            discovered = [
                pair["name"] for pair in discover_behavior_aligned_pairs(root)
            ]
            self.assertTrue(any(name.startswith("category/play/") for name in discovered))
            self.assertFalse(
                any(name.startswith("animal_pooled/") for name in discovered)
            )
            self.assertTrue(
                os.path.isfile(
                    os.path.join(
                        root,
                        "animal_pooled",
                        "play_bout",
                        "all",
                        "onset",
                        "Batch1RatA2",
                        "animal_pooling.json",
                    )
                )
            )
            other_lfp = np.load(
                os.path.join(
                    root,
                    "animal_pooled",
                    "category",
                    "play",
                    "onset",
                    "Batch1RatB2",
                    "lfp.npy",
                )
            )
            np.testing.assert_allclose(other_lfp[0, :, 2:], 100.0)

    def test_animal_pooled_figures_use_the_session_name_pattern(self):
        with tempfile.TemporaryDirectory() as root:
            leaves = [
                (
                    os.path.join(
                        root,
                        "animal_pooled",
                        "category",
                        "play",
                        "onset",
                        "Batch1RatA2",
                    ),
                    b"play-onset",
                ),
                (
                    os.path.join(
                        root,
                        "animal_pooled",
                        "category",
                        "play",
                        "offset",
                        "Batch1RatA2",
                        "probe2",
                    ),
                    b"play-offset",
                ),
                (
                    os.path.join(
                        root,
                        "animal_pooled",
                        "play_bout",
                        "all",
                        "onset",
                        "Batch1RatA2",
                    ),
                    b"bout",
                ),
            ]
            for folder, payload in leaves:
                figure_dir = os.path.join(folder, "csd_1d_onset_baseline", "figures")
                os.makedirs(figure_dir)
                with open(os.path.join(figure_dir, "kcsd_x0.0.png"), "wb") as stream:
                    stream.write(b"column")
                with open(os.path.join(figure_dir, "kcsd_xmean.png"), "wb") as stream:
                    stream.write(payload)
            outside = os.path.join(
                root, "category", "play", "onset", "B1A2 20240101 RatA",
                "csd_1d_onset_baseline", "figures",
            )
            os.makedirs(outside)
            with open(os.path.join(outside, "kcsd_xmean.png"), "wb") as stream:
                stream.write(b"session")

            result = organize_animal_pooled_figures(root)

            by_animal = os.path.join(
                result["output_root"],
                "by_animal",
                "onset",
                "Batch1RatA2__category__play__NPX1__xmean.png",
            )
            by_behavior = os.path.join(
                result["output_root"],
                "by_behavior",
                "offset",
                "category__play__Batch1RatA2__probe2__xmean.png",
            )
            bout = os.path.join(
                result["output_root"],
                "by_behavior",
                "onset",
                "play_bout__all__Batch1RatA2__NPX1__xmean.png",
            )
            self.assertTrue(os.path.isfile(by_animal))
            self.assertTrue(os.path.isfile(by_behavior))
            self.assertTrue(os.path.isfile(bout))
            with open(by_animal, "rb") as stream:
                self.assertEqual(stream.read(), b"play-onset")
            self.assertEqual(result["n_source_figures"], 3)
            self.assertEqual(result["n_copies"], 6)
            self.assertFalse(
                any("kcsd_x0.0" in name for name in os.listdir(os.path.dirname(by_animal)))
            )
            with open(result["index_path"], newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(
                {(row["scope"], row["behavior"], row["probe"]) for row in rows},
                {
                    ("category", "play", "NPX1"),
                    ("category", "play", "probe2"),
                    ("play_bout", "all", "NPX1"),
                },
            )

    def _write_constant_pair(self, root, scope, behavior, session, *, n_events, response):
        onset = np.zeros((n_events, 2, 10), dtype=float)
        onset[:, :, 5:7] = response
        offset = onset.copy()
        for alignment, lfp in (("onset", onset), ("offset", offset)):
            folder = os.path.join(root, scope, behavior, alignment, session)
            os.makedirs(folder, exist_ok=True)
            np.save(os.path.join(folder, "lfp.npy"), lfp)
            np.save(os.path.join(folder, "x_coords.npy"), np.zeros(2))
            np.save(os.path.join(folder, "y_coords.npy"), np.array([0.0, 20.0]))
            np.save(os.path.join(folder, "T0.npy"), np.array(5))
            np.save(os.path.join(folder, "fs.npy"), np.array(1000.0))
            np.save(
                os.path.join(folder, "event_ids.npy"),
                np.arange(1, n_events + 1),
            )


if __name__ == "__main__":
    unittest.main()
