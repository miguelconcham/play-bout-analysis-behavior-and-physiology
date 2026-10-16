import csv
import os
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from behavior_aligned_kcsd import (
    average_kcsd_columns,
    discover_behavior_aligned_pairs,
    epoch_paired_with_onset_baseline,
    load_paired_raw,
    organize_behavior_aligned_figures,
    run_behavior_aligned_kcsd,
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


if __name__ == "__main__":
    unittest.main()
