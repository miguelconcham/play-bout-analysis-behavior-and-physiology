import os
import tempfile
import unittest

import numpy as np

from behavior_aligned_kcsd import (
    discover_behavior_aligned_pairs,
    epoch_paired_with_onset_baseline,
    load_paired_raw,
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


if __name__ == "__main__":
    unittest.main()
