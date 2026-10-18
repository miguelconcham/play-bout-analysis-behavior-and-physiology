import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from area_csd_consistency import (
    animal_key_from_session,
    load_area_limits,
    resample_area_segment,
    select_area_for_depth,
    summarize_area_csd_across_sessions,
    AreaLimit,
)


def _row(area, probe_num, depth_start, depth_end, probe_area, animal="Batch1RatA2"):
    return AreaLimit(
        animal_name=animal,
        area=area,
        probe_num=probe_num,
        depth_start=depth_start,
        depth_end=depth_end,
        probe_area=probe_area,
    )


def _write_mean(leaf, value, y, time_ms):
    folder = os.path.join(leaf, "csd_1d_onset_baseline")
    os.makedirs(folder, exist_ok=True)
    csd = np.full((len(y), len(time_ms)), value, dtype=float)
    np.save(os.path.join(folder, "csd_mean_columns.npy"), csd)
    np.save(os.path.join(folder, "y_mean_columns.npy"), np.asarray(y, dtype=float))
    np.save(os.path.join(folder, "time_mean_columns.npy"), np.asarray(time_ms, dtype=float))


class AreaRuleTests(unittest.TestCase):
    def test_animal_key_matches_matlab_characters(self):
        self.assertEqual(
            animal_key_from_session("B1A2 2024-01-01 RatA notes"),
            "Batch1RatA2",
        )
        self.assertEqual(
            animal_key_from_session("B1XX 2024-01-01 Single2"),
            "Batch1Single2",
        )
        with self.assertRaises(ValueError):
            animal_key_from_session("OnlyTwo Tokens")
        with self.assertRaises(ValueError):
            animal_key_from_session("B1 2024-01-01 RatA")

    def test_npx_filters_and_narrowest_unique_span(self):
        rows = [
            _row("LPAG", 1, 0, 200, "PAG"),
            _row("VTA", 1, 90, 95, "VTA"),
            _row("DR", 2, 0, 80, "OTHER"),
            _row("Wide", 1, 0, 100, "PAG"),
            _row("Narrow", 1, 40, 60, "PAG"),
            _row("TieA", 3, 0, 50, "PAG"),
            _row("TieB", 3, 0, 50, "OTHER"),
        ]
        self.assertEqual(
            select_area_for_depth(50, rows, npx_type=1, probe_num=1).area,
            "Narrow",
        )
        # VTA is the narrowest row, but NPX1 only keeps Probe_Area PAG.
        self.assertEqual(
            select_area_for_depth(92, rows, npx_type=1, probe_num=1).area,
            "Wide",
        )
        self.assertEqual(
            select_area_for_depth(92, rows, npx_type=2, probe_num=1).area,
            "VTA",
        )
        self.assertEqual(
            select_area_for_depth(10, rows, npx_type=2, probe_num=2).area,
            "DR",
        )
        # NPX2 does not filter Probe_Area, so these equal spans are a tie.
        self.assertIsNone(select_area_for_depth(10, rows, npx_type=2, probe_num=3))
        # The MATLAB end is exclusive at depth_end + 1.
        self.assertIsNone(select_area_for_depth(201, rows, npx_type=1, probe_num=1))
        self.assertEqual(
            select_area_for_depth(200, rows, npx_type=1, probe_num=1).area,
            "LPAG",
        )

    def test_resample_does_not_extrapolate(self):
        y = np.array([0.0, 100.0])
        csd = np.array([[1.0, 1.0], [3.0, 3.0]])
        grid = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
        out = resample_area_segment(y, csd, 0.0, 200.0, grid)
        np.testing.assert_allclose(out[:3, 0], [1.0, 2.0, 3.0])
        self.assertTrue(np.all(np.isnan(out[3:])))

    def test_excel_roundtrip(self):
        import openpyxl

        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "Area_limits_GoodLooking.xlsx")
            workbook = openpyxl.Workbook()
            sheet = workbook.active
            sheet.append(
                [
                    "AnimalName",
                    "area",
                    "ProbeNum",
                    "depth_start",
                    "depth_end",
                    "Probe_Area",
                ]
            )
            sheet.append(["Batch1RatA2", "LPAG", 2, 0, 100, "PAG"])
            workbook.save(path)
            workbook.close()
            loaded = load_area_limits(path)
        self.assertEqual(len(loaded), 1)
        self.assertEqual(loaded[0].animal_name, "Batch1RatA2")
        self.assertEqual(loaded[0].probe_num, 2)
        self.assertEqual(loaded[0].area, "LPAG")


class AcrossSessionTests(unittest.TestCase):
    def test_probes_are_averaged_before_sessions_and_signs_are_scored(self):
        time_ms = np.array([-1000.0, 1000.0])
        y = np.array([0.0, 100.0])
        s1 = "B1A2 20240101 RatA"
        s2 = "B1A2 20240202 RatB"
        rows = []
        for animal in ("Batch1RatA2", "Batch1RatB2"):
            for probe_num in (1, 2):
                rows.append(
                    _row("LPAG", probe_num, 0, 100, "PAG", animal=animal)
                )
        with tempfile.TemporaryDirectory() as root:
            _write_mean(
                os.path.join(root, "category", "play", "onset", s1, "probe1"),
                2.0,
                y,
                time_ms,
            )
            _write_mean(
                os.path.join(root, "category", "play", "onset", s1, "probe2"),
                10.0,
                y,
                time_ms,
            )
            _write_mean(
                os.path.join(root, "category", "play", "onset", s2, "probe1"),
                18.0,
                y,
                time_ms,
            )
            _write_mean(
                os.path.join(root, "category", "aggression", "onset", s1, "probe1"),
                2.0,
                y,
                time_ms,
            )
            _write_mean(
                os.path.join(root, "category", "aggression", "onset", s1, "probe2"),
                10.0,
                y,
                time_ms,
            )
            _write_mean(
                os.path.join(root, "category", "aggression", "onset", s2, "probe1"),
                -6.0,
                y,
                time_ms,
            )
            _write_mean(
                os.path.join(root, "category", "play", "offset", s1, "probe1"),
                5.0,
                y,
                time_ms,
            )
            _write_mean(
                os.path.join(root, "behavior", "Pounce", "onset", s1),
                7.0,
                y,
                time_ms,
            )
            ghost = os.path.join(
                root,
                "category",
                "play",
                "onset",
                "B9Z9 20240101 Ghost",
                "probe1",
                "csd_1d_onset_baseline",
            )
            os.makedirs(ghost, exist_ok=True)
            np.save(os.path.join(ghost, "csd_mean_columns.npy"), np.ones((2, 2)))
            np.save(os.path.join(ghost, "y_mean_columns.npy"), y)
            np.save(os.path.join(ghost, "time_mean_columns.npy"), time_ms)
            missing = os.path.join(
                root,
                "category",
                "play",
                "onset",
                "B1A2 20240303 RatC",
                "probe1",
                "csd_1d_onset_baseline",
            )
            os.makedirs(missing, exist_ok=True)
            np.save(os.path.join(missing, "csd_x0.npy"), np.ones((2, 2)))

            result = summarize_area_csd_across_sessions(
                root,
                area_rows=rows,
                save_plots=True,
            )

            play = np.load(
                os.path.join(
                    root,
                    "area_consistency",
                    "category",
                    "play",
                    "onset",
                    "LPAG",
                    "mean_csd.npy",
                )
            )
            # Probe mean (2 and 10 -> 6) then session mean with 18 -> 12.
            # Treating probes as sessions would give (2 + 10 + 18) / 3 = 10.
            np.testing.assert_allclose(play, 12.0)
            sem = np.load(
                os.path.join(
                    root,
                    "area_consistency",
                    "category",
                    "play",
                    "onset",
                    "LPAG",
                    "sem_csd.npy",
                )
            )
            np.testing.assert_allclose(sem, 6.0)
            consistency = np.load(
                os.path.join(
                    root,
                    "area_consistency",
                    "category",
                    "play",
                    "onset",
                    "LPAG",
                    "sign_consistency.npy",
                )
            )
            np.testing.assert_allclose(consistency, 1.0)
            depth_counts = np.load(
                os.path.join(
                    root,
                    "area_consistency",
                    "category",
                    "play",
                    "onset",
                    "LPAG",
                    "n_sessions_by_depth.npy",
                )
            )
            np.testing.assert_array_equal(depth_counts, 2)
            with open(
                os.path.join(
                    root,
                    "area_consistency",
                    "category",
                    "play",
                    "onset",
                    "LPAG",
                    "sessions.txt",
                ),
                encoding="utf-8",
            ) as stream:
                sessions = stream.read().splitlines()
            self.assertEqual(sessions, [s1, s2])

            opposite = np.load(
                os.path.join(
                    root,
                    "area_consistency",
                    "category",
                    "aggression",
                    "onset",
                    "LPAG",
                    "mean_csd.npy",
                )
            )
            np.testing.assert_allclose(opposite, 0.0)
            opposite_consistency = np.load(
                os.path.join(
                    root,
                    "area_consistency",
                    "category",
                    "aggression",
                    "onset",
                    "LPAG",
                    "sign_consistency.npy",
                )
            )
            np.testing.assert_allclose(opposite_consistency, 0.5)

            offset = np.load(
                os.path.join(
                    root,
                    "area_consistency",
                    "category",
                    "play",
                    "offset",
                    "LPAG",
                    "mean_csd.npy",
                )
            )
            np.testing.assert_allclose(offset, 5.0)
            npx1 = np.load(
                os.path.join(
                    root,
                    "area_consistency",
                    "behavior",
                    "Pounce",
                    "onset",
                    "LPAG",
                    "mean_csd.npy",
                )
            )
            np.testing.assert_allclose(npx1, 7.0)
            self.assertTrue(
                os.path.isfile(
                    os.path.join(
                        root,
                        "area_consistency",
                        "category",
                        "play",
                        "onset",
                        "LPAG",
                        "mean_csd.png",
                    )
                )
            )
            self.assertTrue(
                os.path.isfile(
                    os.path.join(
                        root,
                        "area_consistency",
                        "category",
                        "play",
                        "onset",
                        "LPAG",
                        "sign_consistency.png",
                    )
                )
            )
            self.assertFalse(
                os.path.isfile(
                    os.path.join(
                        root,
                        "area_consistency",
                        "category",
                        "play",
                        "offset",
                        "LPAG",
                        "sign_consistency.png",
                    )
                )
            )
            self.assertEqual(result["color_limits"]["category/play/LPAG"], (-12.0, 12.0))
            self.assertEqual(result["n_missing_mean"], 1)
            self.assertTrue(
                any("Batch9Ghost9" in item["reason"] for item in result["skipped"])
            )
            play_summary = next(
                item
                for item in result["groups"]
                if item["scope"] == "category"
                and item["behavior"] == "play"
                and item["alignment"] == "onset"
            )
            self.assertEqual(play_summary["n_stable_bins"], play.size)
            self.assertEqual(play_summary["saved_kcsd_sign_at_max_abs"], "positive")
            opposite_summary = next(
                item
                for item in result["groups"]
                if item["behavior"] == "aggression"
            )
            self.assertEqual(opposite_summary["n_stable_bins"], 0)

    def test_module_does_not_estimate_kcsd(self):
        code_dir = os.path.dirname(os.path.abspath(__file__))
        with open(
            os.path.join(code_dir, "area_csd_consistency.py"), encoding="utf-8"
        ) as stream:
            source = stream.read()
        self.assertNotIn("compute_kcsd_1d_columns", source)
        self.assertNotIn("import kcsd", source)
        self.assertNotIn("csd_1d_pipeline", source)
        env = os.environ.copy()
        env["PYTHONPATH"] = code_dir + os.pathsep + env.get("PYTHONPATH", "")
        completed = subprocess.run(
            [
                sys.executable,
                "-c",
                "import area_csd_consistency, sys; "
                "assert 'kcsd' not in sys.modules; "
                "assert 'csd_1d_pipeline' not in sys.modules",
            ],
            cwd=code_dir,
            env=env,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)


if __name__ == "__main__":
    unittest.main()
