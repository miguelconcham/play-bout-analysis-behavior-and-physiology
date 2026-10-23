import os
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

import delta_animal_pooled_kcsd as delta_mod
from area_csd_consistency import AreaLimit
from csd_1d_pipeline import discover_delta_csd_datasets, plot_column
from delta_animal_pooled_kcsd import (
    estimate_animal_pooled_delta_kcsd,
    parse_delta_leaf,
    pool_delta_evoked,
    region_border_depths,
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


def _write_leaf(folder, response, n_events, y, *, t0=2, fs=1000.0, n_time=4):
    os.makedirs(folder, exist_ok=True)
    lfp = np.zeros((n_events, len(y), n_time), dtype=np.float32)
    lfp[:, :, t0:] = response
    np.save(os.path.join(folder, "lfp.npy"), lfp)
    np.save(os.path.join(folder, "x_coords.npy"), np.zeros(len(y)))
    np.save(os.path.join(folder, "y_coords.npy"), np.asarray(y, dtype=float))
    np.save(os.path.join(folder, "T0.npy"), np.array(t0))
    np.save(os.path.join(folder, "fs.npy"), np.array(fs))


def _write_area_table(path):
    import openpyxl

    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.append(
        ["AnimalName", "area", "ProbeNum", "depth_start", "depth_end", "Probe_Area"]
    )
    sheet.append(["Batch1RatA2", "LPAG", 1, 0, 40, "PAG"])
    sheet.append(["Batch1RatA2", "VTA", 1, 10, 30, "VTA"])
    sheet.append(["Batch1RatA2", "DR", 2, 20, 40, "OTHER"])
    sheet.append(["Batch1RatB2", "LPAG", 1, 0, 100, "PAG"])
    workbook.save(path)
    workbook.close()


class BorderTests(unittest.TestCase):
    def test_unique_edges_follow_the_probe_filter_and_the_plotted_range(self):
        rows = [
            _row("LPAG", 1, 0, 50, "PAG"),
            _row("VLPAG", 1, 50, 100, "PAG"),
            _row("VTA", 1, 10, 30, "VTA"),
            _row("DR", 2, 20, 80, "OTHER"),
            _row("Outside", 1, 200, 300, "PAG"),
        ]
        npx1 = region_border_depths(
            rows, y_min=10, y_max=90, npx_type=1, probe_num=1
        )
        np.testing.assert_allclose(npx1, [50])
        npx2 = region_border_depths(
            rows, y_min=10, y_max=90, npx_type=2, probe_num=2
        )
        np.testing.assert_allclose(npx2, [20, 80])

    def test_plot_draws_white_dashed_lines_at_the_borders(self):
        import matplotlib.axes

        captured = []
        real_axhline = matplotlib.axes.Axes.axhline

        def capture(self, y, *args, **kwargs):
            captured.append((float(y), kwargs.get("color"), kwargs.get("linestyle")))
            return real_axhline(self, y, *args, **kwargs)

        csd = np.array([[1.0, -1.0], [0.5, -0.5], [-1.0, 1.0]])
        time_ms = np.array([-2.0, 0.0, 2.0])
        y = np.array([0.0, 20.0, 40.0])
        with patch("matplotlib.axes.Axes.axhline", capture):
            plot_column(
                csd,
                time_ms,
                np.array([0, 2]),
                y,
                horizontal_lines=[0.0, 40.0],
                show=False,
            )
        self.assertEqual(captured, [(0.0, "white", "--"), (40.0, "white", "--")])


class PoolingTests(unittest.TestCase):
    def test_depth_alignment_matches_the_play_onset_overlap(self):
        time_ms = np.array([0.0, 1.0])
        pooled = pool_delta_evoked(
            [
                {
                    "session": "day1",
                    "n_events": 1,
                    "x": np.array([0.0, 0.0, 0.0]),
                    "y": np.array([0.0, 20.0, 40.0]),
                    "evoked": np.array([[0.0, 0.0], [10.0, 10.0], [20.0, 20.0]]),
                    "time_ms": time_ms,
                },
                {
                    "session": "day2",
                    "n_events": 1,
                    "x": np.array([0.0, 0.0, 0.0]),
                    "y": np.array([10.0, 30.0, 50.0]),
                    "evoked": np.array([[10.0, 10.0], [30.0, 30.0], [50.0, 50.0]]),
                    "time_ms": time_ms,
                },
            ]
        )
        np.testing.assert_allclose(pooled["y"], [10.0, 20.0, 30.0, 40.0])
        np.testing.assert_allclose(pooled["evoked"][:, 0], [7.5, 15.0, 22.5, 30.0])
        self.assertEqual(pooled["aligned_by"], "shared_column_depths")

    def test_peak_and_trough_are_pooled_separately_within_animal_and_probe(self):
        day1 = "B1A2 20240101 RatA"
        day2 = "B1A2 20240202 RatA"
        other = "B1B2 20240101 RatB"
        y = np.array([0.0, 20.0, 40.0])
        with tempfile.TemporaryDirectory() as root:
            _write_leaf(os.path.join(root, "peak", day1), 8.0, 1, y)
            _write_leaf(os.path.join(root, "peak", day2), 0.0, 3, y)
            _write_leaf(os.path.join(root, "peak", other), 100.0, 1, y)
            _write_leaf(os.path.join(root, "trough", day1), 5.0, 2, y)
            _write_leaf(os.path.join(root, "peak", day1, "probe2"), 7.0, 4, y)
            _write_leaf(
                os.path.join(root, "animal_pooled", "peak", "Batch1RatA2"),
                99.0,
                1,
                y,
            )
            area_table = os.path.join(root, "Area_limits_GoodLooking.xlsx")
            _write_area_table(area_table)

            discovered = discover_delta_csd_datasets(root)
            self.assertTrue(all("animal_pooled" not in item["folder"] for item in discovered))
            self.assertIsNone(
                parse_delta_leaf(
                    root, os.path.join(root, "animal_pooled", "peak", "Batch1RatA2")
                )
            )

            results = estimate_animal_pooled_delta_kcsd(
                root,
                area_table,
                t_pre=2.0,
                t_post=2.0,
                fs=1000.0,
                compute_kcsd=False,
                save_plots=False,
                plot=False,
            )
            self.assertTrue(all("error" not in result for result in results))

            peak = np.load(
                os.path.join(root, "animal_pooled", "peak", "Batch1RatA2", "lfp.npy")
            )
            np.testing.assert_allclose(peak[0, :, 2:], 2.0)
            trough = np.load(
                os.path.join(root, "animal_pooled", "trough", "Batch1RatA2", "lfp.npy")
            )
            np.testing.assert_allclose(trough[0, :, 2:], 5.0)
            other_peak = np.load(
                os.path.join(root, "animal_pooled", "peak", "Batch1RatB2", "lfp.npy")
            )
            np.testing.assert_allclose(other_peak[0, :, 2:], 100.0)
            probe = np.load(
                os.path.join(
                    root, "animal_pooled", "peak", "Batch1RatA2", "probe2", "lfp.npy"
                )
            )
            np.testing.assert_allclose(probe[0, :, 2:], 7.0)

    def test_old_session_peak_layout_is_parsed(self):
        session = "B1A2 20240101 RatA"
        with tempfile.TemporaryDirectory() as root:
            folder = os.path.join(root, session, "peak")
            _write_leaf(folder, 4.0, 1, np.array([0.0, 20.0]))
            parsed = parse_delta_leaf(root, folder)
        self.assertEqual(parsed["alignment"], "peak")
        self.assertEqual(parsed["session"], session)
        self.assertEqual(parsed["probe"], "")

    def test_mean_figure_receives_borders_for_the_plotted_probe(self):
        day1 = "B1A2 20240101 RatA"
        y = np.array([0.0, 20.0, 40.0])
        captured = {}

        def fake_run(cfg):
            folder = cfg["folder"]
            evoked = np.load(os.path.join(folder, "lfp.npy"))[0]
            x = np.load(os.path.join(folder, "x_coords.npy"))
            depths = np.load(os.path.join(folder, "y_coords.npy"))
            time_ms = np.load(os.path.join(folder, "time_ms.npy"))
            return {
                "csd_kcsd": {0.0: np.ones((len(depths), 2))},
                "x": x,
                "y": depths,
                "evoked": evoked,
                "time_ms": time_ms,
                "time_indices": np.array([0, 1]),
            }

        def fake_plot(*_args, **kwargs):
            captured.setdefault("calls", []).append(kwargs.get("horizontal_lines"))

        with tempfile.TemporaryDirectory() as root:
            _write_leaf(os.path.join(root, "peak", day1), 8.0, 1, y)
            _write_leaf(os.path.join(root, "peak", day1, "probe2"), 7.0, 1, y)
            area_table = os.path.join(root, "Area_limits_GoodLooking.xlsx")
            _write_area_table(area_table)
            with patch.object(delta_mod, "run_dataset", fake_run), patch.object(
                delta_mod, "plot_column", fake_plot
            ):
                results = estimate_animal_pooled_delta_kcsd(
                    root,
                    area_table,
                    t_pre=2.0,
                    t_post=2.0,
                    fs=1000.0,
                    compute_kcsd=True,
                    save_plots=True,
                    plot=False,
                )
        self.assertTrue(all("error" not in result for result in results))
        by_name = {result["name"]: result["region_borders_um"] for result in results}
        self.assertEqual(by_name["animal_pooled/peak/Batch1RatA2"], [0.0, 40.0])
        self.assertEqual(by_name["animal_pooled/peak/Batch1RatA2/probe2"], [20.0, 40.0])
        self.assertIn([0.0, 40.0], captured["calls"])
        self.assertIn([20.0, 40.0], captured["calls"])


if __name__ == "__main__":
    unittest.main()
