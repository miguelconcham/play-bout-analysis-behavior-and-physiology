import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from extract_delta_phase_vs_lpag_distance import (  # noqa: E402
    align_geometry_to_lfp_rows,
    extract_lpag_ref_session_trace,
    extract_session_traces,
    lpag_center_map_index,
    lpag_ref_jobs,
    pli_wpli_phase_vs_ref,
    pooled_region_borders,
    probe_numbers,
    region_borders,
    interpolate_session_trend_outliers,
    session_animal_id,
    session_trend_outlier_mask,
)
from find_latest_merged_coherence_files import (  # noqa: E402
    find_latest_merged_coherence_files,
    is_lpag_ref_structure_name,
    is_merged_structure_name,
    partner_names_file,
)


def _linear_phase(distance_um, slope=0.002):
    return slope * np.asarray(distance_um, dtype=float)


class GeometryTests(unittest.TestCase):
    def test_animal_id_is_first_folder_token(self):
        self.assertEqual(
            session_animal_id("B1D1 1007 Dual ORIGINALLY PAG"),
            "B1D1",
        )
        self.assertEqual(session_animal_id("B4S2"), "B4S2")

    def test_npx1_is_one_probe_npx2_splits_by_x(self):
        npx1 = probe_numbers(np.array([11.0, 27.0, 43.0, 59.0]))
        np.testing.assert_array_equal(npx1, np.ones(4))
        npx2 = probe_numbers(np.array([8.0, 40.0, 258.0, 790.0]))
        np.testing.assert_array_equal(npx2, np.array([1.0, 1.0, 2.0, 4.0]))

    def test_chanmap_plus_one_geometry_uses_the_lfp_row_column(self):
        # Sites 1..4 have y = 0, 20, 40, 60. GENERATE stores them at
        # chanMap+1, so row 0 is empty and the phase matrix stays 4 x 4.
        y_by_lfp = np.array([0.0, 20.0, 40.0, 60.0])
        channel_map = np.full((5, 2), np.nan)
        channel_map[1:, 0] = 27.0
        channel_map[1:, 1] = y_by_lfp
        areas = ["", "SupCol", "DLPAG", "LPAG", "VLPAG"]
        phase = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                phase[i, j] = _linear_phase(y_by_lfp[i] - y_by_lfp[j])
        aligned_phase, aligned_map, aligned_areas = align_geometry_to_lfp_rows(
            phase, channel_map, areas
        )
        np.testing.assert_allclose(aligned_map[:, 1], y_by_lfp)
        self.assertEqual(aligned_areas[2], "LPAG")
        traces = extract_session_traces(
            aligned_phase, aligned_map, aligned_areas, "B1D1 chanMap plus one"
        )
        self.assertEqual(len(traces), 1)
        # LPAG is the only LPAG site, LFP row 3 (0-based index 2), y = 40.
        self.assertEqual(traces[0].lpag_center_index, 2)
        self.assertAlmostEqual(traces[0].phase_rad[2], 0.0)
        np.testing.assert_allclose(traces[0].phase_rad, _linear_phase(y_by_lfp - 40.0))

    def test_lpag_center_is_closest_lpag_site_to_mean(self):
        y = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0])
        x = np.full(y.shape, 27.0)
        areas = ["SupCol", "SupCol", "DLPAG", "LPAG", "LPAG", "VLPAG", "VLPAG", "DR"]
        # LPAG sites at 60 and 80 µm; mean y = 70, closest is 60 (first min).
        n = y.size
        phase = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                phase[i, j] = _linear_phase(y[i] - y[j])
        traces = extract_session_traces(
            phase, np.column_stack([x, y]), areas, "B1D1 session A"
        )
        self.assertEqual(len(traces), 1)
        trace = traces[0]
        self.assertEqual(trace.animal_id, "B1D1")
        self.assertEqual(trace.probe_num, 1)
        self.assertEqual(trace.lpag_center_index, 3)
        self.assertEqual(trace.lpag_center_y, 60.0)
        np.testing.assert_allclose(trace.distance_um, y - 60.0)
        np.testing.assert_allclose(trace.phase_rad, _linear_phase(y - 60.0))
        self.assertAlmostEqual(trace.phase_rad[3], 0.0)
        self.assertEqual(trace.area[3], "LPAG")
        np.testing.assert_array_equal(np.diff(trace.distance_um) >= 0, True)

    def test_resorts_scrambled_lfp_rows_by_channel_map_depth(self):
        # LFP row order is not anatomical; y jumps around the probe.
        y = np.array([80.0, 0.0, 40.0, 120.0])
        x = np.array([43.0, 11.0, 27.0, 59.0])
        areas = ["LPAG", "SupCol", "DLPAG", "VLPAG"]
        n = y.size
        phase = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                phase[i, j] = _linear_phase(y[i] - y[j])
        traces = extract_session_traces(
            phase, np.column_stack([x, y]), areas, "B1S3 scrambled"
        )
        self.assertEqual(len(traces), 1)
        trace = traces[0]
        np.testing.assert_array_equal(y[trace.channel_index], [0.0, 40.0, 80.0, 120.0])
        np.testing.assert_array_equal(trace.channel_index, [1, 2, 0, 3])
        self.assertTrue(np.all(np.diff(trace.distance_um) >= 0))
        np.testing.assert_allclose(trace.distance_um, y[trace.channel_index] - 80.0)
        np.testing.assert_allclose(
            trace.phase_rad, _linear_phase(y[trace.channel_index] - 80.0)
        )

    def test_region_borders_sit_between_area_changes(self):
        distance = np.array([-60.0, -40.0, -20.0, 0.0, 20.0, 40.0, 60.0, 80.0])
        areas = ["SupCol", "SupCol", "DLPAG", "LPAG", "LPAG", "VLPAG", "VLPAG", "DR"]
        border_d, border_from, border_to = region_borders(distance, areas)
        np.testing.assert_allclose(border_d, [-30.0, -10.0, 30.0, 70.0])
        self.assertEqual(border_from, ["SupCol", "DLPAG", "LPAG", "VLPAG"])
        self.assertEqual(border_to, ["DLPAG", "LPAG", "VLPAG", "DR"])

    def test_skips_sessions_without_lpag(self):
        y = np.array([0.0, 20.0, 40.0])
        x = np.full(y.shape, 27.0)
        areas = ["mPFC", "mPFC", "mPFC"]
        phase = np.eye(3)
        traces = extract_session_traces(
            phase, np.column_stack([x, y]), areas, "B1D1 mPFC"
        )
        self.assertEqual(traces, [])

    def test_npx2_keeps_only_probes_that_have_lpag(self):
        x = np.array([8.0, 8.0, 8.0, 258.0, 258.0, 258.0, 258.0])
        y = np.array([0.0, 40.0, 80.0, 0.0, 40.0, 80.0, 120.0])
        areas = ["SC", "SC", "DR", "DLPAG", "LPAG", "LPAG", "VLPAG"]
        n = y.size
        phase = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                phase[i, j] = 0.01 * (y[i] - y[j])
        traces = extract_session_traces(
            phase, np.column_stack([x, y]), areas, "B2S2 NPX2"
        )
        self.assertEqual(len(traces), 1)
        self.assertEqual(traces[0].probe_num, 2)
        self.assertEqual(traces[0].line_label, "B2S2 NPX2")
        self.assertTrue(np.all(x[traces[0].channel_index] == 258.0))
        # LPAG at y=40 and 80; mean 60; closest first min is y=40 (index 4).
        self.assertEqual(traces[0].lpag_center_index, 4)
        np.testing.assert_allclose(traces[0].distance_um, y[3:] - 40.0)

    def test_two_lpag_probes_become_two_labeled_lines(self):
        x = np.array([8.0, 8.0, 258.0, 258.0])
        y = np.array([10.0, 30.0, 10.0, 50.0])
        areas = ["LPAG", "VLPAG", "DLPAG", "LPAG"]
        phase = np.zeros((4, 4))
        traces = extract_session_traces(
            phase, np.column_stack([x, y]), areas, "B3D2 session"
        )
        self.assertEqual([t.probe_num for t in traces], [1, 2])
        self.assertEqual(traces[0].line_label, "B3D2 session probe 1")
        self.assertEqual(traces[1].line_label, "B3D2 session probe 2")

    def test_unwrap_follows_the_probe_not_the_branch_cut(self):
        y = np.linspace(0.0, 700.0, 8)
        x = np.full(y.shape, 27.0)
        areas = ["SupCol", "DLPAG", "LPAG", "LPAG", "VLPAG", "VLPAG", "DR", "DR"]
        n = y.size
        raw = np.linspace(-3.2, 3.2, n)
        phase = np.zeros((n, n))
        # Column 2 is the LPAG-center column; wrap it through ±pi.
        phase[:, 2] = np.angle(np.exp(1j * raw))
        traces = extract_session_traces(
            phase, np.column_stack([x, y]), areas, "B4S2"
        )
        self.assertEqual(len(traces), 1)
        self.assertGreater(np.ptp(traces[0].phase_unwrapped), np.pi)
        np.testing.assert_allclose(
            np.diff(traces[0].phase_unwrapped),
            np.diff(raw),
            atol=1e-10,
        )

    def test_pooled_borders_collapse_nearby_session_lines(self):
        x = np.full(4, 27.0)
        areas = ["DLPAG", "LPAG", "LPAG", "VLPAG"]
        phase = np.zeros((4, 4))
        a = extract_session_traces(
            phase,
            np.column_stack([x, np.array([0.0, 40.0, 80.0, 120.0])]),
            areas,
            "B1D1 s1",
        )[0]
        b = extract_session_traces(
            phase,
            np.column_stack([x, np.array([0.0, 38.0, 82.0, 120.0])]),
            areas,
            "B1D1 s2",
        )[0]
        distances, labels = pooled_region_borders([a, b], bin_um=10.0)
        self.assertEqual(list(labels), ["DLPAG | LPAG", "LPAG | VLPAG"])
        self.assertEqual(len(distances), 2)


class MergedFileSelectionTests(unittest.TestCase):
    def test_partner_names_follow_the_merged_stem(self):
        self.assertEqual(
            partner_names_file("coherence_structure_merged.mat"),
            "coherence_structure_animal_names_merged.mat",
        )
        self.assertEqual(
            partner_names_file("coherence_structure_all_together.mat"),
            "coherence_structure_animal_names_all_together.mat",
        )

    def test_ignores_unmerged_and_satellite_pieces(self):
        self.assertFalse(is_merged_structure_name("coherence_structure.mat"))
        self.assertFalse(
            is_merged_structure_name("coherence_structure_pli_struct_NPX2.mat")
        )
        self.assertFalse(
            is_merged_structure_name("coherence_structure_pli_struct_mpfc.mat")
        )
        self.assertFalse(
            is_merged_structure_name("coherence_structure_animal_names_merged.mat")
        )
        self.assertTrue(is_merged_structure_name("coherence_structure_merged.mat"))
        self.assertTrue(
            is_merged_structure_name("coherence_structure_all_together.mat")
        )

    def test_picks_the_newest_merged_save(self):
        entries = [
            ("coherence_structure.mat", 100.0),
            ("coherence_structure_pli_struct_NPX2.mat", 300.0),
            ("coherence_structure_merged.mat", 200.0),
            ("coherence_structure_animal_names_merged.mat", 200.0),
            ("coherence_structure_all_together.mat", 250.0),
            ("coherence_structure_animal_names_all_together.mat", 250.0),
        ]
        structure, names = find_latest_merged_coherence_files(entries)
        self.assertEqual(structure, "coherence_structure_all_together.mat")
        self.assertEqual(names, "coherence_structure_animal_names_all_together.mat")

    def test_2026_all_animal_save_is_combined_and_excludes_lpag_ref(self):
        self.assertTrue(is_merged_structure_name("coherence_pli_struct_2026.mat"))
        self.assertFalse(
            is_merged_structure_name("coherence_pli_struct_lpag_ref_2026.mat")
        )
        self.assertTrue(
            is_lpag_ref_structure_name("coherence_pli_struct_lpag_ref_2026.mat")
        )
        self.assertFalse(is_lpag_ref_structure_name("coherence_pli_struct_2026.mat"))
        self.assertEqual(
            partner_names_file("coherence_pli_struct_2026.mat"),
            "coherence_pli_structure_animal_names2026.mat",
        )
        self.assertEqual(
            partner_names_file("coherence_pli_struct_lpag_ref_2026.mat"),
            "coherence_pli_struct_lpag_ref_animal_names_2026.mat",
        )
        entries = [
            ("coherence_structure.mat", 100.0),
            ("coherence_structure_merged.mat", 200.0),
            ("coherence_pli_struct_lpag_ref_2026.mat", 400.0),
            ("coherence_pli_struct_lpag_ref_animal_names_2026.mat", 400.0),
            ("coherence_pli_struct_2026.mat", 300.0),
            ("coherence_pli_structure_animal_names2026.mat", 300.0),
        ]
        structure, names = find_latest_merged_coherence_files(entries)
        self.assertEqual(structure, "coherence_pli_struct_2026.mat")
        self.assertEqual(names, "coherence_pli_structure_animal_names2026.mat")
        structure, names = find_latest_merged_coherence_files(entries, kind="lpag_ref")
        self.assertEqual(structure, "coherence_pli_struct_lpag_ref_2026.mat")
        self.assertEqual(
            names, "coherence_pli_struct_lpag_ref_animal_names_2026.mat"
        )


class LpagRefPipelineTests(unittest.TestCase):
    def test_pairwise_extractor_skips_vector_phase(self):
        y = np.array([0.0, 40.0, 80.0, 120.0])
        x = np.full(y.shape, 27.0)
        areas = ["DLPAG", "LPAG", "LPAG", "VLPAG"]
        traces = extract_session_traces(
            _linear_phase(y - 40.0).reshape(-1, 1),
            np.column_stack([x, y]),
            areas,
            "B1D1 vector",
        )
        self.assertEqual(traces, [])

    def test_lpag_ref_extractor_sorts_scrambled_rows(self):
        y = np.array([80.0, 0.0, 40.0, 120.0])
        x = np.array([43.0, 11.0, 27.0, 59.0])
        areas = ["LPAG", "SupCol", "DLPAG", "VLPAG"]
        ref = 0
        phase = _linear_phase(y - y[ref])
        trace = extract_lpag_ref_session_trace(
            phase, np.column_stack([x, y]), areas, "B1S3 scrambled", ref
        )
        self.assertIsNotNone(trace)
        np.testing.assert_array_equal(y[trace.channel_index], [0.0, 40.0, 80.0, 120.0])
        np.testing.assert_array_equal(trace.channel_index, [1, 2, 0, 3])
        self.assertTrue(np.all(np.diff(trace.distance_um) >= 0))
        np.testing.assert_allclose(trace.distance_um, y[trace.channel_index] - 80.0)
        np.testing.assert_allclose(
            trace.phase_rad, _linear_phase(y[trace.channel_index] - 80.0)
        )
        self.assertAlmostEqual(trace.phase_rad[list(trace.channel_index).index(0)], 0.0)

    def test_npx2_lpag_jobs_keep_only_the_lpag_shank(self):
        x = np.array([8.0, 8.0, 8.0, 258.0, 258.0, 258.0, 258.0])
        y = np.array([0.0, 40.0, 80.0, 0.0, 40.0, 80.0, 120.0])
        chan_map = np.arange(x.size, dtype=float)
        jobs = lpag_ref_jobs(
            2, [(2, 30.0, 90.0)], x, y, chan_map, n_lfp=384
        )
        self.assertEqual(len(jobs), 1)
        job = jobs[0]
        self.assertEqual(job["probe_num"], 2)
        np.testing.assert_array_equal(job["x_coords"], [258.0, 258.0, 258.0, 258.0])
        # LPAG at y=40 and 80; mean 60; closest first min is y=40.
        self.assertEqual(job["ref_in_shank"], 1)
        self.assertEqual(job["detection_lfp_row"], 5.0)
        self.assertEqual(
            lpag_center_map_index(x, y, np.array([False, False, False, False, True, True, False])),
            4,
        )

    def test_phase_vs_ref_recovers_known_lag_and_is_linear_in_nch(self):
        n_time = 4000
        t = np.arange(n_time) / 2500.0
        omega = 2 * np.pi * 3.0
        lag = np.pi / 4
        analytic = np.vstack(
            [
                np.exp(1j * omega * t),
                np.exp(1j * (omega * t + lag)),
                np.exp(1j * (omega * t - lag)),
            ]
        )
        _pli, _wpli, phase_lag = pli_wpli_phase_vs_ref(analytic, 0)
        self.assertAlmostEqual(phase_lag[0], 0.0, places=6)
        self.assertAlmostEqual(phase_lag[1], lag, places=6)
        self.assertAlmostEqual(phase_lag[2], -lag, places=6)
        # Pairwise maps would be nCh x nCh; vs-ref is nCh.
        self.assertEqual(phase_lag.shape, (3,))


class SessionTrendOutlierTests(unittest.TestCase):
    def test_flags_a_point_far_from_the_session_trend(self):
        distance = np.linspace(-200.0, 200.0, 41)
        phase = 0.002 * distance
        phase[20] = 3.0
        outlier = session_trend_outlier_mask(distance, phase)
        self.assertTrue(outlier[20])
        self.assertEqual(int(np.count_nonzero(outlier)), 1)

    def test_keeps_a_smooth_session_trend(self):
        distance = np.linspace(-200.0, 200.0, 41)
        phase = 0.002 * distance + 0.01 * np.sin(distance / 40.0)
        outlier = session_trend_outlier_mask(distance, phase)
        self.assertFalse(np.any(outlier))

    def test_interpolates_the_outlier_onto_the_session_trend(self):
        distance = np.linspace(-200.0, 200.0, 41)
        phase = 0.002 * distance
        phase[20] = 3.0
        filled, outlier = interpolate_session_trend_outliers(distance, phase)
        self.assertTrue(outlier[20])
        self.assertAlmostEqual(filled[20], 0.002 * distance[20], places=6)
        self.assertAlmostEqual(filled[19], phase[19])

    def test_interpolates_when_two_channels_share_a_depth(self):
        distance = np.array(
            [-40.0, -20.0, -20.0, 0.0, 20.0, 20.0, 40.0, 60.0, 80.0, 100.0]
        )
        phase = 0.002 * distance
        phase[4] = 3.0
        filled, outlier = interpolate_session_trend_outliers(distance, phase)
        self.assertTrue(outlier[4])
        self.assertAlmostEqual(filled[4], 0.04, places=5)


if __name__ == "__main__":
    unittest.main()
