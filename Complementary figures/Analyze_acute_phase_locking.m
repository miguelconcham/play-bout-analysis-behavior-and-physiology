%% Analyze_acute_phase_locking
% Load pooled GENERATE_ACUTE_PHASE_LOCKING results. All figures use neurons
% in focus_regions only (pooled; default PAG columns, SC, Raphe).
% Test whether delta trough-locked cells are more often breathing-locked
% than a random (non-trough) neuron. Fisher's exact is the primary test;
% chi-square is reported alongside. Breathing peak vs trough PSTHs are
% shown both by LFP lock and by breathing preferred phase.

%% Paths and parameters
data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
run('\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

saving_folder = fullfile(data_root, 'Analysis results', 'acute phase locking');
freq_range    = [1 5];
freq_tag      = sprintf('%g-%gHz', freq_range(1), freq_range(2));

alpha         = 0.01;          % PPC p-value for "locked"
freq2filter   = 8;             % PSTH smoothing (Hz), as in Fig 3
c_lim         = [-20 20];      % % change color scale
y_lim         = [-15 15];      % mean PSTH y-axis (% change)
angle_bins    = -pi:(pi/32):pi;
breath_raster_n = 1;           % row in the printed breathing-locked table; [] skip; NaN prompt
breath_mod_win  = [-0.4 0.4];  % window (s) for peak–trough around the event
% All plots use neurons in these regions only (pooled; no per-area figures).
% 'PAG' = every PAG column. To plot one area later, set e.g. {'DlPAG'} or {'SupCol'}.
focus_regions = {'PAG', 'SupCol', 'Raphe'};
mod_outlier_k = 1.5;           % Tukey IQR multiplier; drop outliers on delta vs breath % scatter

load(fullfile(saving_folder, ['acute_phase_locking_', freq_tag, '.mat']), 'acute_struct');
if exist(fullfile(saving_folder, 'acute_phase_locking_session_names.mat'), 'file')
    load(fullfile(saving_folder, 'acute_phase_locking_session_names.mat'), 'session_names');
else
    session_names = {acute_struct.session_folder}';
end

%% Concatenate neurons
phase_stat_names = acute_struct(1).phase_stat_names;
psth_time        = acute_struct(1).psth_time;
psth_bins        = mean(diff(psth_time));
smooth_n         = max(1, round((1 / freq2filter) / psth_bins));
x_lim            = [psth_time(1) psth_time(end)];

lfp_stats        = [];
breath_stats     = [];
psth_peak        = [];
psth_trough      = [];
breath_psth_peak   = [];
breath_psth_trough = [];
cluster_ids        = [];
channel_list_all   = [];
sess_idx_all       = [];
session_label      = {};
date_list          = [];
exp_list           = [];
for ns = 1:numel(acute_struct)
    if ~isequal(acute_struct(ns).psth_time, psth_time)
        warning('Session %s has a different PSTH time base; skipped.', acute_struct(ns).session_folder);
        continue
    end
    n = size(acute_struct(ns).lfp_phase_stats, 1);
    lfp_stats        = [lfp_stats; acute_struct(ns).lfp_phase_stats]; %#ok<AGROW>
    breath_stats     = [breath_stats; acute_struct(ns).breathing_phase_stats]; %#ok<AGROW>
    psth_peak        = [psth_peak; acute_struct(ns).lfp_psth_peak]; %#ok<AGROW>
    psth_trough      = [psth_trough; acute_struct(ns).lfp_psth_trough]; %#ok<AGROW>
    breath_psth_peak   = [breath_psth_peak; acute_struct(ns).breathing_psth_peak]; %#ok<AGROW>
    breath_psth_trough = [breath_psth_trough; acute_struct(ns).breathing_psth_trough]; %#ok<AGROW>
    cluster_ids      = [cluster_ids; double(acute_struct(ns).clusters_list(:))]; %#ok<AGROW>
    ch_here = nan(n, 1);
    if isfield(acute_struct, 'channel_list') && ~isempty(acute_struct(ns).channel_list)
        ch_here = double(acute_struct(ns).channel_list(:));
    end
    channel_list_all = [channel_list_all; ch_here]; %#ok<AGROW>
    sess_idx_all     = [sess_idx_all; repmat(ns, n, 1)]; %#ok<AGROW>
    session_label    = [session_label; repmat({acute_struct(ns).session_folder}, n, 1)]; %#ok<AGROW>
    date_list        = [date_list; repmat(acute_struct(ns).date, n, 1)]; %#ok<AGROW>
    exp_list         = [exp_list; repmat(acute_struct(ns).experiment, n, 1)]; %#ok<AGROW>
end

i_ang  = find(strcmp(phase_stat_names, 'PreferedAngle'), 1);
i_ppc  = find(strcmp(phase_stat_names, 'PPC'), 1);
i_ppcv = find(strcmp(phase_stat_names, 'PPCPval'), 1);

lfp_ang    = lfp_stats(:, i_ang);
lfp_ppc    = lfp_stats(:, i_ppc);
lfp_ppcv   = lfp_stats(:, i_ppcv);
br_ppc     = breath_stats(:, i_ppc);
br_ppcv    = breath_stats(:, i_ppcv);
br_ang     = breath_stats(:, i_ang);

valid        = ~isnan(lfp_ppc) & ~isnan(lfp_ppcv);
delta_lock   = valid & lfp_ppcv < alpha;
trough_lock  = delta_lock & (lfp_ang > pi/2 | lfp_ang < -pi/2);
peak_lock    = delta_lock & ~trough_lock;
unlocked     = valid & ~delta_lock;
breath_lock  = valid & ~isnan(br_ppcv) & br_ppcv < alpha;

disp(['Neurons with LFP stats: ', num2str(sum(valid)), ' / ', num2str(numel(valid))])
disp(['Delta locked (PPC p<', num2str(alpha), '): ', num2str(sum(delta_lock))])
disp(['  trough: ', num2str(sum(trough_lock)), '   peak: ', num2str(sum(peak_lock))])
disp(['Breathing locked: ', num2str(sum(breath_lock))])

%% Match NeuronSummary (DATE + EXP + ID)
% Restrict every plot to focus_regions (pooled). Edit that list at the top
% to keep a single area.
brain_region = strings(numel(cluster_ids), 1);
brain_area   = strings(numel(cluster_ids), 1);
in_focus     = false(numel(cluster_ids), 1);
ns_file = fullfile(data_root, 'Acute data', 'population analysis', 'NeuronSummary.xlsx');
if ~exist(ns_file, 'file')
    warning('Analyze_acute_phase_locking:NeuronSummary', ...
        'No NeuronSummary.xlsx; cannot apply focus_regions, keeping all neurons.')
    in_focus(:) = true;
else
    NS = readtable(ns_file);
    n_matched = 0;
    for i = 1:numel(cluster_ids)
        hit = NS.DATE == date_list(i) & NS.EXP == exp_list(i) & NS.ID == cluster_ids(i);
        if ~any(hit)
            brain_region(i) = "unmatched";
            brain_area(i) = "unmatched";
            continue
        end
        rows = find(hit);
        brain_region(i) = strtrim(string(NS.BRAIN_REGION(rows(1))));
        brain_area(i)   = strtrim(string(NS.AREA(rows(1))));
        n_matched = n_matched + 1;
    end
    brain_region(brain_region == "Ohters") = "Others";
    brain_region(brain_region == "Rhape")  = "Raphe";
    brain_region(ismember(lower(brain_region), ["vlpag"])) = "VlPAG";
    brain_region(ismember(lower(brain_region), ["dlpag"])) = "DlPAG";
    brain_region(ismember(lower(brain_region), ["sc", "supcol", "superior colliculus"])) = "SupCol";
    brain_region(brain_region == "") = "unlabeled";
    disp(['NeuronSummary matched: ', num2str(n_matched), ' / ', num2str(numel(cluster_ids))])
    disp('BRAIN_REGION counts among matched cells:')
    [g, names] = findgroups(brain_region);
    disp(table(names, splitapply(@numel, g, g), 'VariableNames', {'region', 'n'}))
    in_focus = neurons_in_focus_list(brain_region, brain_area, focus_regions);
end
focus_label = strjoin(focus_regions, ', ');
disp(['Focus regions (pooled): ', focus_label])
disp(['  in list: ', num2str(sum(in_focus)), ' / ', num2str(numel(in_focus))])
valid      = valid & in_focus;
delta_lock = valid & lfp_ppcv < alpha;
trough_lock = delta_lock & (lfp_ang > pi/2 | lfp_ang < -pi/2);
peak_lock   = delta_lock & ~trough_lock;
unlocked    = valid & ~delta_lock;
breath_lock = valid & ~isnan(br_ppcv) & br_ppcv < alpha;
disp(['After region filter, LFP stats: ', num2str(sum(valid))])
disp(['  delta locked: ', num2str(sum(delta_lock)), ...
    '  trough: ', num2str(sum(trough_lock)), '  peak: ', num2str(sum(peak_lock))])
disp(['  breathing locked: ', num2str(sum(breath_lock))])

%% Preferred-angle distribution (locked cells)
figure('Name', 'Acute delta preferred phase')
subplot(1, 2, 1)
polarhistogram(lfp_ang(delta_lock), angle_bins, 'Normalization', 'count', ...
    'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none')
title({'Preferred phase, delta-locked', [num2str(sum(delta_lock)), ' cells']})

subplot(1, 2, 2)
hold on
histogram(180 * lfp_ang(trough_lock) / pi, -180:11.25:180, 'Normalization', 'count', ...
    'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.6)
histogram(180 * lfp_ang(peak_lock) / pi, -180:11.25:180, 'Normalization', 'count', ...
    'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.6)
xline(-90, ':k'); xline(90, ':k'); xline(0, ':k')
xlabel('Preferred angle (deg)')
ylabel('Count')
legend({'Trough-locked', 'Peak-locked'})
title('Peak vs trough (delta-locked)')
xlim([-180 180])

%% (Fig 3 style) Heatmaps / mean PSTHs of locked vs unlocked — peak-aligned
figure('Name', 'Acute delta PSTH peak-aligned', 'units', 'normalized', 'outerposition', [0.05 0.05 0.4 0.9])
plot_entrainment_psth_stack(psth_time, psth_peak, delta_lock, unlocked, lfp_ang, ...
    smooth_n, x_lim, c_lim, y_lim, 'Peak-aligned')

%% Same layout, trough-aligned
figure('Name', 'Acute delta PSTH trough-aligned', 'units', 'normalized', 'outerposition', [0.45 0.05 0.4 0.9])
plot_entrainment_psth_stack(psth_time, psth_trough, delta_lock, unlocked, lfp_ang, ...
    smooth_n, x_lim, c_lim, y_lim, 'Trough-aligned')

%% Breathing-lock proportions: all neurons vs trough-locked
% Primary comparison: trough-locked vs everyone else (a "random" neuron).
% 2x2 counts; Fisher's exact is preferred with modest n (3 acute sessions).
% Chi-square is shown as well. Also report trough vs peak among delta-locked.
breathing_test = report_breathing_lock_proportions( ...
    valid, trough_lock, peak_lock, delta_lock, breath_lock, focus_label);

%% Breathing peak- and trough-aligned PSTHs
% GENERATE peak vs trough PSTHs are complementary when cells are split by
% *breathing* preferred phase (phase 0 = peak, +/-pi = trough). Grouping by
% LFP delta lock mixes those phases (LFP-trough cells are ~half breath-peak
% and ~half breath-trough), so those two alignments look similar — that is
% not a swapped-event bug. Plot both groupings.

group_col  = {[0.45 0.45 0.45], [0.15 0.35 0.85], [0.85 0.2 0.2]};

lfp_group_mask = {unlocked, trough_lock, peak_lock};
lfp_group_name = {'LFP unlocked', 'LFP trough-lock', 'LFP peak-lock'};

br_trough_pref = breath_lock & (br_ang > pi/2 | br_ang < -pi/2);
br_peak_pref   = breath_lock & ~br_trough_pref;
br_unlocked    = valid & ~breath_lock;
br_group_mask  = {br_unlocked, br_trough_pref, br_peak_pref};
br_group_name  = {'Breath unlocked', 'Breath trough-pref', 'Breath peak-pref'};

disp('--- Breathing phase mix inside LFP groups ---')
disp(sprintf('  LFP-trough cells: n=%d  breath-locked=%d  trough-pref=%d  peak-pref=%d', ...
    sum(trough_lock), sum(trough_lock & breath_lock), ...
    sum(trough_lock & br_trough_pref), sum(trough_lock & br_peak_pref)))
disp(sprintf('  LFP-peak cells:   n=%d  breath-locked=%d  trough-pref=%d  peak-pref=%d', ...
    sum(peak_lock), sum(peak_lock & breath_lock), ...
    sum(peak_lock & br_trough_pref), sum(peak_lock & br_peak_pref)))

breathing_psth_test = struct();
breathing_psth_test.by_lfp_lock = struct( ...
    'peak', compare_breathing_psth_groups(psth_time, breath_psth_peak, ...
        lfp_group_mask, lfp_group_name, group_col, smooth_n, x_lim, 'peak, LFP groups'), ...
    'trough', compare_breathing_psth_groups(psth_time, breath_psth_trough, ...
        lfp_group_mask, lfp_group_name, group_col, smooth_n, x_lim, 'trough, LFP groups'));
breathing_psth_test.by_breath_phase = struct( ...
    'peak', compare_breathing_psth_groups(psth_time, breath_psth_peak, ...
        br_group_mask, br_group_name, group_col, smooth_n, x_lim, 'peak, breath-phase groups'), ...
    'trough', compare_breathing_psth_groups(psth_time, breath_psth_trough, ...
        br_group_mask, br_group_name, group_col, smooth_n, x_lim, 'trough, breath-phase groups'));

plot_breathing_peak_vs_trough(psth_time, breath_psth_peak, breath_psth_trough, ...
    br_group_mask, br_group_name, group_col, smooth_n, x_lim);

%% LFP × breathing lock combinations
% 2×2: LFP trough/peak × breath trough/peak. Gray = LFP-unlocked cells with
% the same breathing lock as that panel (trough-pref or peak-pref), so the
% comparison is breath-lock-matched. Peak- and trough-aligned figures share
% the same y-axis.
y_lim_combo = plot_lfp_breath_lock_combos(psth_time, breath_psth_peak, breath_psth_trough, ...
    unlocked, trough_lock, peak_lock, br_trough_pref, br_peak_pref, ...
    smooth_n, x_lim, [], focus_label, true); %#ok<NASGU>

%% Breathing and delta modulation as % of mean rate
% Amplitude = 100 * (max-min) / (2 * mean PSTH rate) inside breath_mod_win.
breath_pct_at_0 = nan(numel(valid), 1);
breath_pct_amp  = nan(numel(valid), 1);
breath_mean_hz  = nan(numel(valid), 1);
breath_pct_ptp  = nan(numel(valid), 1);
delta_pct_at_0 = nan(numel(valid), 1);
delta_pct_amp  = nan(numel(valid), 1);
delta_mean_hz  = nan(numel(valid), 1);
delta_pct_ptp  = nan(numel(valid), 1);
for i = 1:numel(valid)
    if br_trough_pref(i)
        s_br = psth_mod_of_mean(breath_psth_trough(i, :), psth_time, smooth_n, breath_mod_win);
    elseif br_peak_pref(i)
        s_br = psth_mod_of_mean(breath_psth_peak(i, :), psth_time, smooth_n, breath_mod_win);
    else
        s_tr = psth_mod_of_mean(breath_psth_trough(i, :), psth_time, smooth_n, breath_mod_win);
        s_pk = psth_mod_of_mean(breath_psth_peak(i, :), psth_time, smooth_n, breath_mod_win);
        if isfinite(s_pk.pct_amp) && (~isfinite(s_tr.pct_amp) || s_pk.pct_amp > s_tr.pct_amp)
            s_br = s_pk;
        else
            s_br = s_tr;
        end
    end
    breath_pct_at_0(i) = s_br.pct_at_0;
    breath_pct_amp(i)  = s_br.pct_amp;
    breath_mean_hz(i)  = s_br.mean_hz;
    breath_pct_ptp(i)  = s_br.pct_ptp;

    if trough_lock(i)
        s_d = psth_mod_of_mean(psth_trough(i, :), psth_time, smooth_n, breath_mod_win);
    elseif peak_lock(i)
        s_d = psth_mod_of_mean(psth_peak(i, :), psth_time, smooth_n, breath_mod_win);
    else
        s_tr = psth_mod_of_mean(psth_trough(i, :), psth_time, smooth_n, breath_mod_win);
        s_pk = psth_mod_of_mean(psth_peak(i, :), psth_time, smooth_n, breath_mod_win);
        if isfinite(s_pk.pct_amp) && (~isfinite(s_tr.pct_amp) || s_pk.pct_amp > s_tr.pct_amp)
            s_d = s_pk;
        else
            s_d = s_tr;
        end
    end
    delta_pct_at_0(i) = s_d.pct_at_0;
    delta_pct_amp(i)  = s_d.pct_amp;
    delta_mean_hz(i)  = s_d.mean_hz;
    delta_pct_ptp(i)  = s_d.pct_ptp;
end

lfp_mod_groups = {trough_lock, peak_lock, unlocked};
lfp_mod_names  = {'LFP trough-lock', 'LFP peak-lock', 'LFP unlocked'};
lfp_mod_col    = {[0.15 0.35 0.85], [0.85 0.2 0.2], [0.55 0.55 0.55]};
disp('--- Breathing rate modulation by LFP class ---')
for g = 1:numel(lfp_mod_groups)
    idx = lfp_mod_groups{g} & isfinite(breath_pct_amp);
    if ~any(idx)
        disp(sprintf('  %s: n=0', lfp_mod_names{g}))
        continue
    end
    disp(sprintf(['  %s n=%d  mean rate %.1f Hz  |  ' ...
        'modulation +/-%.1f%% of mean (IQR %.1f-%.1f)'], ...
        lfp_mod_names{g}, sum(idx), median(breath_mean_hz(idx), 'omitnan'), ...
        median(breath_pct_amp(idx), 'omitnan'), ...
        prctile(breath_pct_amp(idx), 25), prctile(breath_pct_amp(idx), 75)))
end

lfp_lock_str = strings(numel(valid), 1);
lfp_lock_str(:) = "unlocked";
lfp_lock_str(trough_lock) = "trough-locked";
lfp_lock_str(peak_lock) = "peak-locked";
breath_pref_str = strings(numel(valid), 1);
breath_pref_str(:) = "breath-unlocked";
breath_pref_str(br_trough_pref) = "breath-trough";
breath_pref_str(br_peak_pref) = "breath-peak";
breathing_rate_mod = table(session_label, cluster_ids, lfp_lock_str, breath_pref_str, ...
    breath_mean_hz, breath_pct_at_0, breath_pct_amp, breath_pct_ptp, ...
    delta_mean_hz, delta_pct_at_0, delta_pct_amp, delta_pct_ptp, ...
    'VariableNames', {'session', 'id', 'LFP_lock', 'breath_pref', ...
    'mean_hz', 'pct_at_event', 'pct_amp_of_mean', 'pct_peak_to_trough', ...
    'delta_mean_hz', 'delta_pct_at_event', 'delta_pct_amp', 'delta_pct_peak_to_trough'});

edges_amp = 0:2:80;
figure('Name', 'Breathing modulation histograms by LFP class')
yl_h = 0;
axs = gobjects(1, 3);
for g = 1:3
    axs(g) = subplot(1, 3, g);
    idx = lfp_mod_groups{g} & isfinite(breath_pct_amp);
    histogram(breath_pct_amp(idx), edges_amp, 'Normalization', 'probability', ...
        'FaceColor', lfp_mod_col{g}, 'FaceAlpha', 0.85, 'EdgeColor', 'none')
    xlabel('% of mean rate')
    if g == 1
        ylabel('Probability')
    end
    med = median(breath_pct_amp(idx), 'omitnan');
    xline(med, ':k')
    title(sprintf('%s  n=%d  median \\pm%.1f%%', lfp_mod_names{g}, sum(idx), med))
    yl_h = max(yl_h, axs(g).YLim(2));
end
linkaxes(axs, 'xy')
ylim([0 max(yl_h, eps)])
sgtitle(sprintf('Breathing modulation (%% of mean)  %s', focus_label))

%% Delta vs breathing modulation and PPC
both_keep = valid & isfinite(breath_pct_amp) & isfinite(delta_pct_amp) ...
    & isfinite(lfp_ppc) & isfinite(br_ppc);
disp(sprintf('%s: %d cells with both modulations.', focus_label, sum(both_keep)))
if ~any(valid)
    warning('Analyze_acute_phase_locking:focusMod', ...
        'No cells in focus_regions; skip delta vs breathing figure.')
else
    figure('Name', ['Delta vs breathing modulation  ', focus_label], ...
        'units', 'normalized', 'outerposition', [0.08 0.08 0.84 0.82])
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact')

    nexttile
    hold on
    leg = {};
    for g = 1:numel(lfp_mod_groups)
        idx = lfp_mod_groups{g} & isfinite(breath_pct_amp);
        histogram(breath_pct_amp(idx), edges_amp, 'Normalization', 'probability', ...
            'FaceColor', lfp_mod_col{g}, 'FaceAlpha', 0.55, 'EdgeColor', 'none')
        leg{end + 1} = sprintf('%s n=%d', lfp_mod_names{g}, sum(idx)); %#ok<AGROW>
    end
    xlabel('Breathing modulation (% of mean rate)')
    ylabel('Probability')
    legend(leg, 'Location', 'best')
    title('Breathing')

    nexttile
    hold on
    for g = 1:numel(lfp_mod_groups)
        idx = lfp_mod_groups{g} & isfinite(delta_pct_amp);
        histogram(delta_pct_amp(idx), edges_amp, 'Normalization', 'probability', ...
            'FaceColor', lfp_mod_col{g}, 'FaceAlpha', 0.55, 'EdgeColor', 'none')
    end
    xlabel('Delta modulation (% of mean rate)')
    ylabel('Probability')
    legend(leg, 'Location', 'best')
    title('Delta LFP')

    nexttile
    hold on
    mod_out = is_xy_tukey_outlier(breath_pct_amp, delta_pct_amp, both_keep, mod_outlier_k);
    mod_ok = both_keep & ~mod_out;
    disp(sprintf('  delta vs breathing %% scatter: dropped %d / %d Tukey outliers (k=%.2g).', ...
        sum(mod_out), sum(both_keep), mod_outlier_k))
    h_mod = gobjects(1, numel(lfp_mod_groups));
    for g = 1:numel(lfp_mod_groups)
        idx = lfp_mod_groups{g} & mod_ok;
        if any(idx)
            h_mod(g) = plot(breath_pct_amp(idx), delta_pct_amp(idx), '.', ...
                'Color', lfp_mod_col{g}, 'MarkerSize', 12);
        else
            h_mod(g) = plot(NaN, NaN, '.', 'Color', lfp_mod_col{g}, 'MarkerSize', 12);
        end
    end
    for g = 1:numel(lfp_mod_groups)
        idx = lfp_mod_groups{g} & mod_ok;
        plot_mean_cross(breath_pct_amp(idx), delta_pct_amp(idx), lfp_mod_col{g});
    end
    xb = breath_pct_amp(mod_ok);
    yb = delta_pct_amp(mod_ok);
    if isempty(xb) || ~any(isfinite(xb)) || ~any(isfinite(yb))
        lims = [0 1];
    else
        lims = [0, max([max(xb(:)), max(yb(:)), 1], [], 'omitnan')];
    end
    plot(lims, lims, ':k')
    xlim(lims); ylim(lims)
    xlabel('Breathing modulation (% of mean rate)')
    ylabel('Delta modulation (% of mean rate)')
    [leg_mod, ttl_mod] = corr_and_paired_rank( ...
        breath_pct_amp, delta_pct_amp, ...
        {trough_lock & mod_ok, peak_lock & mod_ok, unlocked & mod_ok}, ...
        lfp_mod_names, 'breath', 'delta', true);
    legend(h_mod, leg_mod, 'Location', 'best')
    title({'Delta % vs breathing %  (outliers removed; x = group mean)', ttl_mod}, 'FontSize', 9)
    axis square

    nexttile
    hold on
    ppc_keep = valid & isfinite(lfp_ppc) & isfinite(br_ppc);
    h_ppc = gobjects(1, numel(lfp_mod_groups));
    for g = 1:numel(lfp_mod_groups)
        idx = lfp_mod_groups{g} & ppc_keep;
        if any(idx)
            h_ppc(g) = plot(br_ppc(idx), lfp_ppc(idx), '.', ...
                'Color', lfp_mod_col{g}, 'MarkerSize', 12);
        else
            h_ppc(g) = plot(NaN, NaN, '.', 'Color', lfp_mod_col{g}, 'MarkerSize', 12);
        end
    end
    for g = 1:numel(lfp_mod_groups)
        idx = lfp_mod_groups{g} & ppc_keep;
        plot_mean_cross(br_ppc(idx), lfp_ppc(idx), lfp_mod_col{g});
    end
    xlabel('Breathing PPC')
    ylabel('Delta PPC')
    [leg_ppc, ~] = corr_and_paired_rank( ...
        br_ppc, lfp_ppc, ...
        {trough_lock & ppc_keep, peak_lock & ppc_keep, unlocked & ppc_keep}, ...
        lfp_mod_names, 'breath', 'delta', false);
    legend(h_ppc, leg_ppc, 'Location', 'best')
    title('Delta PPC vs breathing PPC  (x = group mean)')
    axis square
    sgtitle(focus_label)
end

%% Example breathing-locked neuron: raster locked to preferred peak or trough
lfp_label = strings(numel(valid), 1);
lfp_label(:) = "unlocked";
lfp_label(trough_lock) = "trough-locked";
lfp_label(peak_lock) = "peak-locked";
br_pref = strings(numel(valid), 1);
br_pref(:) = "breath-unlocked";
br_pref(br_trough_pref) = "breath-trough";
br_pref(br_peak_pref) = "breath-peak";

bidx = find(breath_lock);
lock_rank = 3 * ones(size(bidx));
lock_rank(trough_lock(bidx)) = 1;
lock_rank(peak_lock(bidx)) = 2;
[~, ord_b] = sortrows([lock_rank, cluster_ids(bidx)]);
bidx = bidx(ord_b);

breath_pick_table = table((1:numel(bidx))', session_label(bidx), cluster_ids(bidx), ...
    lfp_label(bidx), br_pref(bidx), brain_region(bidx), brain_area(bidx), ...
    lfp_ppcv(bidx), br_ppcv(bidx), breath_mean_hz(bidx), ...
    breath_pct_at_0(bidx), breath_pct_amp(bidx), ...
    'VariableNames', {'n', 'session', 'id', 'LFP_lock', 'breath_pref', ...
    'region', 'area', 'LFP_PPCp', 'breath_PPCp', 'mean_Hz', ...
    'pct_at_event', 'pct_amp_of_mean'});
disp('Breathing-locked neurons (set breath_raster_n to pick one; NaN = prompt, [] = skip):')
disp(breath_pick_table)

pick_n = breath_raster_n;
if isempty(bidx)
    disp('No breathing-locked neurons; skip example raster.')
elseif ~isempty(pick_n) && isnumeric(pick_n) && isnan(pick_n)
    pick_n = input('Row n from table (0 skip): ');
end
if isempty(bidx) || isempty(pick_n) || ~isfinite(pick_n) || pick_n < 1
    disp('Skip example breathing raster.')
elseif pick_n > numel(bidx)
    warning('breath_raster_n = %g is past the table (%d rows).', pick_n, numel(bidx))
else
    ii = bidx(pick_n);
    area_txt = char(brain_region(ii));
    if strlength(brain_area(ii)) > 0 && brain_area(ii) ~= brain_region(ii) ...
            && brain_area(ii) ~= "unmatched"
        area_txt = [area_txt, ' / ', char(brain_area(ii))];
    end
    title_txt = sprintf('n=%d  id=%g  %s  LFP %s  %s', ...
        pick_n, cluster_ids(ii), area_txt, char(lfp_label(ii)), char(br_pref(ii)));
    if br_trough_pref(ii)
        align_to = 'trough';
    else
        align_to = 'peak';
    end
    ns = sess_idx_all(ii);
    plot_example_breath_raster(acute_struct(ns), cluster_ids(ii), ...
        data_root, title_txt, freq2filter, align_to, breath_mod_win);
end

%% Store counts used in the test
breathing_test.alpha = alpha;
breathing_test.focus_regions = focus_regions;


function plot_entrainment_psth_stack(psth_time, psth, locked, unlocked, angles, ...
    smooth_n, x_lim, c_lim, y_lim, align_label)
% Four-row Fig 3 layout without area columns: locked heatmap, locked mean,
% unlocked heatmap, unlocked mean. Red = peak half-cycle, blue = trough.

    psth_lock = pct_smooth_psth(psth(locked, :), smooth_n);
    psth_unlk = pct_smooth_psth(psth(unlocked, :), smooth_n);
    ang_lock  = angles(locked);
    ang_unlk  = angles(unlocked);
    ang_lock(isnan(ang_lock)) = 0;
    ang_unlk(isnan(ang_unlk)) = 0;

    [sorted_lock, ord_lock] = sort(ang_lock);
    matrix_lock = psth_lock(ord_lock, :);
    y_ticks = 180 * sorted_lock / pi;

    subplot(4, 1, 1)
    if isempty(matrix_lock)
        title([align_label, '  locked (none)'])
    else
        hold on
        imagesc(psth_time, y_ticks, matrix_lock)
        axis xy
        clim(c_lim)
        xlim(x_lim)
        ylim([-180 180])
        yticks([-180 -90 0 90 180])
        plot([psth_time(1) psth_time(end)], [0 0], 'w')
        plot([psth_time(1) psth_time(end)], [90 90], 'w')
        plot([psth_time(1) psth_time(end)], [-90 -90], 'w')
        title({[align_label, '  locked'], [num2str(size(matrix_lock, 1)), ' cells']})
        ylabel('Preferred angle (deg)')
    end

    subplot(4, 1, 2)
    hold on
    if ~isempty(matrix_lock)
        sel = robust_rows(matrix_lock);
        is_trough = sorted_lock > pi/2 | sorted_lock < -pi/2;
        n_tr = sum(is_trough & sel);
        n_pk = sum(~is_trough & sel);
        plot_mean_ci(psth_time, matrix_lock, is_trough & sel, 'b')
        plot_mean_ci(psth_time, matrix_lock, ~is_trough & sel, 'r')
        title(sprintf('Locked mean  trough n=%d  peak n=%d', n_tr, n_pk))
    else
        title('Locked mean (none)')
    end
    xlim(x_lim)
    ylim(y_lim)
    ylabel('% change')

    [sorted_unlk, ord_unlk] = sort(ang_unlk);
    matrix_unlk = psth_unlk(ord_unlk, :);

    subplot(4, 1, 3)
    if isempty(matrix_unlk)
        title([align_label, '  unlocked (none)'])
    else
        hold on
        imagesc(psth_time, 180 * sorted_unlk / pi, matrix_unlk)
        axis xy
        clim(c_lim)
        xlim(x_lim)
        ylim([-180 180])
        yticks([-180 -90 0 90 180])
        plot([psth_time(1) psth_time(end)], [0 0], 'w')
        plot([psth_time(1) psth_time(end)], [90 90], 'w')
        plot([psth_time(1) psth_time(end)], [-90 -90], 'w')
        title({[align_label, '  unlocked'], [num2str(size(matrix_unlk, 1)), ' cells']})
        ylabel('Preferred angle (deg)')
    end

    subplot(4, 1, 4)
    hold on
    if ~isempty(matrix_unlk)
        sel = robust_rows(matrix_unlk);
        is_neg = sorted_unlk > -pi & sorted_unlk < 0;
        plot_mean_ci(psth_time, matrix_unlk, is_neg & sel, 'b')
        plot_mean_ci(psth_time, matrix_unlk, ~is_neg & sel, 'r')
        title('Unlocked mean')
    else
        title('Unlocked mean (none)')
    end
    xlim(x_lim)
    ylim(y_lim)
    xlabel('Time from event (s)')
    ylabel('% change')
end


function psth_pct = pct_smooth_psth(psth, smooth_n)
    psth_pct = psth;
    if isempty(psth_pct)
        return
    end
    if size(psth_pct, 2) == 1
        psth_pct = psth_pct';
    end
    for j = 1:size(psth_pct, 1)
        y = psth_pct(j, :);
        if smooth_n > 1
            y = smooth(y, smooth_n)';
        end
        mu = mean(y, 'omitnan');
        if mu == 0 || isnan(mu)
            psth_pct(j, :) = nan(size(y));
        else
            psth_pct(j, :) = 100 * (y - mu) / mu;
        end
    end
end


function sel = robust_rows(M)
    if isempty(M)
        sel = false(0, 1);
        return
    end
    mu = mean(M, 2, 'omitnan');
    sd = std(M, 0, 2, 'omitnan');
    sd(sd == 0 | isnan(sd)) = 1;
    z = abs((M - mu) ./ sd);
    sel = max(z, [], 2, 'omitnan') < 8;
end


function plot_mean_ci(t, M, idx, col, ls)
    if nargin < 5 || isempty(ls)
        ls = '-';
    end
    if nnz(idx) < 1
        return
    end
    X = M(idx, :);
    X = X(all(isfinite(X), 2), :);
    if isempty(X)
        return
    end
    if size(X, 1) == 1
        plot(t, X, 'Color', col, 'LineStyle', ls, 'LineWidth', 2)
        return
    end
    [~, ~, ci] = ttest(X);
    if size(ci, 2) == numel(t)
        fill([t fliplr(t)], [ci(1, :) fliplr(ci(2, :))], col, 'FaceAlpha', 0.25, 'EdgeColor', 'none')
        hold on
        plot(t, mean(X, 1, 'omitnan'), 'Color', col, 'LineStyle', ls, 'LineWidth', 2)
    end
end


function s = pct_str(p, k, n)
    s = sprintf('%.1f%%  (%d / %d)', 100 * p, k, n);
end


function [p, ci] = binofit_safe(k, n)
    if n < 1
        p = NaN;
        ci = [NaN NaN];
    else
        [p, ci] = binofit(k, n);
    end
end


function out = report_breathing_lock_proportions(valid, trough_lock, peak_lock, delta_lock, breath_lock, name_tag)
    n_all     = sum(valid);
    n_trough  = sum(trough_lock);
    n_peak    = sum(peak_lock);
    n_delta   = sum(delta_lock);

    p_all     = sum(breath_lock & valid) / max(n_all, 1);
    p_trough  = sum(breath_lock & trough_lock) / max(n_trough, 1);
    p_peak    = sum(breath_lock & peak_lock) / max(n_peak, 1);
    p_delta   = sum(breath_lock & delta_lock) / max(n_delta, 1);

    [phat_all, ci_all] = binofit_safe(sum(breath_lock & valid), n_all);
    [phat_tr, ci_tr]   = binofit_safe(sum(breath_lock & trough_lock), n_trough);
    [phat_pk, ci_pk]   = binofit_safe(sum(breath_lock & peak_lock), n_peak);
    [phat_del, ci_del] = binofit_safe(sum(breath_lock & delta_lock), n_delta);

    N_trough_vs_rest = [sum(trough_lock & breath_lock), sum(trough_lock & ~breath_lock); ...
                        sum(~trough_lock & valid & breath_lock), sum(~trough_lock & valid & ~breath_lock)];
    p_fisher = NaN;
    p_fisher_right = NaN;
    p_chi = NaN;
    chi2_tr = NaN;
    odds_ratio = NaN;
    if all(sum(N_trough_vs_rest, 2) > 0) && all(sum(N_trough_vs_rest, 1) > 0)
        [~, p_fisher, stats_fisher] = fishertest(N_trough_vs_rest);
        [~, p_fisher_right] = fishertest(N_trough_vs_rest, 'Tail', 'right');
        odds_ratio = stats_fisher.OddsRatio;
        [~, chi2_tr, p_chi] = crosstab(trough_lock(valid), breath_lock(valid));
    end

    N_trough_vs_peak = [sum(trough_lock & breath_lock), sum(trough_lock & ~breath_lock); ...
                        sum(peak_lock & breath_lock), sum(peak_lock & ~breath_lock)];
    p_fisher_tp = NaN;
    if all(sum(N_trough_vs_peak, 2) > 0) && all(sum(N_trough_vs_peak, 1) > 0)
        [~, p_fisher_tp] = fishertest(N_trough_vs_peak);
    end

    if nargin < 6 || isempty(name_tag)
        hdr = '--- Breathing lock enrichment ---';
        fig_name = 'Breathing lock proportions';
        title_prefix = '';
    else
        hdr = ['--- Breathing lock enrichment  |  ', name_tag, ' ---'];
        fig_name = ['Breathing lock proportions  |  ', name_tag];
        title_prefix = [name_tag, '  |  '];
    end
    disp(hdr)
    disp(['  all neurons:          ', pct_str(p_all, sum(breath_lock & valid), n_all)])
    disp(['  delta-locked:         ', pct_str(p_delta, sum(breath_lock & delta_lock), n_delta)])
    disp(['  trough (delta-lock):  ', pct_str(p_trough, sum(breath_lock & trough_lock), n_trough)])
    disp(['  peak (delta-lock):    ', pct_str(p_peak, sum(breath_lock & peak_lock), n_peak)])
    disp('2x2 trough vs rest (rows: trough / rest; cols: breath+ / breath-):')
    disp(N_trough_vs_rest)
    disp(['  Fisher two-sided p = ', num2str(p_fisher), '   OR = ', num2str(odds_ratio)])
    disp(['  Fisher right-tail p (trough more often breath-locked) = ', num2str(p_fisher_right)])
    disp(['  Chi-square p = ', num2str(p_chi), '   chi2 = ', num2str(chi2_tr)])
    disp(['  Fisher trough vs peak p = ', num2str(p_fisher_tp)])

    figure('Name', fig_name)
    vals = [phat_all; phat_del; phat_tr; phat_pk];
    cis  = [ci_all; ci_del; ci_tr; ci_pk];
    bar(vals, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
    hold on
    errorbar(1:4, vals, vals - cis(:, 1), cis(:, 2) - vals, 'k', 'LineStyle', 'none', 'LineWidth', 1)
    set(gca, 'XTick', 1:4, 'XTickLabel', { ...
        sprintf('All (n=%d)', n_all), ...
        sprintf('Delta-locked (n=%d)', n_delta), ...
        sprintf('Trough (n=%d)', n_trough), ...
        sprintf('Peak (n=%d)', n_peak)})
    ylabel('Proportion breathing-locked')
    ylim([0 1])
    title({[title_prefix, 'Fisher (trough vs rest) p = ', num2str(p_fisher, '%.3g')], ...
           ['Chi-square p = ', num2str(p_chi, '%.3g'), '  |  trough vs peak p = ', num2str(p_fisher_tp, '%.3g')]})

    out = struct();
    out.name = name_tag;
    out.N_trough_vs_rest = N_trough_vs_rest;
    out.p_fisher = p_fisher;
    out.p_fisher_right = p_fisher_right;
    out.odds_ratio = odds_ratio;
    out.p_chi2 = p_chi;
    out.chi2 = chi2_tr;
    out.p_fisher_trough_vs_peak = p_fisher_tp;
    out.proportions = table(vals, cis(:, 1), cis(:, 2), ...
        'VariableNames', {'p', 'ci_lo', 'ci_hi'}, ...
        'RowNames', {'all', 'delta_locked', 'trough', 'peak'});
end


function Z = zscore_rows_smooth(X, smooth_n)
    Z = nan(size(X));
    if isempty(X)
        return
    end
    if size(X, 2) == 1
        X = X';
        Z = nan(size(X));
    end
    for j = 1:size(X, 1)
        y = X(j, :);
        if nargin >= 2 && smooth_n > 1
            y = smooth(y, smooth_n)';
        end
        sd = std(y, 0, 'omitnan');
        mu = mean(y, 'omitnan');
        if ~(sd > 0) || isnan(sd)
            continue
        end
        Z(j, :) = (y - mu) / sd;
    end
end


function y_lim_share = plot_lfp_breath_lock_combos(psth_time, psth_peak, psth_trough, ...
    lfp_unlocked, lfp_trough, lfp_peak, br_trough, br_peak, smooth_n, x_lim, ...
    y_lim_in, name_suffix, do_plot)
% Four LFP × breathing lock combinations. Gray = LFP-unlocked with the same
% breathing preference as that panel.
    if nargin < 11
        y_lim_in = [];
    end
    if nargin < 12 || isempty(name_suffix)
        name_suffix = '';
    end
    if nargin < 13 || isempty(do_plot)
        do_plot = true;
    end
    z_peak = zscore_rows_smooth(psth_peak, smooth_n);
    z_trough = zscore_rows_smooth(psth_trough, smooth_n);
    gray = [0.62 0.62 0.62];
    col_tr = [0.15 0.35 0.85];
    col_pk = [0.85 0.2 0.2];
    combos = {
        lfp_trough & br_trough, lfp_unlocked & br_trough, 'LFP trough  ×  breath trough', col_tr
        lfp_trough & br_peak,   lfp_unlocked & br_peak,   'LFP trough  ×  breath peak',   col_pk
        lfp_peak & br_trough,   lfp_unlocked & br_trough, 'LFP peak  ×  breath trough',   col_tr
        lfp_peak & br_peak,     lfp_unlocked & br_peak,   'LFP peak  ×  breath peak',     col_pk
        };
    align = {z_peak, z_trough};
    align_name = {'peak', 'trough'};
    if isempty(y_lim_in)
        y_all = [];
        for a = 1:2
            Z = align{a};
            for c = 1:4
                mu_u = mean(Z(combos{c, 2} & all(isfinite(Z), 2), :), 1, 'omitnan');
                mu_c = mean(Z(combos{c, 1} & all(isfinite(Z), 2), :), 1, 'omitnan');
                y_all = [y_all; mu_u(:); mu_c(:)]; %#ok<AGROW>
            end
        end
        y_all = y_all(isfinite(y_all));
        if isempty(y_all)
            y_lim_share = [-1 1];
        else
            pad = 0.1 * max(range(y_all), 1);
            y_lim_share = [min(y_all) - pad, max(y_all) + pad];
        end
    else
        y_lim_share = y_lim_in;
    end
    if ~do_plot
        return
    end
    if isempty(name_suffix)
        fig_tag = '';
        title_tag = '';
    else
        fig_tag = ['  ', name_suffix];
        title_tag = [name_suffix, '  |  '];
    end
    axs = gobjects(2, 4);
    for a = 1:2
        Z = align{a};
        figure('Name', ['LFP × breath lock  breathing-', align_name{a}, ' aligned', fig_tag], ...
            'units', 'normalized', 'outerposition', [0.08 0.12 0.72 0.78])
        for c = 1:4
            axs(a, c) = subplot(2, 2, c);
            hold on
            plot_mean_ci(psth_time, Z, combos{c, 2}, gray)
            plot_mean_ci(psth_time, Z, combos{c, 1}, combos{c, 4})
            xline(0, ':k')
            xlim(x_lim)
            ylim(y_lim_share)
            n_c = sum(combos{c, 1});
            n_u = sum(combos{c, 2});
            title(sprintf('%s   n=%d   (gray n=%d)', combos{c, 3}, n_c, n_u))
            if c > 2
                xlabel(['Time from breathing ', align_name{a}, ' (s)'])
            end
            if c == 1 || c == 3
                ylabel('Z-scored rate')
            end
            if c == 1
                legend({'LFP unlocked, same breath lock', 'LFP + breath lock'}, 'Location', 'best')
            end
        end
        sgtitle(sprintf(['%sBreathing %s-aligned.  Gray = LFP-unlocked with the ' ...
            'same breathing lock as the panel'], title_tag, align_name{a}))
    end
    linkaxes(axs(:), 'y');
end


function plot_breathing_peak_vs_trough(psth_time, psth_peak, psth_trough, group_mask, group_name, group_col, smooth_n, x_lim)
% One figure: peak-aligned vs trough-aligned means for the same grouping.
% Breath-phase groups should invert (peak-pref high at peak t=0, trough-pref
% high at trough t=0).
    z_peak = zscore_rows_smooth(psth_peak, smooth_n);
    z_trough = zscore_rows_smooth(psth_trough, smooth_n);
    [~, i0] = min(abs(psth_time));

    figure('Name', 'Breathing peak vs trough alignment (breath-phase groups)', ...
        'units', 'normalized', 'outerposition', [0.1 0.15 0.7 0.7])
    subplot(1, 2, 1)
    hold on
    for g = 1:numel(group_mask)
        plot_mean_ci(psth_time, z_peak, group_mask{g} & all(isfinite(z_peak), 2), group_col{g})
    end
    xline(0, ':k')
    xlim(x_lim)
    xlabel('Time from breathing peak (s)')
    ylabel('Z-scored rate')
    legend(group_name, 'Location', 'best')
    title(sprintf('Peak-aligned   (t=0 z: trough-pref %.2f, peak-pref %.2f)', ...
        mean(z_peak(group_mask{2}, i0), 'omitnan'), mean(z_peak(group_mask{3}, i0), 'omitnan')))

    subplot(1, 2, 2)
    hold on
    for g = 1:numel(group_mask)
        plot_mean_ci(psth_time, z_trough, group_mask{g} & all(isfinite(z_trough), 2), group_col{g})
    end
    xline(0, ':k')
    xlim(x_lim)
    xlabel('Time from breathing trough (s)')
    ylabel('Z-scored rate')
    legend(group_name, 'Location', 'best')
    title(sprintf('Trough-aligned   (t=0 z: trough-pref %.2f, peak-pref %.2f)', ...
        mean(z_trough(group_mask{2}, i0), 'omitnan'), mean(z_trough(group_mask{3}, i0), 'omitnan')))
end


function out = compare_breathing_psth_groups(psth_time, psth, group_mask, group_name, group_col, smooth_n, x_lim, align_label)
    breath_z = zscore_rows_smooth(psth, smooth_n);
    [~, i0] = min(abs(psth_time));
    z0 = breath_z(:, i0);

    figure('Name', ['Breathing ', align_label, ' PSTH by delta-lock group'])
    subplot(1, 2, 1)
    hold on
    z_groups = cell(1, 3);
    for g = 1:3
        idx = group_mask{g} & isfinite(z0);
        z_groups{g} = z0(idx);
        plot_mean_ci(psth_time, breath_z, idx, group_col{g})
    end
    xline(0, ':k')
    xlim(x_lim)
    xlabel(['Time from breathing ', align_label, ' (s)'])
    ylabel('Z-scored rate')
    legend(group_name, 'Location', 'best')
    title(['Breathing-', align_label, ' PSTH (z-scored per neuron)'])

    vals_kw = [];
    lab_kw  = [];
    for g = 1:3
        v = z_groups{g};
        v = v(isfinite(v));
        vals_kw = [vals_kw; v(:)];
        lab_kw  = [lab_kw; g * ones(numel(v), 1)];
    end
    p_kw = NaN;
    if numel(unique(lab_kw)) >= 2 && numel(vals_kw) >= 3
        p_kw = kruskalwallis(vals_kw, lab_kw, 'off');
    end
    pair_p = nan(3, 3);
    for i = 1:3
        for j = i + 1:3
            a = z_groups{i}; a = a(isfinite(a));
            b = z_groups{j}; b = b(isfinite(b));
            if numel(a) >= 2 && numel(b) >= 2
                pair_p(i, j) = ranksum(a, b);
                pair_p(j, i) = pair_p(i, j);
            end
        end
    end
    pair_p_bonf = min(1, pair_p * 3);

    disp(['--- Breathing-', align_label, ' PSTH (z at t=0) ---'])
    disp(['  Kruskal-Wallis p = ', num2str(p_kw)])
    for g = 1:3
        disp(sprintf('  %s: n=%d  mean z=%.3f', group_name{g}, numel(z_groups{g}), mean(z_groups{g}, 'omitnan')))
    end
    disp('  Pairwise rank-sum (raw / Bonferroni):')
    disp(sprintf('    trough vs unlocked: p=%.3g / %.3g', pair_p(1, 2), pair_p_bonf(1, 2)))
    disp(sprintf('    peak vs unlocked:   p=%.3g / %.3g', pair_p(1, 3), pair_p_bonf(1, 3)))
    disp(sprintf('    trough vs peak:     p=%.3g / %.3g', pair_p(2, 3), pair_p_bonf(2, 3)))

    subplot(1, 2, 2)
    hold on
    for g = 1:3
        xj = g + 0.08 * randn(size(z_groups{g}));
        scatter(xj, z_groups{g}, 14, group_col{g}, 'filled', 'MarkerFaceAlpha', 0.45)
        plot(g, mean(z_groups{g}, 'omitnan'), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 7)
    end
    set(gca, 'XTick', 1:3, 'XTickLabel', group_name)
    ylabel(sprintf('Z-scored rate at t = %.3f s', psth_time(i0)))
    title({['Kruskal-Wallis p = ', num2str(p_kw, '%.3g')], ...
           sprintf('trough vs unlock p=%.3g   peak vs unlock p=%.3g   trough vs peak p=%.3g', ...
           pair_p_bonf(1, 2), pair_p_bonf(1, 3), pair_p_bonf(2, 3))})

    out = struct();
    out.align = align_label;
    out.p_kruskalwallis = p_kw;
    out.pairwise_ranksum = pair_p;
    out.pairwise_bonferroni = pair_p_bonf;
    out.z_at_event = z0;
    out.group_names = group_name;
end


function s = psth_mod_of_mean(psth_row, psth_time, smooth_n, win)
% Percent of mean rate modulated in win: amp = 100*(max-min)/(2*mean).
    s = struct('pct_at_0', NaN, 'pct_amp', NaN, 'mean_hz', NaN, 'pct_ptp', NaN);
    y = psth_row(:)';
    t = psth_time(:)';
    if isempty(y) || numel(y) ~= numel(t) || ~any(isfinite(y))
        return
    end
    if nargin >= 3 && ~isempty(smooth_n) && smooth_n > 1
        y = smooth(y, smooth_n)';
    end
    if nargin < 4 || isempty(win)
        in_win = true(size(t));
    else
        in_win = t >= win(1) & t <= win(2);
    end
    yw = y(in_win);
    if ~any(isfinite(yw))
        return
    end
    mu = mean(yw, 'omitnan');
    s.mean_hz = mu;
    if ~(isfinite(mu) && mu > 0)
        return
    end
    [~, i0] = min(abs(t));
    s.pct_at_0 = 100 * (y(i0) - mu) / mu;
    ptp = max(yw, [], 'omitnan') - min(yw, [], 'omitnan');
    s.pct_ptp = 100 * ptp / mu;
    s.pct_amp = 100 * ptp / (2 * mu);
end


function in_focus = neurons_in_focus_list(brain_region, brain_area, focus_regions)
% Token 'PAG' keeps every PAG column. Anything else matches region or area.
    n = numel(brain_region);
    in_focus = false(n, 1);
    if nargin < 3 || isempty(focus_regions)
        in_focus(:) = true;
        return
    end
    region = lower(strtrim(string(brain_region)));
    area = lower(strtrim(string(brain_area)));
    for k = 1:numel(focus_regions)
        tok = char(lower(strtrim(string(focus_regions{k}))));
        if strcmp(tok, 'pag')
            hit = contains(region, 'pag') | contains(area, 'pag');
        elseif strcmp(tok, 'raphe')
            hit = ismember(region, ["raphe", "rhape", "dr"]) ...
                | ismember(area, ["raphe", "rhape", "dr"]) ...
                | contains(region, 'raphe') | contains(area, 'raphe');
        elseif any(strcmp(tok, {'sc', 'supcol', 'superior colliculus'}))
            hit = ismember(region, ["supcol", "sc"]) | ismember(area, ["supcol", "sc"]) ...
                | contains(region, 'supcol') | contains(area, 'supcol');
        else
            hit = region == string(tok) | area == string(tok);
        end
        in_focus = in_focus | hit;
    end
end


function flag = is_xy_tukey_outlier(x, y, keep, k)
% True where x or y is outside [Q1-k*IQR, Q3+k*IQR] among keep.
    if nargin < 4 || isempty(k)
        k = 1.5;
    end
    flag = false(size(x));
    if nargin < 3 || isempty(keep)
        keep = true(size(x));
    end
    keep = keep & isfinite(x) & isfinite(y);
    idx = find(keep);
    if numel(idx) < 4
        return
    end
    ox = tukey_flag(x(idx), k);
    oy = tukey_flag(y(idx), k);
    flag(idx(ox | oy)) = true;
end


function out = tukey_flag(v, k)
    v = v(:);
    out = false(size(v));
    ok = isfinite(v);
    if nnz(ok) < 4
        return
    end
    q = prctile(v(ok), [25 75]);
    iqr = q(2) - q(1);
    if ~(isfinite(iqr) && iqr > 0)
        return
    end
    lo = q(1) - k * iqr;
    hi = q(2) + k * iqr;
    out(ok) = v(ok) < lo | v(ok) > hi;
end


function plot_mean_cross(x, y, col)
    ok = isfinite(x) & isfinite(y);
    if ~any(ok)
        return
    end
    plot(mean(x(ok)), mean(y(ok)), 'x', 'Color', col, ...
        'MarkerSize', 14, 'LineWidth', 2.5)
end


function [leg, ttl] = corr_and_paired_rank(x, y, group_masks, group_names, xname, yname, do_rank)
% Legend: Spearman r_s per group. Title: Wilcoxon signed-rank (y vs x) if do_rank.
    if nargin < 7 || isempty(do_rank)
        do_rank = true;
    end
    n_g = numel(group_masks);
    leg = cell(1, n_g);
    ttl_bits = {};
    xv = x(:);
    yv = y(:);
    for g = 1:n_g
        idx = group_masks{g}(:) & isfinite(xv) & isfinite(yv);
        n = sum(idx);
        r = NaN;
        p_corr = NaN;
        if n >= 3
            [r, p_corr] = corr(xv(idx), yv(idx), 'Type', 'Spearman', 'Rows', 'complete');
        end
        if isfinite(r)
            leg{g} = sprintf('%s  n=%d  r_s=%.2f p=%.3g', group_names{g}, n, r, p_corr);
        else
            leg{g} = sprintf('%s  n=%d', group_names{g}, n);
        end
        if do_rank && n >= 2
            p_sr = NaN;
            try
                p_sr = signrank(xv(idx), yv(idx));
            catch
            end
            dmed = median(yv(idx) - xv(idx), 'omitnan');
            if dmed > 0
                who = sprintf('%s>%s', yname, xname);
            elseif dmed < 0
                who = sprintf('%s>%s', xname, yname);
            else
                who = 'tied';
            end
            if isfinite(p_sr)
                ttl_bits{end + 1} = sprintf('%s %s p=%.3g', group_names{g}, who, p_sr); %#ok<AGROW>
            end
        end
    end
    if do_rank
        if isempty(ttl_bits)
            ttl = 'Wilcoxon signed-rank: n too small';
        else
            ttl = ['Wilcoxon signed-rank  ', strjoin(ttl_bits, '  |  ')];
        end
    else
        ttl = '';
    end
end


function plot_example_breath_raster(sess, cluster_id, data_root, title_txt, ...
    freq2filter, align_to, win)
% Reload spikes + breathing and raster-lock to preferred peak or trough.
    if nargin < 7 || isempty(win)
        win = [-0.4 0.4];
    end
    if nargin < 6 || isempty(align_to)
        align_to = 'peak';
    end
    if nargin < 5 || isempty(freq2filter)
        freq2filter = 8;
    end

    files_folder = fullfile(data_root, 'Acute data', 'population analysis');
    merged_file = fullfile(files_folder, 'MERGED_SESSIONS.mat');
    calls_file = fullfile(files_folder, 'NEW_ALL_CALLS_TOGETHER.mat');
    npx_dir = '';
    if isfield(sess, 'npx_data_dir') && ~isempty(sess.npx_data_dir) ...
            && exist(fullfile(sess.npx_data_dir, 'spike_times.npy'), 'file')
        npx_dir = sess.npx_data_dir;
    else
        cand = fullfile(data_root, 'Acute data', sess.session_folder);
        if exist(fullfile(cand, 'spike_times.npy'), 'file')
            npx_dir = cand;
        end
    end
    if isempty(npx_dir) || ~exist(merged_file, 'file')
        warning('Analyze_acute_phase_locking:exampleRaster', ...
            'Cannot load spikes/breathing for session %s; skip example raster.', ...
            sess.session_folder)
        return
    end

    spike_times = readNPY(fullfile(npx_dir, 'spike_times.npy'));
    spike_clusters = readNPY(fullfile(npx_dir, 'spike_clusters.npy'));
    spikes_npx = double(spike_times(spike_clusters == cluster_id)) / 30000;
    if isempty(spikes_npx)
        warning('Analyze_acute_phase_locking:exampleRaster', ...
            'No spikes for cluster %g in %s.', cluster_id, sess.session_folder)
        return
    end

    load(merged_file);
    if exist(calls_file, 'file')
        load(calls_file);
    end
    if ~exist('ALL_BREATHING', 'var') || ~exist('SYNCH_MODELS', 'var')
        warning('Analyze_acute_phase_locking:exampleRaster', ...
            'MERGED_SESSIONS is missing ALL_BREATHING / SYNCH_MODELS.')
        return
    end
    syn_i = find([SYNCH_MODELS.date] == sess.date ...
        & [SYNCH_MODELS.experiment] == sess.experiment, 1);
    br_i = find([ALL_BREATHING.date] == sess.date ...
        & [ALL_BREATHING.experiment] == sess.experiment, 1);
    if isempty(syn_i) || isempty(br_i)
        warning('Analyze_acute_phase_locking:exampleRaster', ...
            'No SYNCH_MODELS / ALL_BREATHING for %s.', sess.session_folder)
        return
    end
    spikes_audio = predict(SYNCH_MODELS(syn_i).PAG_NPX, spikes_npx(:));

    breath_fs = 1000;
    breath_peak_distance = 0.1;
    br_vals = zscore(ALL_BREATHING(br_i).values(:));
    t_br = ((1:numel(br_vals)) / breath_fs)';
    t_br_audio = predict(ALL_BREATHING(br_i).SynchModel, t_br);
    t_br_audio = t_br_audio(:);
    std_amp_br = std(br_vals);
    [~, pk_idx] = findpeaks(br_vals, 'MinPeakProminence', 0.5 * std_amp_br, ...
        'MinPeakDistance', breath_peak_distance * breath_fs);
    [~, tr_idx] = findpeaks(-br_vals, 'MinPeakProminence', 0.5 * std_amp_br, ...
        'MinPeakDistance', breath_peak_distance * breath_fs);
    peak_times = t_br_audio(pk_idx);
    trough_times = t_br_audio(tr_idx);

    intervals = sess.valid_intervals;
    if isfield(sess, 'peak_margin_left') && ~isempty(sess.peak_margin_left)
        left_m = sess.peak_margin_left;
        right_m = sess.peak_margin_right;
    else
        left_m = -win(1);
        right_m = win(2);
    end
    peak_times = peak_times(times_in_green(peak_times, intervals, right_m, left_m));
    trough_times = trough_times(times_in_green(trough_times, intervals, right_m, left_m));
    spikes_green = spikes_audio(times_in_green(spikes_audio, intervals, 0, 0));

    if strcmpi(align_to, 'trough')
        events = trough_times;
        ev_name = 'breathing trough';
    else
        events = peak_times;
        ev_name = 'breathing peak';
    end
    if isempty(events)
        warning('Analyze_acute_phase_locking:exampleRaster', ...
            'No %s events in green windows for %s.', ev_name, sess.session_folder)
        return
    end

    if isfield(sess, 'psth_edges') && ~isempty(sess.psth_edges)
        edges = sess.psth_edges;
        edges = edges(edges >= win(1) - 1e-9 & edges <= win(2) + 1e-9);
        if numel(edges) < 3
            edges = win(1):0.01:win(2);
        end
    else
        edges = win(1):0.01:win(2);
    end
    bin_w = mean(diff(edges));
    t_psth = edges(2:end) - 0.5 * bin_w;
    C = zeros(numel(events), numel(edges) - 1);
    figure('Name', ['Example breathing raster  ', title_txt], ...
        'units', 'normalized', 'outerposition', [0.2 0.15 0.45 0.7])
    subplot(2, 1, 1)
    hold on
    for j = 1:numel(events)
        t0 = events(j);
        rel = spikes_green(spikes_green >= t0 + win(1) & spikes_green <= t0 + win(2)) - t0;
        if ~isempty(rel)
            plot([rel(:) rel(:)]', (j + [-0.4 0.4])' * ones(1, numel(rel)), 'k')
        end
        C(j, :) = histcounts(rel, edges);
    end
    xline(0, ':k')
    xlim(win)
    ylim([0.5, numel(events) + 0.5])
    ylabel('Breath event')
    title({title_txt, sprintf('Locked to %s  n=%d events', ev_name, numel(events))})

    subplot(2, 1, 2)
    hold on
    rate = mean(C, 1) / bin_w;
    smooth_n = max(1, round((1 / freq2filter) / bin_w));
    if smooth_n > 1
        rate_s = smooth(rate, smooth_n)';
    else
        rate_s = rate;
    end
    plot(t_psth, rate_s, 'k', 'LineWidth', 1.5)
    yline(mean(rate_s, 'omitnan'), ':k')
    xline(0, ':k')
    xlim(win)
    xlabel(['Time from ', ev_name, ' (s)'])
    ylabel('Rate (Hz)')
end


function mask = times_in_green(t, intervals, right_margin, left_margin)
    if nargin < 3
        right_margin = 0;
    end
    if nargin < 4
        left_margin = 0;
    end
    t = t(:);
    mask = false(size(t));
    if isempty(intervals)
        return
    end
    for k = 1:size(intervals, 1)
        mask = mask | (t >= intervals(k, 1) + left_margin ...
            & t <= intervals(k, 2) - right_margin);
    end
end