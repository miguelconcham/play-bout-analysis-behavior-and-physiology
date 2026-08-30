%% Analyze_acute_phase_locking
% Load pooled GENERATE_ACUTE_PHASE_LOCKING results. Plot delta phase
% distributions and peak/trough PSTHs (Fig 3-style heatmaps, no area split).
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
session_label      = {};
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
    session_label    = [session_label; repmat({acute_struct(ns).session_folder}, n, 1)]; %#ok<AGROW>
end

i_ang  = find(strcmp(phase_stat_names, 'PreferedAngle'), 1);
i_ppc  = find(strcmp(phase_stat_names, 'PPC'), 1);
i_ppcv = find(strcmp(phase_stat_names, 'PPCPval'), 1);

lfp_ang    = lfp_stats(:, i_ang);
lfp_ppc    = lfp_stats(:, i_ppc);
lfp_ppcv   = lfp_stats(:, i_ppcv);
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

disp('--- Breathing lock enrichment ---')
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

figure('Name', 'Breathing lock proportions')
vals = [phat_all; phat_del; phat_tr; phat_pk];
cis  = [ci_all; ci_del; ci_tr; ci_pk];
b = bar(vals, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
hold on
errorbar(1:4, vals, vals - cis(:, 1), cis(:, 2) - vals, 'k', 'LineStyle', 'none', 'LineWidth', 1)
set(gca, 'XTick', 1:4, 'XTickLabel', { ...
    sprintf('All (n=%d)', n_all), ...
    sprintf('Delta-locked (n=%d)', n_delta), ...
    sprintf('Trough (n=%d)', n_trough), ...
    sprintf('Peak (n=%d)', n_peak)})
ylabel('Proportion breathing-locked')
ylim([0 1])
title({['Fisher (trough vs rest) p = ', num2str(p_fisher, '%.3g')], ...
       ['Chi-square p = ', num2str(p_chi, '%.3g'), '  |  trough vs peak p = ', num2str(p_fisher_tp, '%.3g')]})

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
plot_lfp_breath_lock_combos(psth_time, breath_psth_peak, breath_psth_trough, ...
    unlocked, trough_lock, peak_lock, br_trough_pref, br_peak_pref, ...
    smooth_n, x_lim);

%% Store counts used in the test
breathing_test = struct();
breathing_test.alpha = alpha;
breathing_test.N_trough_vs_rest = N_trough_vs_rest;
breathing_test.p_fisher = p_fisher;
breathing_test.p_fisher_right = p_fisher_right;
breathing_test.odds_ratio = odds_ratio;
breathing_test.p_chi2 = p_chi;
breathing_test.chi2 = chi2_tr;
breathing_test.p_fisher_trough_vs_peak = p_fisher_tp;
breathing_test.proportions = table(vals, cis(:, 1), cis(:, 2), ...
    'VariableNames', {'p', 'ci_lo', 'ci_hi'}, ...
    'RowNames', {'all', 'delta_locked', 'trough', 'peak'});


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


function plot_lfp_breath_lock_combos(psth_time, psth_peak, psth_trough, ...
    lfp_unlocked, lfp_trough, lfp_peak, br_trough, br_peak, smooth_n, x_lim)
% Four LFP × breathing lock combinations. Gray = LFP-unlocked with the same
% breathing preference as that panel.
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
    axs = gobjects(2, 4);
    for a = 1:2
        Z = align{a};
        figure('Name', ['LFP × breath lock  breathing-', align_name{a}, ' aligned'], ...
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
        sgtitle(sprintf(['Breathing %s-aligned.  Gray = LFP-unlocked with the ' ...
            'same breathing lock as the panel'], align_name{a}))
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
