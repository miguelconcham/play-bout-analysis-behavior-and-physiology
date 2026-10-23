%% Analyze_psth_full_spectrogram_fooof_sliding
% Slide one window of duration T along the play-bout ONSET spectrogram and
% run FOOOF at each window center (raw mean spectrum over trials × time in
% the window). Saves full aperiodic and periodic parameters for every
% session and every window center.
%
% Uses play_bout_onset [trials x freqs x time], not time-wrapped data.
% Time axis from hist_range / bin_size (same as Figure 2.m).
%
% Editable options:
%   T              : window duration (s), centered on each t_win
%   dt             : step between window centers (s)
%   time_lim       : optional [t_min t_max] to restrict centers ( [] = all )
%   AperiodicMode, PeakWidthLimits, MaxNPeaks, ... (FOOOF name-value)
%
% Outputs (workspace):
%   fooof_sliding_results
%     .f, .time, .t_win, .T, .dt, .session_list, .session_animal
%     .aperiodic_fit   {session} -> [nFreq x nWin]
%     .periodic_fit    {session} -> [nFreq x nWin]   (Gaussian peak_fit)
%     .periodic_data   {session} -> [nFreq x nWin]   (logP - aperiodic_fit)
%     .full_fit        {session} -> [nFreq x nWin]
%     .flat_spectrum   {session} -> [nFreq x nWin]
%     .mean_power      {session} -> [nFreq x nWin]   (linear mean P)
%     .log_power       {session} -> [nFreq x nWin]   (log10 mean P)
%     .aperiodic_params{session} -> [nWin x nApParams]
%         fixed: [offset, exponent]
%         knee : [offset, knee, exponent]
%     .peak_params     {session}{win} -> [nPeaks x 3]  [CF, PW, BW]
%     .gaussian_params {session}{win} -> [nPeaks x 3]  [CF, height, SD]
%     .r_squared       {session} -> [1 x nWin]
%     .error           {session} -> [1 x nWin]
%     .fooof_out       {session}{win} -> full fooof_matlab struct
%
% Optional save:
%   set save_results = true to write a .mat under saving_folder.

this_file = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Event aligned LFP\Analyze_psth_full_spectrogram_fooof_sliding.m';
repo_root = fileparts(fileparts(fileparts(fileparts(this_file))));
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');

run(fullfile(repo_root, 'add_repo_paths.m'));

%% Load precomputed spectrogram PSTH
disp('loading')
load(fullfile(saving_folder, 'psth_structure_delta_full_spectrogram.mat'), 'psth_structure');
load(fullfile(saving_folder, 'animal_names_delta_full_spectrogram.mat'), 'animal_names');
f = psth_structure(1).f(:);
disp('ready')

%% Time axis (same as Figure 2 / Analyze_psth_full_spectrogram_fooof)
bin_size    = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
psth_ranges = psth_structure(1).hist_range;
time        = psth_ranges(1):bin_size:psth_ranges(2) + bin_size;
n_time      = size(psth_structure(1).play_bout_onset, 3);
time        = time(1:n_time);

%% Sliding-window options (edit these)
T  = 5;          % window duration (s), centered
dt = 0.5;        % step between centers (s)
time_lim = [];   % e.g. [-10 10], or [] for full range where a full window fits

% FOOOF options
aperiodic_mode   = 'fixed';     % 'fixed' or 'knee'
peak_width_lims  = [0.5 12];
max_n_peaks      = Inf;
min_peak_height  = 0;
peak_threshold   = 2;

% Periodic delta-band summary (subplot 3)
delta_band = [1 5];              % Hz — mean periodic power in this band
baseline_range_z = [-5 -2];      % s (window centers) for baseline z-score

save_results = false;
save_name = sprintf('fooof_sliding_T%g_dt%g.mat', T, dt);

%% Session list
session_list = [1:11, 13];
animal_index = [1 1 1 2 2 2 3 3 3 4 5 5 6 6];
nSessions = numel(session_list);
nFreq = numel(f);
session_animal = animal_index(session_list).';

%% Window centers: only positions where [t-T/2, t+T/2] is inside time
t_min_full = time(1) + T/2;
t_max_full = time(end) - T/2;
if ~isempty(time_lim)
    t_min_full = max(t_min_full, time_lim(1));
    t_max_full = min(t_max_full, time_lim(2));
end
if t_max_full <= t_min_full
    error('No valid window centers for T = %g with current time axis / time_lim.', T);
end
t_win = t_min_full:dt:t_max_full;
nWin = numel(t_win);

fprintf('Sliding FOOOF: T = %g s, dt = %g s, %d centers from %.2f to %.2f s\n', ...
    T, dt, nWin, t_win(1), t_win(end));

%% Preallocate per-session containers
aperiodic_fit_all   = cell(nSessions, 1);
periodic_fit_all    = cell(nSessions, 1);
periodic_data_all   = cell(nSessions, 1);  % logP - aperiodic_fit (data, not Gaussian)
log_power_all       = cell(nSessions, 1);
full_fit_all        = cell(nSessions, 1);
flat_spectrum_all   = cell(nSessions, 1);
mean_power_all      = cell(nSessions, 1);
aperiodic_params_all = cell(nSessions, 1);
peak_params_all     = cell(nSessions, 1);
gaussian_params_all = cell(nSessions, 1);
r_squared_all       = cell(nSessions, 1);
error_all           = cell(nSessions, 1);
fooof_out_all       = cell(nSessions, 1);

nApParams = 2 + strcmpi(aperiodic_mode, 'knee');  % fixed:2, knee:3

%% Per session × window center
for i = 1:nSessions
    fn = session_list(i);
    fprintf('Session %d / %d (psth index %d)\n', i, nSessions, fn);

    onset = psth_structure(fn).play_bout_onset;  % [trials x freqs x time]

    ap_fit = nan(nFreq, nWin);
    pk_fit = nan(nFreq, nWin);
    pk_data = nan(nFreq, nWin);
    log_P = nan(nFreq, nWin);
    fl_fit = nan(nFreq, nWin);
    ft_fit = nan(nFreq, nWin);
    mean_P = nan(nFreq, nWin);
    ap_params = nan(nWin, nApParams);
    peak_params_win = cell(nWin, 1);
    gauss_params_win = cell(nWin, 1);
    r2_win = nan(1, nWin);
    err_win = nan(1, nWin);
    fooof_win = cell(nWin, 1);

    for w = 1:nWin
        t0 = t_win(w);
        win_idx = time >= (t0 - T/2) & time < (t0 + T/2);
        if ~any(win_idx)
            continue
        end

        % Mean LINEAR power over trials and times in the window
        P = squeeze(mean(mean(onset(:, :, win_idx), 1, 'omitnan'), 3, 'omitnan'));
        P = max(P(:), eps);
        mean_P(:, w) = P;

        try
            out = fooof_matlab(f, P, ...
                'AperiodicMode', aperiodic_mode, ...
                'PeakWidthLimits', peak_width_lims, ...
                'MaxNPeaks', max_n_peaks, ...
                'MinPeakHeight', min_peak_height, ...
                'PeakThreshold', peak_threshold);
        catch ME
            warning('FOOOF failed session %d, t = %.3f s: %s', i, t0, ME.message);
            continue
        end

        fooof_win{w} = out;
        ap_fit(:, w) = out.aperiodic_fit(:);
        pk_fit(:, w) = out.peak_fit(:);   % Gaussian periodic model
        log_P(:, w) = out.log_power(:);
        % Aperiodic-corrected spectrum (real data residual, not Gaussian)
        pk_data(:, w) = out.log_power(:) - out.aperiodic_fit(:);
        fl_fit(:, w) = out.flat_spectrum(:);
        ft_fit(:, w) = out.full_fit(:);

        ap_row = out.aperiodic_params(:).';
        ap_params(w, 1:numel(ap_row)) = ap_row;

        peak_params_win{w} = out.peak_params;       % [CF PW BW]
        gauss_params_win{w} = out.gaussian_params;  % [CF height SD]
        r2_win(w) = out.r_squared;
        err_win(w) = out.error;
    end

    aperiodic_fit_all{i} = ap_fit;
    periodic_fit_all{i} = pk_fit;
    periodic_data_all{i} = pk_data;
    log_power_all{i} = log_P;
    full_fit_all{i} = ft_fit;
    flat_spectrum_all{i} = fl_fit;
    mean_power_all{i} = mean_P;
    aperiodic_params_all{i} = ap_params;
    peak_params_all{i} = peak_params_win;
    gaussian_params_all{i} = gauss_params_win;
    r_squared_all{i} = r2_win;
    error_all{i} = err_win;
    fooof_out_all{i} = fooof_win;
end

%% Package results
fooof_sliding_results = struct();
fooof_sliding_results.f = f;
fooof_sliding_results.time = time;
fooof_sliding_results.t_win = t_win(:);
fooof_sliding_results.T = T;
fooof_sliding_results.dt = dt;
fooof_sliding_results.time_lim = time_lim;
fooof_sliding_results.bin_size = bin_size;
fooof_sliding_results.session_list = session_list;
fooof_sliding_results.session_animal = session_animal;
fooof_sliding_results.animal_names = animal_names;
fooof_sliding_results.settings.AperiodicMode = aperiodic_mode;
fooof_sliding_results.settings.PeakWidthLimits = peak_width_lims;
fooof_sliding_results.settings.MaxNPeaks = max_n_peaks;
fooof_sliding_results.settings.MinPeakHeight = min_peak_height;
fooof_sliding_results.settings.PeakThreshold = peak_threshold;

fooof_sliding_results.aperiodic_fit = aperiodic_fit_all;
fooof_sliding_results.periodic_fit = periodic_fit_all;     % Gaussian peak_fit
fooof_sliding_results.periodic_data = periodic_data_all;   % logP - aperiodic_fit
fooof_sliding_results.log_power = log_power_all;
fooof_sliding_results.full_fit = full_fit_all;
fooof_sliding_results.flat_spectrum = flat_spectrum_all;
fooof_sliding_results.mean_power = mean_power_all;

% Parameters indexed by session and window center
fooof_sliding_results.aperiodic_params = aperiodic_params_all;  % {sess}: [nWin x nAp]
fooof_sliding_results.peak_params = peak_params_all;            % {sess}{win}: [nPeaks x 3]
fooof_sliding_results.gaussian_params = gaussian_params_all;    % {sess}{win}: [nPeaks x 3]
fooof_sliding_results.r_squared = r_squared_all;
fooof_sliding_results.error = error_all;
fooof_sliding_results.fooof_out = fooof_out_all;

fooof_sliding_results.note = [ ...
    'Sliding centered window FOOOF on raw mean onset spectrogram; ' ...
    'parameters saved per session and per t_win center; ' ...
    'periodic_data = log_power - aperiodic_fit at each window'];

%% Quick look: aperiodic params + baseline-z-scored aperiodic-corrected delta power
ap_offset = nan(nSessions, nWin);
ap_exponent = nan(nSessions, nWin);
for i = 1:nSessions
    ap = aperiodic_params_all{i};
    if isempty(ap)
        continue
    end
    ap_offset(i, :) = ap(:, 1).';
    % Session-wise z-score of offset (as in local modified version)
    ap_offset(i, :) = (ap_offset(i, :) - mean(ap_offset(i, :), 'omitnan')) / ...
        std(ap_offset(i, :), 0, 'omitnan');
    ap_exponent(i, :) = ap(:, nApParams).';  % last col = exponent (fixed or knee)
end

% Mean aperiodic-corrected spectrum in delta_band, then baseline z-score
% (data - aperiodic, NOT the Gaussian peak_fit)
delta_freq_idx = f >= delta_band(1) & f <= delta_band(2);
if ~any(delta_freq_idx)
    error('delta_band [%g %g] does not overlap frequency vector f.', ...
        delta_band(1), delta_band(2));
end
baseline_idx_z = t_win >= baseline_range_z(1) & t_win < baseline_range_z(2);
if ~any(baseline_idx_z)
    error('baseline_range_z [%g %g] has no window centers in t_win.', ...
        baseline_range_z(1), baseline_range_z(2));
end

periodic_delta = nan(nSessions, nWin);
for i = 1:nSessions
    pk = periodic_data_all{i};   % [nFreq x nWin], logP - aperiodic_fit
    if isempty(pk)
        continue
    end
    periodic_delta(i, :) = mean(pk(delta_freq_idx, :), 1, 'omitnan');
    base = periodic_delta(i, baseline_idx_z);
    mu_b = mean(base, 'omitnan');
    sd_b = std(base, 0, 'omitnan');
    if ~(isfinite(sd_b) && sd_b > 0)
        sd_b = 1;
    end
    periodic_delta(i, :) = (periodic_delta(i, :) - mu_b) / sd_b;
end

fooof_sliding_results.delta_band = delta_band;
fooof_sliding_results.baseline_range_z = baseline_range_z;
fooof_sliding_results.ap_offset_z = ap_offset;
fooof_sliding_results.ap_exponent = ap_exponent;
fooof_sliding_results.periodic_delta_z = periodic_delta;

%% Quick look figure: mean ± t-test CI + top significance bars (p < 0.05)
alpha_sig = 0.05;

[~, p_off, ci_off] = ttest(ap_offset);
[~, p_exp, ci_exp] = ttest(ap_exponent);
[~, p_del, ci_del] = ttest(periodic_delta);

fooof_sliding_results.ttest.ap_offset.p = p_off;
fooof_sliding_results.ttest.ap_offset.ci = ci_off;
fooof_sliding_results.ttest.ap_exponent.p = p_exp;
fooof_sliding_results.ttest.ap_exponent.ci = ci_exp;
fooof_sliding_results.ttest.periodic_delta.p = p_del;
fooof_sliding_results.ttest.periodic_delta.ci = ci_del;

figure('Name', 'Sliding FOOOF aperiodic params + periodic delta (ttest CI)');
subplot(3, 1, 1)
plot_trace_ci_sigbars(t_win, ap_offset, ...
    mean(ap_offset, 1, 'omitnan'), ci_off, p_off, alpha_sig);
ylabel('Offset (z)')
title(sprintf('Aperiodic offset (T = %g s, dt = %g s)', T, dt))

subplot(3, 1, 2)
plot_trace_ci_sigbars(t_win, ap_exponent, ...
    mean(ap_exponent, 1, 'omitnan'), ci_exp, p_exp, alpha_sig);
ylabel('Exponent')
title('Aperiodic exponent')

subplot(3, 1, 3)
plot_trace_ci_sigbars(t_win, periodic_delta, ...
    mean(periodic_delta, 1, 'omitnan'), ci_del, p_del, alpha_sig);
xlabel('Window center (s from play-bout onset)')
ylabel(sprintf('Data−aperiodic %g–%g Hz (z)', delta_band(1), delta_band(2)))
title(sprintf('Aperiodic-corrected delta power (baseline z [%g %g] s)', ...
    baseline_range_z(1), baseline_range_z(2)))

%% Extra: LME mean ± CI with animal as random intercept (session nested in animal)
% Same pattern as Figure 2 / Analyze_psth_full_spectrogram:
%   Power ~ 1 + (1|Animal)
% Applied at each window center to:
%   1) aperiodic offset (z)
%   2) aperiodic exponent
%   3) aperiodic-corrected delta power (z)

measures = {ap_offset, ap_exponent, periodic_delta};
measure_names = {'ap_offset', 'ap_exponent', 'periodic_delta'};
measure_ylabels = { ...
    'Offset (z)', ...
    'Exponent', ...
    sprintf('Data−aperiodic %g–%g Hz (z)', delta_band(1), delta_band(2))};
measure_titles = { ...
    'LME aperiodic offset', ...
    'LME aperiodic exponent', ...
    'LME aperiodic-corrected delta power'};

lme_results = struct();

figure('Name', 'Sliding FOOOF LME (animal random intercept)');
for m = 1:numel(measures)
    A = measures{m};
    beta_est = nan(1, nWin);
    SE = nan(1, nWin);
    CI_low = nan(1, nWin);
    CI_high = nan(1, nWin);
    pValue = nan(1, nWin);

    for w = 1:nWin
        y = A(:, w);
        ok = isfinite(y);
        if nnz(ok) < 3
            continue
        end
        tbl = table(y(ok), categorical(session_animal(ok)), ...
            'VariableNames', {'Power', 'Animal'});
        lme = fitlme(tbl, 'Power ~ 1 + (1|Animal)');
        coef = lme.Coefficients;
        ci95 = coefCI(lme, 'Alpha', 0.05);
        beta_est(w) = coef.Estimate(1);
        SE(w) = coef.SE(1);
        CI_low(w) = ci95(1, 1);
        CI_high(w) = ci95(1, 2);
        pValue(w) = coef.pValue(1);
    end

    lme_results.(measure_names{m}).beta_est = beta_est;
    lme_results.(measure_names{m}).SE = SE;
    lme_results.(measure_names{m}).CI_low_95 = CI_low;
    lme_results.(measure_names{m}).CI_high_95 = CI_high;
    lme_results.(measure_names{m}).pValue = pValue;
    lme_results.(measure_names{m}).formula = 'Power ~ 1 + (1|Animal)';

    subplot(3, 1, m)
    plot_lme_ci_sigbars(t_win, beta_est, CI_low, CI_high, pValue, alpha_sig);
    ylabel(measure_ylabels{m})
    title(measure_titles{m})
    if m == numel(measures)
        xlabel('Window center (s from play-bout onset)')
    end
end

fooof_sliding_results.lme = lme_results;
fooof_sliding_results.lme.session_animal = session_animal;
fooof_sliding_results.lme.alpha = alpha_sig;

%% Quick look: example session aperiodic-corrected spectrum heatmap
example_sess = 1;
figure('Name', 'Sliding FOOOF aperiodic-corrected spectrum (example session)');
imagesc(t_win, f, periodic_data_all{example_sess})
axis xy
colorbar
xlabel('Window center (s)')
ylabel('Frequency (Hz)')
title(sprintf('logP − aperiodic (session %d)', session_list(example_sess)))

%% Optional save
if save_results
    out_path = fullfile(saving_folder, save_name);
    disp(['saving ' out_path])
    save(out_path, 'fooof_sliding_results', '-v7.3');
end

disp('Sliding-window FOOOF complete.')
disp('Results in workspace: fooof_sliding_results')
disp('  periodic_data{session}: [nFreq x nWin] = logP - aperiodic_fit')
disp('  periodic_fit{session}:  Gaussian peak_fit (still stored)')
disp('  periodic_delta_z: baseline-z-scored mean(periodic_data) in delta_band')
disp('  ttest / lme summaries stored under fooof_sliding_results.ttest and .lme')
disp('  indexed by t_win (window centers)')

%% ========================================================================
% Local helpers
% ========================================================================
function plot_trace_ci_sigbars(t_win, session_mat, mean_trace, ci, p, alpha_sig)
% Session traces (grey), mean (black), t-test CI shade, top sig bars (p < alpha).

hold on
plot(t_win, session_mat.', 'Color', [0.8 0.8 0.8])
ok = isfinite(ci(1, :)) & isfinite(ci(2, :));
if any(ok)
    fill([t_win(ok) fliplr(t_win(ok))], ...
        [ci(1, ok) fliplr(ci(2, ok))], ...
        'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
end
plot(t_win, mean_trace, 'k', 'LineWidth', 2)
xline(0, ':k')
yline(0, ':k')
axis tight
add_top_sig_bars(t_win, p < alpha_sig)
box off
end

function plot_lme_ci_sigbars(t_win, beta_est, CI_low, CI_high, pValue, alpha_sig)
% LME intercept estimate ± 95% CI, top sig bars (p < alpha).

hold on
ok = isfinite(CI_low) & isfinite(CI_high) & isfinite(beta_est);
if any(ok)
    fill([t_win(ok) fliplr(t_win(ok))], ...
        [CI_low(ok) fliplr(CI_high(ok))], ...
        'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
end
plot(t_win, beta_est, 'k', 'LineWidth', 2)
xline(0, ':k')
yline(0, ':k')
axis tight
add_top_sig_bars(t_win, pValue < alpha_sig)
box off
end

function add_top_sig_bars(t_win, sig)
% Shaded rectangles for contiguous significant periods in upper 10% of ylim.

if ~any(sig)
    return
end
sig = sig(:).';
sig(~isfinite(sig)) = false;
y_lim = ylim;
y_lo = y_lim(1) + 0.9 * (y_lim(2) - y_lim(1));
y_hi = y_lim(2);

d = diff([false, logical(sig), false]);
starts = find(d == 1);
ends = find(d == -1) - 1;
for j = 1:numel(starts)
    i0 = starts(j);
    i1 = ends(j);
    if i0 == i1
        % single bin: give a small width from neighboring spacing
        if numel(t_win) > 1
            half = median(diff(t_win)) / 2;
        else
            half = 0.05;
        end
        x = [t_win(i0)-half, t_win(i0)+half, t_win(i0)+half, t_win(i0)-half];
    else
        x = t_win([i0 i1 i1 i0]);
    end
    fill(x, [y_lo y_lo y_hi y_hi], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.55);
end
end
