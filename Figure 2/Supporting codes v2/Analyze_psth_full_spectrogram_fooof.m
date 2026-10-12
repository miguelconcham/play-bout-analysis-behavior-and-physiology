%% Analyze_psth_full_spectrogram_fooof
% Per-session FOOOF decomposition of jointly z-scored mean power spectra
% before vs after play-bout ONSET, then plot the change in the periodic
% (peak) component only.
%
% Uses play_bout_onset (not time-wrapped play_bout_tw_this).
% Time indexes are taken from hist_range / bin_size, same as Figure 2.m.
%
% Joint z-scoring (critical):
%   For each session, mean and SD are taken from the concatenated
%   log10 spectra [before; after], so both conditions share one
%   normalization and the before/after difference is preserved.
%
% FOOOF input note:
%   fooof_matlab requires linear power > 0. After joint z-scoring of
%   log10 power, we pass P = 10.^logP_z so FOOOF refits that same
%   z-scored log spectrum.
%
% Stored per condition (before / after):
%   peak_fit, flat_spectrum, aperiodic_fit, peak_params
%
% Extra: aperiodic-subtracted delta onset PSTH (Figure 2 style).
%   Edit delta_band and aperiodic_template ('before'|'after'|'combined').
%
this_file = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Analyze_psth_full_spectrogram_fooof.m';
repo_root = fileparts(fileparts(fileparts(this_file)));
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');

% Shared helpers (includes Custom functions/fooof_matlab.m)
run(fullfile(repo_root, 'add_repo_paths.m'));

%% Load precomputed spectrogram PSTH
disp('loading')
load(fullfile(saving_folder, 'psth_structure_delta_full_spectrogram.mat'), 'psth_structure');
load(fullfile(saving_folder, 'animal_names_delta_full_spectrogram.mat'), 'animal_names');
f = psth_structure(1).f(:);
disp('ready')

%% Time axis from stored hist_range / spectrogram hop (same as Figure 2.m)
bin_size    = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
psth_ranges = psth_structure(1).hist_range;
time        = psth_ranges(1):bin_size:psth_ranges(2) + bin_size;
n_time      = size(psth_structure(1).play_bout_onset, 3);
time        = time(1:n_time);

%% Time ranges relative to play-bout onset (edit these, not bin numbers)
before_range = [-5 0];   % seconds before onset
after_range  = [0 5];    % seconds after onset

before_idx = time >= before_range(1) & time < before_range(2);
after_idx  = time >= after_range(1)  & time < after_range(2);

if ~any(before_idx) || ~any(after_idx)
    error('No time samples in before_range or after_range. Check ranges vs hist_range.');
end

%% Session loop settings
session_list = [1:11, 13];   % same sessions as Analyze_psth_full_spectrogram
animal_index = [1 1 1 2 2 2 3 3 3 4 5 5 6 6];

nSessions = numel(session_list);
nFreq     = numel(f);

% Spectra / FOOOF fits (one row per session)
periodic_before   = nan(nSessions, nFreq);
periodic_after    = nan(nSessions, nFreq);
flat_before       = nan(nSessions, nFreq);
flat_after        = nan(nSessions, nFreq);
aperiodic_before  = nan(nSessions, nFreq);
aperiodic_after   = nan(nSessions, nFreq);
logP_z_before     = nan(nSessions, nFreq);
logP_z_after      = nan(nSessions, nFreq);

% peak_params: variable number of peaks -> cell per session
% each entry is [CF, PW, BW] as returned by fooof_matlab
peak_params_before = cell(nSessions, 1);
peak_params_after  = cell(nSessions, 1);

fooof_before = cell(nSessions, 1);
fooof_after  = cell(nSessions, 1);
session_animal = nan(nSessions, 1);

%% Per-session: onset PSTH -> joint z-score -> FOOOF
for i = 1:nSessions
    fn = session_list(i);
    session_animal(i) = animal_index(fn);

    % Onset spectrogram PSTH: [trials x freqs x time]
    onset = psth_structure(fn).play_bout_onset;

    % Mean LINEAR power over trials and selected time windows
    before_pow = squeeze(mean(mean(onset(:, :, before_idx), 1, 'omitnan'), 3, 'omitnan'));
    after_pow  = squeeze(mean(mean(onset(:, :, after_idx),  1, 'omitnan'), 3, 'omitnan'));

    before_pow = before_pow(:);
    after_pow  = after_pow(:);

    % Guard against non-positive power (FOOOF / log10)
    before_pow = max(before_pow, eps);
    after_pow  = max(after_pow, eps);

    % --- Joint z-score on log10 power (before and after together) ---
    log_before = log10(before_pow);
    log_after  = log10(after_pow);
    joint = [log_before; log_after];
    mu = mean(joint);
    sd = std(joint);
    if sd == 0
        sd = 1;
    end
    log_before_z = (log_before - mu) / sd;
    log_after_z  = (log_after  - mu) / sd;

    logP_z_before(i, :) = log_before_z.';
    logP_z_after(i, :)  = log_after_z.';

    % Linear power for FOOOF = 10^(z-scored log spectrum)
    P_before = 10.^log_before_z;
    P_after  = 10.^log_after_z;

    out_before = fooof_matlab(f, P_before, ...
        'AperiodicMode', 'fixed', ...
        'PeakWidthLimits', [0.5 12]);
    out_after = fooof_matlab(f, P_after, ...
        'AperiodicMode', 'fixed', ...
        'PeakWidthLimits', [0.5 12]);

    fooof_before{i} = out_before;
    fooof_after{i}  = out_after;

    periodic_before(i, :)  = out_before.peak_fit.';
    periodic_after(i, :)   = out_after.peak_fit.';
    flat_before(i, :)      = out_before.flat_spectrum.';
    flat_after(i, :)       = out_after.flat_spectrum.';
    aperiodic_before(i, :) = out_before.aperiodic_fit.';
    aperiodic_after(i, :)  = out_after.aperiodic_fit.';
    peak_params_before{i}  = out_before.peak_params;
    peak_params_after{i}   = out_after.peak_params;
end

%% Periodic-component change (after - before)
periodic_diff = periodic_after - periodic_before;
mean_periodic_diff = mean(periodic_diff, 1, 'omitnan');
sem_periodic_diff  = std(periodic_diff, 0, 1, 'omitnan') ./ ...
    sqrt(sum(isfinite(periodic_diff), 1));

%% Plot: change in periodic power only
figure('Name', 'FOOOF periodic change (after - before play onset)');
hold on
plot(f, periodic_diff.', 'Color', [0.85 0.7 0.85], 'LineWidth', 0.5);
fill([f.' fliplr(f.')], ...
    [mean_periodic_diff + sem_periodic_diff, ...
     fliplr(mean_periodic_diff - sem_periodic_diff)], ...
    'm', 'EdgeColor', 'none', 'FaceAlpha', 0.25);
plot(f, mean_periodic_diff, 'm', 'LineWidth', 2);
yline(0, ':k');
axis tight
xlabel('Frequency (Hz)')
ylabel('\Delta periodic power (log_{10}, z-scored spectrum)')
title(sprintf('Periodic only: after [%g %g] − before [%g %g] s', ...
    after_range(1), after_range(2), before_range(1), before_range(2)))
box off

%% Optional reference: full z-scored spectrum difference vs periodic-only
figure('Name', 'Z-scored spectrum change (reference)');
hold on
spec_diff = logP_z_after - logP_z_before;
plot(f, mean(spec_diff, 1, 'omitnan'), 'k', 'LineWidth', 1.5);
plot(f, mean_periodic_diff, 'm', 'LineWidth', 1.5);
yline(0, ':k');
legend({'Full z-scored spectrum', 'Periodic only'}, 'Location', 'best')
xlabel('Frequency (Hz)')
ylabel('\Delta (after − before)')
title('Full spectrum vs periodic-only change')
box off

%% Extra analysis: LME on periodic difference (animal as random intercept)
% Same pattern as Analyze_psth_full_spectrogram.m:
%   Power ~ 1 + (1|Animal)
% One observation per session; sessions nested in animals via session_animal.
A = periodic_diff;

beta_est    = nan(nFreq, 1);
SE          = nan(nFreq, 1);
CI_low_95   = nan(nFreq, 1);
CI_high_95  = nan(nFreq, 1);
CI_low_99   = nan(nFreq, 1);
CI_high_99  = nan(nFreq, 1);
pValue      = nan(nFreq, 1);

for freq_n = 1:nFreq
    tbl = table(A(:, freq_n), categorical(session_animal), ...
        'VariableNames', {'Power', 'Animal'});
    lme = fitlme(tbl, 'Power ~ 1 + (1|Animal)');

    coef  = lme.Coefficients;
    ci_95 = coefCI(lme, 'Alpha', 0.05);
    ci_99 = coefCI(lme, 'Alpha', 0.01);

    beta_est(freq_n)   = coef.Estimate(1);
    SE(freq_n)         = coef.SE(1);
    CI_low_95(freq_n)  = ci_95(1, 1);
    CI_high_95(freq_n) = ci_95(1, 2);
    CI_low_99(freq_n)  = ci_99(1, 1);
    CI_high_99(freq_n) = ci_99(1, 2);
    pValue(freq_n)     = coef.pValue(1);
end

sigFreq = find(pValue < 0.01);
fprintf('LME significant frequencies (p < 0.01): %d bins\n', numel(sigFreq));

%% Plot LME intercept estimate for periodic Δ (extra; does not replace SEM plot)
figure('Name', 'FOOOF periodic LME (animal random intercept)');
hold on
plot(f, beta_est, 'k', 'LineWidth', 2);
yline(0, 'k--');
if ~isempty(sigFreq)
    plot(f(sigFreq), beta_est(sigFreq), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
end
xlabel('Frequency (Hz)')
ylabel('LME estimate of \Delta periodic power')
title('Periodic \Delta: LME intercept + (1|Animal), p < 0.01 marked')
box off

%% Plot LME estimate vs session-mean periodic difference
figure('Name', 'FOOOF periodic LME vs session mean');
hold on
plot(f, beta_est, 'k', 'LineWidth', 2)
plot(f, mean_periodic_diff, ':k')
if ~isempty(sigFreq)
    plot(f(sigFreq), beta_est(sigFreq), 'ro', 'MarkerFaceColor', 'r')
end
yline(0, 'k--');
legend({'LME estimate', 'Session-mean \Delta', 'p < 0.01'}, 'Location', 'best')
xlabel('Frequency (Hz)')
ylabel('\Delta periodic power')
title('LME vs session-mean periodic difference')
box off

%% ========================================================================
% Extra analysis: aperiodic-subtracted delta onset PSTH (Figure 2 style)
% ========================================================================
% Subtract a session aperiodic spectrum from the time-resolved onset
% spectrogram (log10 power), then average residual over delta_band to get
% a 1D onset PSTH, z-scored/smoothed like Figure 2 delta power plots.
%
% Editable options:
%   delta_band            : frequency band for the PSTH (Hz)
%   aperiodic_template    : 'before' | 'after' | 'combined'
%                           default 'before'
%   baseline_range_psth   : baseline window for trial z-scoring (s)
%   smooth_wind_psth      : movmean window (samples), as in Figure 2
%
% IMPORTANT:
%   Templates are FOOOF aperiodics fit on RAW (non z-scored) mean linear
%   spectra so subtraction matches log10(power) of play_bout_onset.

delta_band = [1 5];                 % Hz — change freely
aperiodic_template = 'before';       % 'before' | 'after' | 'combined'
baseline_range_psth = [-5 -2];      % seconds, Figure 2 style
smooth_wind_psth = 20;
psth_xlim = [-3 3];

delta_freq_idx = f >= delta_band(1) & f <= delta_band(2);
if ~any(delta_freq_idx)
    error('delta_band [%g %g] does not overlap frequency vector f.', ...
        delta_band(1), delta_band(2));
end

baseline_idx_psth = time > baseline_range_psth(1) & time < baseline_range_psth(2);
if ~any(baseline_idx_psth)
    error('baseline_range_psth [%g %g] has no samples in time.', ...
        baseline_range_psth(1), baseline_range_psth(2));
end

all_delta_psth = [];
all_delta_animal = [];
all_delta_session = [];
aperiodic_template_fit = nan(nSessions, nFreq);
aperiodic_raw_before = nan(nSessions, nFreq);
aperiodic_raw_after = nan(nSessions, nFreq);
aperiodic_raw_combined = nan(nSessions, nFreq);

for i = 1:nSessions
    fn = session_list(i);
    onset = psth_structure(fn).play_bout_onset;   % [trials x freqs x time]
    nTrials = size(onset, 1);

    % Mean LINEAR spectra for raw-scale FOOOF templates
    before_pow = squeeze(mean(mean(onset(:, :, before_idx), 1, 'omitnan'), 3, 'omitnan'));
    after_pow  = squeeze(mean(mean(onset(:, :, after_idx),  1, 'omitnan'), 3, 'omitnan'));
    combined_pow = squeeze(mean(mean( ...
        onset(:, :, before_idx | after_idx), 1, 'omitnan'), 3, 'omitnan'));

    before_pow = max(before_pow(:), eps);
    after_pow  = max(after_pow(:), eps);
    combined_pow = max(combined_pow(:), eps);

    out_ap_before = fooof_matlab(f, before_pow, ...
        'AperiodicMode', 'fixed', 'PeakWidthLimits', [0.5 12]);
    out_ap_after = fooof_matlab(f, after_pow, ...
        'AperiodicMode', 'fixed', 'PeakWidthLimits', [0.5 12]);
    out_ap_combined = fooof_matlab(f, combined_pow, ...
        'AperiodicMode', 'fixed', 'PeakWidthLimits', [0.5 12]);

    aperiodic_raw_before(i, :) = out_ap_before.aperiodic_fit.';
    aperiodic_raw_after(i, :) = out_ap_after.aperiodic_fit.';
    aperiodic_raw_combined(i, :) = out_ap_combined.aperiodic_fit.';

    switch lower(aperiodic_template)
        case 'before'
            ap = out_ap_before.aperiodic_fit(:);
        case 'after'
            ap = out_ap_after.aperiodic_fit(:);
        case 'combined'
            ap = out_ap_combined.aperiodic_fit(:);
        otherwise
            error('aperiodic_template must be ''before'', ''after'', or ''combined''.');
    end
    aperiodic_template_fit(i, :) = ap.';

    % Subtract aperiodic from time-wise log10 power
    log_onset = log10(max(onset, eps));                 % [trials x freq x time]
    residual = log_onset - reshape(ap, [1, nFreq, 1]);   % broadcast over trials/time

    % Delta-band mean residual -> onset PSTH [trials x time]
    delta_psth = squeeze(mean(residual(:, delta_freq_idx, :), 2, 'omitnan'));
    if isvector(delta_psth)
        delta_psth = delta_psth(:).';
        if size(delta_psth, 2) ~= n_time
            delta_psth = delta_psth.';
        end
    end

    % Trial-wise baseline z-score + smooth (Figure 2)
    for trial = 1:nTrials
        base = delta_psth(trial, baseline_idx_psth);
        mu_b = mean(base, 'omitnan');
        sd_b = std(base, 0, 'omitnan');
        if ~(isfinite(sd_b) && sd_b > 0)
            sd_b = 1;
        end
        delta_psth(trial, :) = (delta_psth(trial, :) - mu_b) / sd_b;
        delta_psth(trial, :) = movmean(delta_psth(trial, :), smooth_wind_psth);
    end

    all_delta_psth = [all_delta_psth; delta_psth]; %#ok<AGROW>
    all_delta_animal = [all_delta_animal; repmat(session_animal(i), nTrials, 1)]; %#ok<AGROW>
    all_delta_session = [all_delta_session; repmat(i, nTrials, 1)]; %#ok<AGROW>
end

%% Plot aperiodic-subtracted delta PSTH (mean ± SEM across trials)
mean_delta_psth = mean(all_delta_psth, 1, 'omitnan');
sem_delta_psth  = std(all_delta_psth, 0, 1, 'omitnan') ./ ...
    sqrt(sum(isfinite(all_delta_psth), 1));

figure('Name', 'Aperiodic-subtracted delta onset PSTH');
hold on
fill([time fliplr(time)], ...
    [mean_delta_psth + sem_delta_psth, fliplr(mean_delta_psth - sem_delta_psth)], ...
    'k', 'EdgeColor', 'none', 'FaceAlpha', 0.25);
plot(time, mean_delta_psth, 'k', 'LineWidth', 1.5);
xline(0, ':k');
yline(0, ':k');
xlim(psth_xlim);
xlabel('Time from play-bout onset (s)')
ylabel(sprintf('Z-scored residual %g–%g Hz power', delta_band(1), delta_band(2)))
title(sprintf('Aperiodic-subtracted delta PSTH (template = %s)', aperiodic_template))
box off

%% Per-animal mean aperiodic-subtracted delta PSTH
animal_ids = unique(all_delta_animal);
figure('Name', 'Aperiodic-subtracted delta PSTH by animal');
hold on
animal_mean_delta = nan(numel(animal_ids), n_time);
for a = 1:numel(animal_ids)
    rows = all_delta_animal == animal_ids(a);
    animal_mean_delta(a, :) = mean(all_delta_psth(rows, :), 1, 'omitnan');
    plot(time, animal_mean_delta(a, :), 'LineWidth', 1);
end
plot(time, mean(animal_mean_delta, 1, 'omitnan'), 'k', 'LineWidth', 2);
xline(0, ':k');
yline(0, ':k');
xlim(psth_xlim);
xlabel('Time from play-bout onset (s)')
ylabel(sprintf('Z-scored residual %g–%g Hz power', delta_band(1), delta_band(2)))
title(sprintf('Per-animal means (template = %s)', aperiodic_template))
box off

%% Extra: LME across time on aperiodic-subtracted delta PSTH (animal random intercept)
beta_est_psth   = nan(n_time, 1);
pValue_psth     = nan(n_time, 1);
CI_low_95_psth  = nan(n_time, 1);
CI_high_95_psth = nan(n_time, 1);

for t = 1:n_time
    y = all_delta_psth(:, t);
    ok = isfinite(y);
    if nnz(ok) < 3
        continue
    end
    tbl = table(y(ok), categorical(all_delta_animal(ok)), ...
        'VariableNames', {'Power', 'Animal'});
    lme_t = fitlme(tbl, 'Power ~ 1 + (1|Animal)');
    coef = lme_t.Coefficients;
    ci_95 = coefCI(lme_t, 'Alpha', 0.05);
    beta_est_psth(t) = coef.Estimate(1);
    pValue_psth(t) = coef.pValue(1);
    CI_low_95_psth(t) = ci_95(1, 1);
    CI_high_95_psth(t) = ci_95(1, 2);
end

sigTime = find(pValue_psth < 0.01);

figure('Name', 'Aperiodic-subtracted delta PSTH LME');
hold on
plot(time, beta_est_psth, 'k', 'LineWidth', 2);
yline(0, 'k--');
if ~isempty(sigTime)
    plot(time(sigTime), beta_est_psth(sigTime), 'ro', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 4);
end
xline(0, ':k');
xlim(psth_xlim);
xlabel('Time from play-bout onset (s)')
ylabel('LME estimate')
title(sprintf('Delta residual PSTH LME (1|Animal), template = %s', aperiodic_template))
box off

%% Package results (in-memory; save manually if desired)
fooof_play_results = struct();
fooof_play_results.f = f;
fooof_play_results.time = time;
fooof_play_results.before_range = before_range;
fooof_play_results.after_range = after_range;
fooof_play_results.before_idx = before_idx;
fooof_play_results.after_idx = after_idx;
fooof_play_results.session_list = session_list;
fooof_play_results.session_animal = session_animal;
fooof_play_results.animal_names = animal_names;

fooof_play_results.logP_z_before = logP_z_before;
fooof_play_results.logP_z_after = logP_z_after;

fooof_play_results.periodic_before = periodic_before;       % peak_fit
fooof_play_results.periodic_after = periodic_after;
fooof_play_results.flat_spectrum_before = flat_before;
fooof_play_results.flat_spectrum_after = flat_after;
fooof_play_results.aperiodic_fit_before = aperiodic_before;
fooof_play_results.aperiodic_fit_after = aperiodic_after;
fooof_play_results.peak_params_before = peak_params_before; % [CF PW BW]
fooof_play_results.peak_params_after = peak_params_after;

fooof_play_results.periodic_diff = periodic_diff;
fooof_play_results.fooof_before = fooof_before;  % full FOOOF structs
fooof_play_results.fooof_after = fooof_after;

% LME outputs (animal as random intercept), same layout as full-spectrogram Analyze
fooof_play_results.lme.beta_est = beta_est;
fooof_play_results.lme.SE = SE;
fooof_play_results.lme.CI_low_95 = CI_low_95;
fooof_play_results.lme.CI_high_95 = CI_high_95;
fooof_play_results.lme.CI_low_99 = CI_low_99;
fooof_play_results.lme.CI_high_99 = CI_high_99;
fooof_play_results.lme.pValue = pValue;
fooof_play_results.lme.sigFreq = sigFreq;
fooof_play_results.lme.formula = 'Power ~ 1 + (1|Animal)';

% Aperiodic-subtracted delta onset PSTH
fooof_play_results.delta_psth.delta_band = delta_band;
fooof_play_results.delta_psth.aperiodic_template = aperiodic_template;
fooof_play_results.delta_psth.baseline_range = baseline_range_psth;
fooof_play_results.delta_psth.smooth_wind = smooth_wind_psth;
fooof_play_results.delta_psth.psth = all_delta_psth;
fooof_play_results.delta_psth.animal = all_delta_animal;
fooof_play_results.delta_psth.session = all_delta_session;
fooof_play_results.delta_psth.animal_mean = animal_mean_delta;
fooof_play_results.delta_psth.aperiodic_raw_before = aperiodic_raw_before;
fooof_play_results.delta_psth.aperiodic_raw_after = aperiodic_raw_after;
fooof_play_results.delta_psth.aperiodic_raw_combined = aperiodic_raw_combined;
fooof_play_results.delta_psth.aperiodic_template_fit = aperiodic_template_fit;
fooof_play_results.delta_psth.lme.beta_est = beta_est_psth;
fooof_play_results.delta_psth.lme.pValue = pValue_psth;
fooof_play_results.delta_psth.lme.CI_low_95 = CI_low_95_psth;
fooof_play_results.delta_psth.lme.CI_high_95 = CI_high_95_psth;
fooof_play_results.delta_psth.lme.sigTime = sigTime;

fooof_play_results.note = [ ...
    'Onset PSTH (play_bout_onset); time from hist_range/bin_size; ' ...
    'joint z-score of log10 mean spectra [before; after] per session; ' ...
    'FOOOF on 10.^z; stores peak_fit, flat_spectrum, aperiodic_fit, peak_params; ' ...
    'extra LME on periodic_diff with animal random intercept; ' ...
    'extra aperiodic-subtracted delta onset PSTH (raw-scale FOOOF template)'];

disp('FOOOF periodic before/after (onset) analysis complete.')
disp('Aperiodic-subtracted delta PSTH complete.')
disp('Results stored in workspace variable: fooof_play_results')
