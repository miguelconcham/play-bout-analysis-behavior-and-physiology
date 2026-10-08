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
fooof_play_results.note = [ ...
    'Onset PSTH (play_bout_onset); time from hist_range/bin_size; ' ...
    'joint z-score of log10 mean spectra [before; after] per session; ' ...
    'FOOOF on 10.^z; stores peak_fit, flat_spectrum, aperiodic_fit, peak_params'];

disp('FOOOF periodic before/after (onset) analysis complete.')
disp('Results stored in workspace variable: fooof_play_results')
