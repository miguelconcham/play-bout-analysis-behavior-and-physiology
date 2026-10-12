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
%     .periodic_fit    {session} -> [nFreq x nWin]   (peak_fit)
%     .full_fit        {session} -> [nFreq x nWin]
%     .flat_spectrum   {session} -> [nFreq x nWin]
%     .mean_power      {session} -> [nFreq x nWin]   (linear mean P)
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

this_file = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Analyze_psth_full_spectrogram_fooof_sliding.m';
repo_root = fileparts(fileparts(fileparts(this_file)));
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
        pk_fit(:, w) = out.peak_fit(:);
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
fooof_sliding_results.periodic_fit = periodic_fit_all;
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
    'parameters saved per session and per t_win center'];

%% Quick look: mean aperiodic offset / exponent across sessions
ap_offset = nan(nSessions, nWin);
ap_exponent = nan(nSessions, nWin);
for i = 1:nSessions
    ap = aperiodic_params_all{i};
    if isempty(ap)
        continue
    end
    ap_offset(i, :) = ap(:, 1).';
    ap_exponent(i, :) = ap(:, nApParams).';  % last col = exponent (fixed or knee)
end

figure('Name', 'Sliding FOOOF aperiodic params');
subplot(2, 1, 1)
hold on
plot(t_win, ap_offset.', 'Color', [0.8 0.8 0.8])
plot(t_win, mean(ap_offset, 1, 'omitnan'), 'k', 'LineWidth', 2)
xline(0, ':k')
ylabel('Offset')
title(sprintf('Aperiodic offset (T = %g s, dt = %g s)', T, dt))
box off

subplot(2, 1, 2)
hold on
plot(t_win, ap_exponent.', 'Color', [0.8 0.8 0.8])
plot(t_win, mean(ap_exponent, 1, 'omitnan'), 'k', 'LineWidth', 2)
xline(0, ':k')
xlabel('Window center (s from play-bout onset)')
ylabel('Exponent')
title('Aperiodic exponent')
box off

%% Quick look: example session periodic fit heatmap
example_sess = 1;
figure('Name', 'Sliding FOOOF periodic fit (example session)');
imagesc(t_win, f, periodic_fit_all{example_sess})
axis xy
colorbar
xlabel('Window center (s)')
ylabel('Frequency (Hz)')
title(sprintf('Periodic (peak) fit — session %d', session_list(example_sess)))

%% Optional save
if save_results
    out_path = fullfile(saving_folder, save_name);
    disp(['saving ' out_path])
    save(out_path, 'fooof_sliding_results', '-v7.3');
end

disp('Sliding-window FOOOF complete.')
disp('Results in workspace: fooof_sliding_results')
disp('  aperiodic_params{session}(win,:)  and  peak_params{session}{win}')
disp('  indexed by t_win (window centers)')
