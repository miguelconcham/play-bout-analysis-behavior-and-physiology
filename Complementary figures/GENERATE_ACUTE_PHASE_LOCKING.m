function acute_struct = GENERATE_ACUTE_PHASE_LOCKING(npx_data_dir, Hd_freq, saving_folder)
% GENERATE_ACUTE_PHASE_LOCKING  Green-window LFP and breathing phase locking.
%
%   acute_struct = GENERATE_ACUTE_PHASE_LOCKING(npx_data_dir, Hd_freq, saving_folder)
%
% Inputs
%   npx_data_dir   Acute session folder named YYYYMMDD_N (continuous.dat + kilosort).
%   Hd_freq        digitalFilter for the LFP band (designfilt bandpass).
%   saving_folder  Folder for LFP/breathing interval figures (created if needed).
%
% Other parameters (PSTH bins, surrogates, call-free windows, breathing) are
% set inside this function. Returns a struct; does not write the results .mat.
%
% LFP: per channel with units, uniformized Hilbert phase on green samples,
%   MVL/PPC + surrogate p-values, PSTH at oscillation peaks and troughs.
% Breathing: same stats/PSTHs on z-scored breathing, independent of Hd_freq.
% Spike times are mapped to audio with predict(PAG_NPX, spikes/30000).
% Green windows last 2-10 s (start = last call + 3 s, end = next stim).
% PSTH events are inset by the PSTH half-window on both sides so a peak
% near the left green border does not get an empty pre-event count.
%
% Example
%   acute_struct = GENERATE_ACUTE_PHASE_LOCKING(pwd, Hd_freq, saving_folder);

if nargin < 1 || isempty(npx_data_dir)
    error('GENERATE_ACUTE_PHASE_LOCKING:npx_data_dir', 'Provide the acute session folder.');
end
if nargin < 3 || isempty(saving_folder)
    saving_folder = '';
end

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
files_folder = [data_root, '\Acute data\population analysis\'];

sr_LFP = 2500;
if nargin < 2 || isempty(Hd_freq)
    Hd_freq = designfilt('bandpassfir', ...
        'FilterOrder', 500, ...
        'CutoffFrequency1', 1, ...
        'CutoffFrequency2', 5, ...
        'SampleRate', sr_LFP);
end

%% Parameters
n_rand = 500;
breath_fs = 1000;
call_free_interval = 3;
min_interval_length = 2;
max_interval_length = 10;
range2plot = [-10 1];
breath_peak_distance = 0.1;
psth_bins = 0.01;
psth_range = round(1.25 * [-1 1] ./ Hd_freq.CutoffFrequency1, 2);
psth_edges = psth_range(1):psth_bins:psth_range(2);
psth_time = psth_edges(2:end) - 0.5 * psth_bins;
n_psth_bins = numel(psth_edges) - 1;
peak_margin_left = -psth_range(1);
peak_margin_right = psth_range(2);
phase_stat_names = {'PreferedAngle', 'MVL', 'MVLPval', 'PPC', 'PPCPval', 'MeanRate'};
freq_range = [Hd_freq.CutoffFrequency1 Hd_freq.CutoffFrequency2];

if ~isempty(saving_folder) && ~exist(saving_folder, 'dir')
    mkdir(saving_folder);
end

%% Session from folder name (YYYYMMDD_N)
[~, session_folder] = fileparts(npx_data_dir);
session_tok = regexp(session_folder, '^(\d+)_(\d+)$', 'tokens', 'once');
if isempty(session_tok)
    error('GENERATE_ACUTE_PHASE_LOCKING:session', ...
        'Folder "%s" is not an Acute data session (expected YYYYMMDD_N).', session_folder);
end
date_2analyze = str2double(session_tok{1});
experiment_2analyze = str2double(session_tok{2});
disp(['Session: date=', num2str(date_2analyze), ...
    ', experiment=', num2str(experiment_2analyze), ...
    ', freq=', num2str(freq_range(1)), '-', num2str(freq_range(2)), ' Hz'])
disp(['Green windows: ', num2str(min_interval_length), '-', num2str(max_interval_length), ...
    ' s after ', num2str(call_free_interval), ' s quiet. PSTH inset L/R = ', ...
    num2str(peak_margin_left), '/', num2str(peak_margin_right), ' s'])

%% Population tables (stim, calls, synch, breathing)
load([files_folder, 'MERGED_SESSIONS.mat']);
load([files_folder, 'NEW_ALL_CALLS_TOGETHER.mat']);

STIM_BEG_END = ALL_STIMULATION_TIMES(ALL_STIMULATION_TIMES.Date==date_2analyze ...
    & ALL_STIMULATION_TIMES.Experiment==experiment_2analyze, {'StimStart','StimEnd'});
CALL_BEG_END = ALL_CALLS_TOGETHER(ALL_CALLS_TOGETHER.Date==date_2analyze ...
    & ALL_CALLS_TOGETHER.Experiment==experiment_2analyze, {'BeginTimes','EndTimes'});
Synch_model = SYNCH_MODELS(find([SYNCH_MODELS.date]==date_2analyze ...
    & [SYNCH_MODELS.experiment]==experiment_2analyze));
if isempty(Synch_model)
    error('GENERATE_ACUTE_PHASE_LOCKING:synch', 'No SYNCH_MODELS entry for this session.');
end
LFP_synch_model = Synch_model.PAG_NPX;

session_index_breathin = find([ALL_BREATHING.date]==date_2analyze ...
    & [ALL_BREATHING.experiment]==experiment_2analyze);
if isempty(session_index_breathin)
    error('GENERATE_ACUTE_PHASE_LOCKING:breathing', 'No ALL_BREATHING entry for this session.');
end

%% LFP and spike sorting
file_pointer = fopen(fullfile(npx_data_dir, 'continuous.dat'), 'r');
LFP = fread(file_pointer, 'int16');
fclose(file_pointer);
LFP = reshape(LFP, 384, numel(LFP)/384);

spike_times = readNPY(fullfile(npx_data_dir, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(npx_data_dir, 'spike_clusters.npy'));
cluster_data = tdfread(fullfile(npx_data_dir, 'cluster_info.tsv'));
keep = ismember(cluster_data.group, {'good', 'mua'});
clusters_list = cluster_data.id(keep);
channel_list = cluster_data.ch(keep);
n_neu = numel(clusters_list);
lfp_ch_1based = double(channel_list) + 1;
channels_with_neurons = unique(lfp_ch_1based)';

lfp_time = ((1:size(LFP, 2)) / sr_LFP)';
lfp_time_audio = predict(LFP_synch_model, lfp_time);

%% Call-free intervals (green windows)
last_call_before_stim_start = nan(size(STIM_BEG_END, 1) - 1, 1);
for stim = 2:size(STIM_BEG_END, 1)
    last_call = max(max(CALL_BEG_END.EndTimes(CALL_BEG_END.EndTimes < STIM_BEG_END.StimStart(stim))), ...
        STIM_BEG_END.StimEnd(stim - 1));
    if isempty(last_call)
        last_call = 0;
    end
    last_call_before_stim_start(stim - 1) = last_call;
end
no_stim_intervals = [(last_call_before_stim_start + call_free_interval) STIM_BEG_END.StimStart(2:end)];
valid_intervals = no_stim_intervals( ...
    no_stim_intervals(:, 2) - no_stim_intervals(:, 1) > min_interval_length ...
    & no_stim_intervals(:, 2) - no_stim_intervals(:, 1) < max_interval_length, :);

valid_lfp_indexes = false(size(lfp_time_audio));
[~, order] = sort(valid_intervals(:, 2) - valid_intervals(:, 1));
for j = 1:size(valid_intervals, 1)
    valid_lfp_indexes(lfp_time_audio >= valid_intervals(j, 1) ...
        & lfp_time_audio < valid_intervals(j, 2)) = true;
end
valid_lfp_indexes = valid_lfp_indexes(:);

plot_interval_stack(LFP(1, :), lfp_time_audio, valid_intervals, order, ...
    CALL_BEG_END, range2plot, peak_margin_left, peak_margin_right, ...
    saving_folder, 'LFP_period_selection', false);

%% LFP phase locking and peak/trough PSTH (one band = Hd_freq)
lfp_phase_stats = nan(n_neu, 6);
lfp_psth_peak = nan(n_neu, n_psth_bins);
lfp_psth_trough = nan(n_neu, n_psth_bins);
green_dur = sum(valid_intervals(:, 2) - valid_intervals(:, 1));
green_lfp = valid_lfp_indexes;
peak_distance_samp = sr_LFP / Hd_freq.CutoffFrequency2;

disp(['LFP channels with neurons: ', num2str(numel(channels_with_neurons))])
for ch_i = 1:numel(channels_with_neurons)
    ch_n = channels_with_neurons(ch_i);
    neu_idx = find(lfp_ch_1based == ch_n);
    disp(['Channel ', num2str(ch_n), ' (', num2str(ch_i), '/', ...
        num2str(numel(channels_with_neurons)), '), ', ...
        num2str(numel(neu_idx)), ' neuron(s)'])

    filtered_LFP = filtfilt(Hd_freq, double(LFP(ch_n, :)));
    hilbert_LFP = hilbert(filtered_LFP);
    phase_raw = angle(hilbert_LFP);
    phase_raw = phase_raw(:);
    amp = abs(hilbert_LFP);
    amp = amp(:);
    if numel(green_lfp) ~= numel(phase_raw) || sum(green_lfp) < 10
        continue
    end
    phase_u = nan(size(phase_raw));
    phase_u(green_lfp) = circular_uniformize(phase_raw(green_lfp));
    pool_phases = phase_u(green_lfp);

    std_amp = std(amp(green_lfp));
    [~, peak_idx] = findpeaks(filtered_LFP(:)', 'MinPeakProminence', 0.5 * std_amp, ...
        'MinPeakDistance', peak_distance_samp);
    [~, trough_idx] = findpeaks(-filtered_LFP(:)', 'MinPeakProminence', 0.5 * std_amp, ...
        'MinPeakDistance', peak_distance_samp);
    peak_audio = predict(LFP_synch_model, lfp_time(peak_idx));
    trough_audio = predict(LFP_synch_model, lfp_time(trough_idx));
    peak_audio = peak_audio(times_in_intervals(peak_audio, valid_intervals, peak_margin_right, peak_margin_left));
    trough_audio = trough_audio(times_in_intervals(trough_audio, valid_intervals, peak_margin_right, peak_margin_left));
    t_lfp_green = lfp_time(green_lfp);
    ph_lfp_green = phase_u(green_lfp);

    for ni = 1:numel(neu_idx)
        dn = neu_idx(ni);
        id = clusters_list(dn);
        spikes_npx = double(spike_times(spike_clusters == id)) / 30000;
        if isempty(spikes_npx)
            continue
        end
        spikes_audio = predict(LFP_synch_model, spikes_npx(:));
        in_green = times_in_intervals(spikes_audio, valid_intervals, 0);
        spikes_green_npx = spikes_npx(in_green);
        spikes_green_audio = spikes_audio(in_green);
        this_phases = interp1(t_lfp_green(:), ph_lfp_green(:), spikes_green_npx(:), 'linear', NaN);
        this_phases = this_phases(~isnan(this_phases));

        stats = phase_lock_stats(this_phases, pool_phases, n_rand);
        stats(6) = numel(this_phases) / max(green_dur, eps);
        lfp_phase_stats(dn, :) = stats;
        lfp_psth_peak(dn, :) = mean_event_psth(spikes_green_audio, peak_audio, psth_edges);
        lfp_psth_trough(dn, :) = mean_event_psth(spikes_green_audio, trough_audio, psth_edges);
    end
end

%% Breathing peaks/troughs and interval figure
breathinv_values = ALL_BREATHING(session_index_breathin).values;
time_br = ((1:numel(breathinv_values)) / breath_fs)';
breathing_time_audio = predict(ALL_BREATHING(session_index_breathin).SynchModel, time_br);
breathing_time_audio = breathing_time_audio(:);
breathinv_values = zscore(breathinv_values(:));
std_amp_br = std(breathinv_values);
[~, breathing_peaks] = findpeaks(breathinv_values, 'MinPeakProminence', 0.5 * std_amp_br, ...
    'MinPeakDistance', breath_peak_distance * breath_fs);
[~, breathing_troughs] = findpeaks(-breathinv_values, 'MinPeakProminence', 0.5 * std_amp_br, ...
    'MinPeakDistance', breath_peak_distance * breath_fs);
peak_times_all = breathing_time_audio(breathing_peaks);
trough_times_all = breathing_time_audio(breathing_troughs);
in_green_peaks = times_in_intervals(peak_times_all, valid_intervals, peak_margin_right, peak_margin_left);
in_green_troughs = times_in_intervals(trough_times_all, valid_intervals, peak_margin_right, peak_margin_left);

plot_interval_stack(breathinv_values, breathing_time_audio, valid_intervals, order, ...
    CALL_BEG_END, range2plot, peak_margin_left, peak_margin_right, saving_folder, ...
    'Breathing_period_selection', true, peak_times_all, trough_times_all);

breathing_peaks = breathing_peaks(in_green_peaks);
breathing_troughs = breathing_troughs(in_green_troughs);
disp(['Breathing peaks kept: ', num2str(numel(breathing_peaks)), ' / ', num2str(numel(in_green_peaks))])
disp(['Breathing troughs kept: ', num2str(numel(breathing_troughs)), ' / ', num2str(numel(in_green_troughs))])

%% Breathing phase locking and peak/trough PSTH
phase_br_raw = angle(hilbert(breathinv_values));
green_br = times_in_intervals(breathing_time_audio, valid_intervals, 0);
phase_br_u = nan(size(phase_br_raw));
if sum(green_br) >= 10
    phase_br_u(green_br) = circular_uniformize(phase_br_raw(green_br));
end
pool_br = phase_br_u(green_br);
peak_times_audio = breathing_time_audio(breathing_peaks);
trough_times_audio = breathing_time_audio(breathing_troughs);
t_br_green = breathing_time_audio(green_br);
ph_br_green = phase_br_u(green_br);

breathing_phase_stats = nan(n_neu, 6);
breathing_psth_peak = nan(n_neu, n_psth_bins);
breathing_psth_trough = nan(n_neu, n_psth_bins);

for dn = 1:n_neu
    id = clusters_list(dn);
    spikes_npx = double(spike_times(spike_clusters == id)) / 30000;
    if isempty(spikes_npx)
        continue
    end
    spikes_audio = predict(LFP_synch_model, spikes_npx(:));
    in_green = times_in_intervals(spikes_audio, valid_intervals, 0);
    spikes_green_audio = spikes_audio(in_green);
    this_phases = interp1(t_br_green(:), ph_br_green(:), spikes_green_audio(:), 'linear', NaN);
    this_phases = this_phases(~isnan(this_phases));
    stats = phase_lock_stats(this_phases, pool_br, n_rand);
    stats(6) = numel(this_phases) / max(green_dur, eps);
    breathing_phase_stats(dn, :) = stats;
    breathing_psth_peak(dn, :) = mean_event_psth(spikes_green_audio, peak_times_audio, psth_edges);
    breathing_psth_trough(dn, :) = mean_event_psth(spikes_green_audio, trough_times_audio, psth_edges);
end

%% Output struct (caller saves later)
acute_struct = struct();
acute_struct.session_folder = session_folder;
acute_struct.date = date_2analyze;
acute_struct.experiment = experiment_2analyze;
acute_struct.freq_range = freq_range;
acute_struct.npx_data_dir = npx_data_dir;
acute_struct.saving_folder = saving_folder;
acute_struct.clusters_list = clusters_list;
acute_struct.channel_list = channel_list;
acute_struct.valid_intervals = valid_intervals;
acute_struct.phase_stat_names = phase_stat_names;
acute_struct.psth_time = psth_time;
acute_struct.psth_edges = psth_edges;
acute_struct.psth_range = psth_range;
acute_struct.peak_margin_left = peak_margin_left;
acute_struct.peak_margin_right = peak_margin_right;
acute_struct.n_rand = n_rand;
acute_struct.lfp_phase_stats = lfp_phase_stats;
acute_struct.lfp_psth_peak = lfp_psth_peak;
acute_struct.lfp_psth_trough = lfp_psth_trough;
acute_struct.breathing_phase_stats = breathing_phase_stats;
acute_struct.breathing_psth_peak = breathing_psth_peak;
acute_struct.breathing_psth_trough = breathing_psth_trough;
end

function plot_interval_stack(signal, t_audio, valid_intervals, order, CALL_BEG_END, ...
    range2plot, peak_margin_left, peak_margin_right, saving_folder, fig_stem, ...
    show_peaks, peak_times_all, trough_times_all)
% Stacked traces aligned to stim onset; green = analysis window.
% Red dots = peaks, blue dots = troughs (optional). Dotted lines mark the
% PSTH inset (events too close to either green border are not used for PSTH).
    if nargin < 11
        show_peaks = false;
    end
    if nargin < 13
        trough_times_all = [];
    end
    if nargin < 7 || isempty(peak_margin_left)
        peak_margin_left = 0;
    end
    signal = signal(:);
    t_audio = t_audio(:);
    figure('units', 'normalized', 'outerposition', [0 0 .25 1]);
    hold on
    for j = 1:size(valid_intervals, 1)
        stim_start = valid_intervals(order(j), 2);
        interval_start = valid_intervals(order(j), 1);
        last_calls = CALL_BEG_END.EndTimes(CALL_BEG_END.EndTimes < stim_start ...
            & CALL_BEG_END.EndTimes > stim_start + range2plot(1));
        idx = t_audio < stim_start + range2plot(2) & t_audio > (stim_start + range2plot(1));
        seg = signal(idx);
        seg_min = min(seg);
        seg_rng = range(seg);
        if seg_rng == 0
            seg_rng = 1;
        end
        y = (seg - seg_min) / seg_rng;
        plot(t_audio(idx) - stim_start, y + j - 1, 'k')
        if show_peaks
            peaks_in_plot = peak_times_all > stim_start + range2plot(1) ...
                & peak_times_all < stim_start + range2plot(2);
            peak_t = peak_times_all(peaks_in_plot);
            peak_amp = interp1(t_audio, signal, peak_t);
            peak_y = (peak_amp - seg_min) / seg_rng + j - 1;
            kept_here = false(size(peak_t));
            if ~isempty(peak_t)
                kept_here = times_in_intervals(peak_t, valid_intervals(order(j), :), peak_margin_right, peak_margin_left);
            end
            plot(peak_t(kept_here) - stim_start, peak_y(kept_here), '.r', 'MarkerSize', 10)
            plot(peak_t(~kept_here) - stim_start, peak_y(~kept_here), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 8)
        end
        if ~isempty(trough_times_all)
            troughs_in_plot = trough_times_all > stim_start + range2plot(1) ...
                & trough_times_all < stim_start + range2plot(2);
            trough_t = trough_times_all(troughs_in_plot);
            trough_amp = interp1(t_audio, signal, trough_t);
            trough_y = (trough_amp - seg_min) / seg_rng + j - 1;
            kept_tr = false(size(trough_t));
            if ~isempty(trough_t)
                kept_tr = times_in_intervals(trough_t, valid_intervals(order(j), :), peak_margin_right, peak_margin_left);
            end
            plot(trough_t(kept_tr) - stim_start, trough_y(kept_tr), '.b', 'MarkerSize', 10)
            plot(trough_t(~kept_tr) - stim_start, trough_y(~kept_tr), '.', 'Color', [0.5 0.6 0.85], 'MarkerSize', 8)
        end
        if ~isempty(last_calls)
            plot([last_calls'; last_calls'] - stim_start, [0 1] + j - 1, 'r')
        end
        fill([interval_start stim_start stim_start interval_start] - stim_start, ...
            [0 0 1 1] + j - 1, 'g', 'FaceAlpha', .5, 'Edgecolor', 'none')
        left_cut = interval_start - stim_start + peak_margin_left;
        plot([left_cut left_cut], [0 1] + j - 1, ':k')
        plot([-peak_margin_right -peak_margin_right], [0 1] + j - 1, ':k')
    end
    axis tight
    xlabel('Time from stim onset (s)')
    if ~isempty(trough_times_all)
        title({strrep(fig_stem, '_', ' '), 'red = peaks,  blue = troughs'})
    else
        title(strrep(fig_stem, '_', ' '))
    end
    if ~isempty(saving_folder)
        saveas(gcf, fullfile(saving_folder, [fig_stem, '.jpg']))
        saveas(gcf, fullfile(saving_folder, [fig_stem, '.svg']))
    end
end

function mask = times_in_intervals(t, intervals, right_margin, left_margin)
    if nargin < 3
        right_margin = 0;
    end
    if nargin < 4
        left_margin = 0;
    end
    t = t(:);
    mask = false(size(t));
    for k = 1:size(intervals, 1)
        mask = mask | (t >= intervals(k, 1) + left_margin ...
            & t <= intervals(k, 2) - right_margin);
    end
end

function stats = phase_lock_stats(spike_phases, pool_phases, n_rand)
    stats = nan(1, 6);
    spike_phases = spike_phases(~isnan(spike_phases(:)));
    pool_phases = pool_phases(~isnan(pool_phases(:)));
    if isempty(spike_phases)
        return
    end
    stats(1) = circ_mean(spike_phases);
    stats(2) = circ_r(spike_phases);
    if numel(spike_phases) < 2
        return
    end
    stats(4) = compute_ppc_fast(spike_phases);
    n_spk = numel(spike_phases);
    if numel(pool_phases) < n_spk || n_rand < 1
        return
    end
    rand_mvl = nan(n_rand, 1);
    rand_ppc = nan(n_rand, 1);
    replace_flag = numel(pool_phases) < n_spk;
    for nr = 1:n_rand
        samp = randsample(pool_phases, n_spk, replace_flag);
        rand_mvl(nr) = circ_r(samp);
        rand_ppc(nr) = compute_ppc_fast(samp);
    end
    stats(3) = sum(rand_mvl > stats(2)) / n_rand;
    stats(5) = sum(rand_ppc > stats(4)) / n_rand;
end

function mu = mean_event_psth(spikes, events, edges)
    n_bins = numel(edges) - 1;
    bin_w = edges(2) - edges(1);
    if isempty(events) || isempty(spikes)
        mu = nan(1, n_bins);
        return
    end
    C = zeros(numel(events), n_bins);
    lo = edges(1);
    hi = edges(end);
    spikes = spikes(:);
    for j = 1:numel(events)
        t0 = events(j);
        rel = spikes(spikes >= t0 + lo & spikes <= t0 + hi) - t0;
        C(j, :) = histcounts(rel, edges);
    end
    mu = mean(C, 1) / bin_w;
end
