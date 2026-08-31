function acute_struct = GENERATE_STIM_CALL_RESPONSES(npx_data_dir, Hd_freq, saving_folder)
% GENERATE_STIM_CALL_RESPONSES  Call and green-window PSTHs plus phases.
%
%   acute_struct = GENERATE_STIM_CALL_RESPONSES(npx_data_dir, Hd_freq, saving_folder)
%
% Same session identity fields as GENERATE_ACUTE_PHASE_LOCKING (date,
% experiment, clusters_list, channel_list, valid_intervals) so results can
% be joined later. Call times are already audio; spikes are mapped with
% predict(PAG_NPX, spike_times/30000).
%
% PSTHs (10 ms bins, ±0.5 s), from a session-wide 10 ms rate that is
% z-scored after setting samples within stim_blank of StimEnd to NaN:
%   call onset / offset: drop calls whose onset is within stim_blank after
%     StimEnd (too close to stimulation end)
%   green-window: aligned at green_start + 0.5 s (looks at [0 1] s from
%   green onset, to check a slow stim-related trend)
%
% Phase at each green-window *onset* (not the +0.5 s align time):
%   LFP Hilbert (Hd_freq) on that neuron's channel
%   breathing Hilbert (one trace for the session)

if nargin < 1 || isempty(npx_data_dir)
    error('GENERATE_STIM_CALL_RESPONSES:npx_data_dir', 'Provide the acute session folder.');
end
if nargin < 3 || isempty(saving_folder)
    saving_folder = '';
end

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
files_folder = [data_root, '\Acute data\population analysis\'];
load([files_folder, 'MERGED_SESSIONS.mat']);
load([files_folder, 'NEW_ALL_CALLS_TOGETHER.mat']);

sr_LFP = 2500;
if nargin < 2 || isempty(Hd_freq)
    Hd_freq = designfilt('bandpassfir', ...
        'FilterOrder', 500, ...
        'CutoffFrequency1', 1, ...
        'CutoffFrequency2', 5, ...
        'SampleRate', sr_LFP);
end

%% Parameters (same green-window rule as GENERATE_ACUTE_PHASE_LOCKING)
breath_fs = 1000;
call_free_interval = 3;
min_interval_length = 2;
max_interval_length = 10;
psth_bins = 0.01;
psth_range = [-0.5 0.5];
green_align_offset = 0.5;
% Stim artifact window (s). Change this one value (default 0.1 = 100 ms):
%   - drop calls whose onset is in [StimEnd, StimEnd + stim_blank)
%   - set rate samples within stim_blank of StimEnd to NaN before z-score
stim_blank = 0.1;
stim_offset_nan = stim_blank;
psth_edges = psth_range(1):psth_bins:psth_range(2);
psth_time = psth_edges(2:end) - 0.5 * psth_bins;
n_psth_bins = numel(psth_edges) - 1;
freq_range = [Hd_freq.CutoffFrequency1 Hd_freq.CutoffFrequency2];
phase_stat_names = {'PreferedAngle', 'MVL', 'MVLPval', 'PPC', 'PPCPval', 'MeanRate'};

%% Session from folder name (YYYYMMDD_N)
[~, session_folder] = fileparts(npx_data_dir);
session_tok = regexp(session_folder, '^(\d+)_(\d+)$', 'tokens', 'once');
if isempty(session_tok)
    error('GENERATE_STIM_CALL_RESPONSES:session', ...
        'Folder "%s" is not an Acute data session (expected YYYYMMDD_N).', session_folder);
end
date_2analyze = str2double(session_tok{1});
experiment_2analyze = str2double(session_tok{2});
disp(['STIM/CALL session: date=', num2str(date_2analyze), ...
    ', experiment=', num2str(experiment_2analyze)])

%% Population tables
STIM_BEG_END = ALL_STIMULATION_TIMES(ALL_STIMULATION_TIMES.Date==date_2analyze ...
    & ALL_STIMULATION_TIMES.Experiment==experiment_2analyze, {'StimStart','StimEnd'});
CALL_BEG_END = ALL_CALLS_TOGETHER(ALL_CALLS_TOGETHER.Date==date_2analyze ...
    & ALL_CALLS_TOGETHER.Experiment==experiment_2analyze, {'BeginTimes','EndTimes'});
Synch_model = SYNCH_MODELS(find([SYNCH_MODELS.date]==date_2analyze ...
    & [SYNCH_MODELS.experiment]==experiment_2analyze));
if isempty(Synch_model)
    error('GENERATE_STIM_CALL_RESPONSES:synch', 'No SYNCH_MODELS entry for this session.');
end
LFP_synch_model = Synch_model.PAG_NPX;

session_index_breathin = find([ALL_BREATHING.date]==date_2analyze ...
    & [ALL_BREATHING.experiment]==experiment_2analyze);
if isempty(session_index_breathin)
    error('GENERATE_STIM_CALL_RESPONSES:breathing', 'No ALL_BREATHING entry for this session.');
end

%% LFP and kilosort
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

lfp_time = ((1:size(LFP, 2)) / sr_LFP)';
lfp_time_audio = predict(LFP_synch_model, lfp_time);

%% Green windows (same rule as GENERATE_ACUTE_PHASE_LOCKING)
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
n_win = size(valid_intervals, 1);
green_onset = valid_intervals(:, 1);
green_align = green_onset + green_align_offset;
disp(['Green windows: ', num2str(n_win)])

%% Call times (already audio): drop onsets too close after StimEnd
stim_start = STIM_BEG_END.StimStart(:);
stim_end = STIM_BEG_END.StimEnd(:);
call_onset_all = CALL_BEG_END.BeginTimes(:);
call_offset_all = CALL_BEG_END.EndTimes(:);
near_stim_end = false(size(call_onset_all));
for s = 1:numel(stim_end)
    near_stim_end = near_stim_end | (call_onset_all >= stim_end(s) ...
        & call_onset_all < stim_end(s) + stim_blank);
end
keep_call = ~near_stim_end;
call_onset = call_onset_all(keep_call);
call_offset = call_offset_all(keep_call);
disp(['Calls: ', num2str(numel(call_onset)), ' / ', num2str(numel(call_onset_all)), ...
    ' (onset >= ', num2str(1000 * stim_blank), ' ms after StimEnd)'])

%% Breathing phase at green onset (shared across neurons)
breathinv_values = zscore(ALL_BREATHING(session_index_breathin).values(:));
time_br = ((1:numel(breathinv_values)) / breath_fs)';
breathing_time_audio = predict(ALL_BREATHING(session_index_breathin).SynchModel, time_br);
breathing_time_audio = breathing_time_audio(:);
phase_br = angle(hilbert(breathinv_values));
green_onset_breath_phase = nan(1, n_win);
if n_win > 0
    green_onset_breath_phase = interp1(breathing_time_audio, phase_br(:), green_onset, 'linear', NaN)';
end

%% LFP phase at green onset, per channel then per neuron
green_onset_lfp_phase = nan(n_neu, n_win);
channels_with_neurons = unique(lfp_ch_1based)';
for ch_i = 1:numel(channels_with_neurons)
    ch_n = channels_with_neurons(ch_i);
    neu_idx = find(lfp_ch_1based == ch_n);
    if n_win < 1
        continue
    end
    filtered_LFP = filtfilt(Hd_freq, double(LFP(ch_n, :)));
    phase_lfp = angle(hilbert(filtered_LFP));
    phase_lfp = phase_lfp(:);
    ph_at = interp1(lfp_time_audio(:), phase_lfp, green_onset, 'linear', NaN);
    green_onset_lfp_phase(neu_idx, :) = repmat(ph_at(:)', numel(neu_idx), 1);
end

%% Spike PSTHs from session rate, NaN near stim offset, then z-score
call_psth_onset = nan(n_neu, n_psth_bins);
call_psth_offset = nan(n_neu, n_psth_bins);
green_psth = nan(n_neu, n_psth_bins);
mean_rate = nan(n_neu, 1);

if ~isempty(lfp_time_audio)
    rec_t0 = lfp_time_audio(1);
    rec_t1 = lfp_time_audio(end);
    rec_dur = rec_t1 - rec_t0;
else
    rec_t0 = NaN;
    rec_t1 = NaN;
    rec_dur = NaN;
end

for dn = 1:n_neu
    id = clusters_list(dn);
    spikes_npx = double(spike_times(spike_clusters == id)) / 30000;
    if isempty(spikes_npx)
        continue
    end
    spikes_audio = predict(LFP_synch_model, spikes_npx(:));
    if rec_dur > 0
        mean_rate(dn) = numel(spikes_audio) / rec_dur;
    end
    [rate, rate_t] = binned_rate(spikes_audio, rec_t0, rec_t1, psth_bins);
    rate = nan_near_times(rate, rate_t, stim_end, stim_offset_nan);
    rate_z = nan_zscore(rate);
    call_psth_onset(dn, :) = mean_snippets(rate_z, rate_t, call_onset, psth_time);
    call_psth_offset(dn, :) = mean_snippets(rate_z, rate_t, call_offset, psth_time);
    green_psth(dn, :) = mean_snippets(rate_z, rate_t, green_align, psth_time);
end

%% Output (identity fields match GENERATE_ACUTE_PHASE_LOCKING)
acute_struct = struct();
acute_struct.analysis = 'stim_call_responses';
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
acute_struct.green_align_offset = green_align_offset;
acute_struct.green_align_times = green_align;
acute_struct.n_rand = NaN;
acute_struct.stim_blank = stim_blank;
acute_struct.stim_offset_nan = stim_offset_nan;
acute_struct.psth_normalized = 'zscore_omit_stim_offset';
acute_struct.stim_start = stim_start;
acute_struct.stim_end = stim_end;
acute_struct.call_onset_times = call_onset;
acute_struct.call_offset_times = call_offset;
acute_struct.n_calls_kept = numel(call_onset);
acute_struct.n_calls_total = numel(call_onset_all);
acute_struct.call_psth_onset = call_psth_onset;
acute_struct.call_psth_offset = call_psth_offset;
acute_struct.green_psth = green_psth;
acute_struct.green_onset_lfp_phase = green_onset_lfp_phase;
acute_struct.green_onset_breath_phase = green_onset_breath_phase;
acute_struct.mean_rate = mean_rate;
end

function [rate, t] = binned_rate(spikes, t0, t1, bin_w)
    if ~(isfinite(t0) && isfinite(t1) && t1 > t0)
        t = [];
        rate = [];
        return
    end
    edges = t0:bin_w:t1;
    if numel(edges) < 2
        t = [];
        rate = [];
        return
    end
    t = edges(2:end) - 0.5 * bin_w;
    rate = histcounts(spikes(:), edges) / bin_w;
end


function rate = nan_near_times(rate, t, t_mark, radius)
    if isempty(rate) || isempty(t_mark)
        return
    end
    t = t(:)';
    near = false(size(rate));
    for k = 1:numel(t_mark)
        if isfinite(t_mark(k))
            near = near | abs(t - t_mark(k)) < radius;
        end
    end
    rate(near) = NaN;
end


function z = nan_zscore(x)
    z = nan(size(x));
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');
    if ~(sd > 0) || isnan(sd)
        return
    end
    z = (x - mu) / sd;
end


function mu = mean_snippets(sig, t, events, rel_t)
    n_bins = numel(rel_t);
    if isempty(sig) || isempty(t) || isempty(events)
        mu = nan(1, n_bins);
        return
    end
    dt = t(2) - t(1);
    t0 = t(1);
    n_t = numel(t);
    X = nan(numel(events), n_bins);
    events = events(:);
    rel_t = rel_t(:)';
    for j = 1:numel(events)
        idx = round((events(j) + rel_t - t0) / dt) + 1;
        ok = idx >= 1 & idx <= n_t;
        row = nan(1, n_bins);
        row(ok) = sig(idx(ok));
        X(j, :) = row;
    end
    mu = mean(X, 1, 'omitnan');
end
