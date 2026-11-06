function session = GENERATE_ACG_DELTA_COMPONENT_PHASE(npx_data_dir, Hd_freq, lag_range)
%GENERATE_ACG_DELTA_COMPONENT_PHASE  Phase of spikes in one ACG side peak.
%
%   session = GENERATE_ACG_DELTA_COMPONENT_PHASE(npx_data_dir)
%   session = GENERATE_ACG_DELTA_COMPONENT_PHASE(npx_data_dir, Hd_freq, lag_range)
%
% For each good/mua neuron, firing runs are the same rate runs used by
% GENERATE_SPIKE_TRAIN_PARSING. Inside each run, a pair whose lag falls in
% lag_range (default [0.2 0.5] s, the shaded delta side peak) contributes
% both spikes once. The symmetric autocorrelogram would list that pair
% twice; this function does not. Spikes between runs are in neither sample.
%
% Component spikes are the unique members of those pairs. Complement spikes
% are the other spikes in the same runs. Both samples are summarized on the
% circularly uniformized Hilbert phase of Hd_freq (default delta, 1-5 Hz)
% on the neuron's channel. Trough is +/- pi and peak is 0 in that coordinate,
% matching delta_lock_group. The lock label itself is applied in Analyze.
%
% Example:
%   session = GENERATE_ACG_DELTA_COMPONENT_PHASE(npx_dir, Hd_freq, [0.2 0.5]);

if nargin < 1 || isempty(npx_data_dir)
    error('GENERATE_ACG_DELTA_COMPONENT_PHASE:npx_data_dir', ...
        'Provide the NPX session folder.');
end

data_root        = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
area_limit_table = [data_root, '\Area_limits_GoodLooking.xlsx'];
area2analyse     = 'PAG';
area2analyze     = 'PAG';
sr_LFP           = 2500;

if nargin < 2 || isempty(Hd_freq)
    Hd_freq = designfilt('bandpassfir', ...
        'FilterOrder', 2000, ...
        'CutoffFrequency1', 1, ...
        'CutoffFrequency2', 5, ...
        'SampleRate', sr_LFP, ...
        'DesignMethod', 'window', ...
        'Window', 'hamming');
end
if nargin < 3 || isempty(lag_range)
    lag_range = [0.2 0.5];
end
lag_range = double(lag_range(:)');
if numel(lag_range) ~= 2 || ~all(isfinite(lag_range)) || lag_range(1) < 0 || ~(lag_range(2) > lag_range(1))
    error('GENERATE_ACG_DELTA_COMPONENT_PHASE:lag', ...
        'lag_range must be [min max] seconds, with 0 <= min < max.');
end

% Same firing-run parameters as GENERATE_SPIKE_TRAIN_PARSING. The side
% peak is defined on those runs, so a different detector would select a
% different set of spikes from the autocorrelograms.
bin_size          = 0.001;
norm_range        = [-3.5 3.5];
conv_function     = 500 * normpdf(linspace(norm_range(1), norm_range(2), 101));
conv_function     = conv_function / sum(normpdf(linspace(norm_range(1), norm_range(2), 101)));
rate_wind         = 0.001;
treshold          = 2;
min_run_length    = 0.3;
min_run_spikes    = 4;
max_coverage      = 0.9;
run_length_limits = [min_run_length 10];
n_phase_bins      = 32;

parts = strsplit(char(npx_data_dir), {'\', '/'});
animal_code = parts{end};
if isempty(animal_code) && numel(parts) >= 2
    animal_code = parts{end - 1};
end
animal_code_params = strsplit(animal_code, ' ');
if numel(animal_code_params) < 3
    error('GENERATE_ACG_DELTA_COMPONENT_PHASE:animal', ...
        'Session folder "%s" does not look like "B# ... Single/Dual ...".', animal_code);
end
animal_batch    = animal_code_params{1};
repeated_animal = animal_code_params{3};

%% Channel map
hard_coded_x_coords_NPX2 = [8 40; 258 290; 508 540; 758 790];
load([npx_data_dir, '\', 'chann_map_', area2analyse, '.mat'], 'chanMap', 'xcoords', 'ycoords')
if any(ismember(xcoords, hard_coded_x_coords_NPX2))
    NPX_Type = 2;
else
    NPX_Type = 1;
    if ~ismember(192, chanMap)
        pos_191 = find(chanMap == 191);
        pos_193 = find(chanMap == 193);
        if pos_193 == pos_191 + 1
            xcoords = [xcoords; NaN];
            xcoords(pos_193 + 1:end) = xcoords(pos_193:end - 1);
            xcoords(pos_193) = 43;
            ycoords = [ycoords; NaN];
            ycoords(pos_193 + 1:end) = ycoords(pos_193:end - 1);
            ycoords(pos_193) = 1900;
            chanMap = [chanMap; NaN];
            chanMap(pos_193 + 1:end) = chanMap(pos_193:end - 1);
            chanMap(pos_193) = 192;
        else
            disp('Inconsistent ChannelMap')
            session = [];
            return
        end
    end
end

disp('Loading Channel Map')
areas_by_channel = cell(384, 1);
channel_map      = nan(384, 2);
area_limit       = readtable(area_limit_table);
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName, this_animal), :);

if NPX_Type == 1
    for ch_n = 1:384
        ch = chanMap(ch_n);
        channel_map(ch, 1) = xcoords(ch_n);
        channel_map(ch, 2) = ycoords(ch_n);
        areas_by_channel{ch} = area_limit.area{ycoords(ch_n) >= area_limit.depth_start & ycoords(ch_n) < area_limit.depth_end + 1 & ismember(area_limit.Probe_Area, area2analyse)};
    end
else
    for ch_n = 1:384
        probe_n = find(any(ismember(hard_coded_x_coords_NPX2, xcoords(ch_n)), 2));
        ch = chanMap(ch_n);
        channel_map(ch, 1) = xcoords(ch_n);
        channel_map(ch, 2) = ycoords(ch_n);
        areas_by_channel{ch} = area_limit.area{ycoords(ch_n) >= area_limit.depth_start & ycoords(ch_n) < area_limit.depth_end + 1 & area_limit.ProbeNum == probe_n & ismember(area_limit.Probe_Area, 'PAG')};
    end
end

%% Spikes and LFP
spike_times    = double(readNPY([npx_data_dir, '\spike_times_', area2analyze, '.npy'])) / 30000;
spike_clusters = readNPY([npx_data_dir, '\spike_clusters_', area2analyze, '.npy']);
cluster_info   = readtable([npx_data_dir, '\cluster_info_', area2analyze, '.tsv'], 'FileType', 'text', 'Delimiter', '\t');
cluster_info   = cluster_info(ismember(cluster_info.group, {'mua', 'good'}), :);
if isempty(cluster_info)
    session = [];
    return
end
these_neurons_areas = areas_by_channel(cluster_info.ch + 1);
cluster_info.area   = these_neurons_areas;
cluster_info.ch     = cluster_info.ch + 1;
channels_with_spikes = unique(cluster_info.ch)';
channels_with_spikes = channels_with_spikes(isfinite(channels_with_spikes));

disp('LOADING LFP')
if NPX_Type == 2
    load([npx_data_dir, '\', 'LFP_', area2analyse, '.mat'], 'LFP')
elseif NPX_Type == 1
    file_pointer = fopen([npx_data_dir, '\', 'LFP_', area2analyse, '.dat'], 'r');
    LFP = fread(file_pointer, 'int16');
    LFP = reshape(LFP, 384, numel(LFP) / 384);
    fclose(file_pointer);
end
disp('LFP LOADED')
PAG_LFP = double(LFP);
clear LFP
lfp_time = (1:size(PAG_LFP, 2)) / sr_LFP;

n_neu = height(cluster_info);
phase_probe = component_phase_summary([], n_phase_bins);
phase_edges = phase_probe.edges;
cluster_id = nan(n_neu, 1);
channel = nan(n_neu, 1);
area = repmat({''}, n_neu, 1);
ks_group = repmat({''}, n_neu, 1);
n_session_spikes = zeros(n_neu, 1);
n_run_spikes = zeros(n_neu, 1);
n_component = zeros(n_neu, 1);
n_complement = zeros(n_neu, 1);
n_pairs = zeros(n_neu, 1);
n_acg_entries = zeros(n_neu, 1);
component_angle = nan(n_neu, 1);
component_mvl = nan(n_neu, 1);
component_ppc = nan(n_neu, 1);
component_proj_peak = nan(n_neu, 1);
component_proj_trough = nan(n_neu, 1);
component_p_peak = nan(n_neu, 1);
component_p_trough = nan(n_neu, 1);
component_frac_peak = nan(n_neu, 1);
complement_angle = nan(n_neu, 1);
complement_mvl = nan(n_neu, 1);
complement_ppc = nan(n_neu, 1);
complement_proj_peak = nan(n_neu, 1);
complement_proj_trough = nan(n_neu, 1);
complement_p_peak = nan(n_neu, 1);
complement_p_trough = nan(n_neu, 1);
complement_frac_peak = nan(n_neu, 1);
component_hist = zeros(n_neu, n_phase_bins);
complement_hist = zeros(n_neu, n_phase_bins);

n_channels = numel(channels_with_spikes);
disp(['Processing ', num2str(n_channels), ' channels'])
for ch_i = 1:n_channels
    ch_n = channels_with_spikes(ch_i);
    if ch_n < 1 || ch_n > size(PAG_LFP, 1)
        continue
    end
    disp(['Channel ', num2str(ch_i), ' of ', num2str(n_channels), ' (ch #', num2str(ch_n), ')'])
    filtered_osc = filtfilt(Hd_freq.Coefficients, 1, PAG_LFP(ch_n, :));
    uniform_phases = circular_uniformize(angle(hilbert(filtered_osc)));

    neurons_this_channel = find(ismember(cluster_info.ch, ch_n))';
    for nn = neurons_this_channel
        this_id = cluster_info.cluster_id(nn);
        this_spikes = spike_times(spike_clusters == this_id);
        this_spikes = this_spikes(this_spikes >= lfp_time(1) & this_spikes <= lfp_time(end));
        this_phases = interp1(lfp_time, uniform_phases, this_spikes);
        valid_spk = isfinite(this_phases) & isfinite(this_spikes);
        this_spikes = this_spikes(valid_spk);
        this_phases = this_phases(valid_spk);

        cluster_id(nn) = double(this_id);
        channel(nn) = ch_n;
        area{nn} = cluster_area_as_char(cluster_info.area, nn);
        ks_group{nn} = cluster_group_as_char(cluster_info.group, nn);
        n_session_spikes(nn) = numel(this_spikes);
        if numel(this_spikes) < 2
            continue
        end

        runs = neuron_firing_runs(this_spikes, conv_function, rate_wind, ...
            treshold, min_run_spikes, run_length_limits, max_coverage, lfp_time(end));
        if isempty(runs)
            selected = acg_delta_component_spikes(this_spikes, [], [], lag_range);
        else
            selected = acg_delta_component_spikes( ...
                this_spikes, runs.RunStartTime, runs.RunEndTime, lag_range);
        end
        component_stats = component_phase_summary(this_phases(selected.mask), n_phase_bins);
        complement_stats = component_phase_summary( ...
            this_phases(selected.in_run & ~selected.mask), n_phase_bins);

        n_run_spikes(nn) = selected.n_run_spikes;
        n_component(nn) = selected.n_spikes;
        n_complement(nn) = complement_stats.n;
        n_pairs(nn) = selected.n_pairs;
        n_acg_entries(nn) = selected.n_acg_entries;
        component_angle(nn) = component_stats.angle;
        component_mvl(nn) = component_stats.mvl;
        component_ppc(nn) = component_stats.ppc;
        component_proj_peak(nn) = component_stats.proj_peak;
        component_proj_trough(nn) = component_stats.proj_trough;
        component_p_peak(nn) = component_stats.p_peak;
        component_p_trough(nn) = component_stats.p_trough;
        component_frac_peak(nn) = component_stats.frac_peak;
        complement_angle(nn) = complement_stats.angle;
        complement_mvl(nn) = complement_stats.mvl;
        complement_ppc(nn) = complement_stats.ppc;
        complement_proj_peak(nn) = complement_stats.proj_peak;
        complement_proj_trough(nn) = complement_stats.proj_trough;
        complement_p_peak(nn) = complement_stats.p_peak;
        complement_p_trough(nn) = complement_stats.p_trough;
        complement_frac_peak(nn) = complement_stats.frac_peak;
        component_hist(nn, :) = component_stats.hist;
        complement_hist(nn, :) = complement_stats.hist;
    end
end

neurons = table(cluster_id, channel, area, ks_group, n_session_spikes, n_run_spikes, ...
    n_component, n_complement, n_pairs, n_acg_entries, ...
    component_angle, component_mvl, component_ppc, ...
    component_proj_peak, component_proj_trough, component_p_peak, component_p_trough, component_frac_peak, ...
    complement_angle, complement_mvl, complement_ppc, ...
    complement_proj_peak, complement_proj_trough, complement_p_peak, complement_p_trough, complement_frac_peak);

session = struct();
session.animal_code = animal_code;
session.lag_range = lag_range;
session.osc_range = [Hd_freq.CutoffFrequency1, Hd_freq.CutoffFrequency2];
session.phase_edges = phase_edges;
session.neurons = neurons;
session.component_hist = component_hist;
session.complement_hist = complement_hist;
session.channel_map = channel_map;
session.NPX_Type = NPX_Type;
session.config = struct( ...
    'bin_size', bin_size, ...
    'rate_wind', rate_wind, ...
    'treshold', treshold, ...
    'min_run_length', min_run_length, ...
    'min_run_spikes', min_run_spikes, ...
    'max_coverage', max_coverage, ...
    'run_length_limits', run_length_limits, ...
    'n_phase_bins', n_phase_bins, ...
    'sr_LFP', sr_LFP);
end


function runs = neuron_firing_runs(this_spikes, conv_function, rate_wind, ...
    treshold, min_run_spikes, run_length_limits, max_coverage, lfp_end)
last_relevant_event = max([lfp_end; this_spikes(:)]);
all_rate = zeros(ceil(last_relevant_event / rate_wind), 1);
spike_index = unique(ceil(this_spikes / rate_wind));
spike_index(spike_index < 1 | spike_index > numel(all_rate)) = [];
all_rate(spike_index) = 1;
all_rate = conv(all_rate', conv_function);
all_rate = all_rate(51:end - 50);
runs = find_runs(all_rate', treshold, rate_wind, this_spikes, ...
    min_run_spikes, run_length_limits, max_coverage);
if isempty(runs) || all(isnan(runs.RunStartTime))
    runs = [];
end
end


function ct = cluster_group_as_char(group_col, nn)
ct = group_col(nn, :);
if iscell(ct)
    ct = ct{1};
end
ct = char(string(ct));
end


function this_area = cluster_area_as_char(area_col, nn)
this_area = area_col{nn};
if iscell(this_area)
    if isempty(this_area)
        this_area = '';
    else
        this_area = this_area{1};
    end
end
this_area = char(string(this_area));
end
