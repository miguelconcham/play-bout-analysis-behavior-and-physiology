function parsing_struct = GENERATE_SPIKE_TRAIN_PARSING(npx_data_dir, Hd_freq)
% GENERATE_SPIKE_TRAIN_PARSING  Firing-run oscillation phase and power per neuron.
%
%   parsing_struct = GENERATE_SPIKE_TRAIN_PARSING(npx_data_dir)
%   parsing_struct = GENERATE_SPIKE_TRAIN_PARSING(npx_data_dir, Hd_freq)
%
% Loads channel map, PAG spikes, and PAG LFP from the session folder.
% For each good/mua neuron, detects firing runs and stores oscillation
% phase/power quantities (band set by Hd_freq):
%   - cycle-aligned PSTHs and sine fits
%   - phase entrainment: MVL and PPC vs circularly shifted spike trains
%   - LFP range / envelope / frequency and spike-triggered envelope
%
% phase_entrainment columns:
%   [MVL, MVL_p, MVL_p_maxshift, preferred_phase, PPC, PPC_p]
%
% Example:
%   parsing_struct = GENERATE_SPIKE_TRAIN_PARSING(npx_dir, Hd_freq);

if nargin < 1 || isempty(npx_data_dir)
    error('GENERATE_SPIKE_TRAIN_PARSING:npx_data_dir', 'Provide the NPX session folder.');
end

data_root        = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
area_limit_table = [data_root, '\Area_limits_GoodLooking.xlsx'];
area2analyse     = 'PAG';
area2analyze     = 'PAG';

sr_LFP = 2500;
if nargin < 2 || isempty(Hd_freq)
    Hd_freq = designfilt('bandpassfir', ...
        'FilterOrder', 2500, ...
        'CutoffFrequency1', 6, ...
        'CutoffFrequency2', 12, ...
        'SampleRate', sr_LFP);
end

% Fits, filtfilt, and related linear algebra emit RCOND / rank-deficient
% warnings on short runs; restore the caller's warning state on exit.
warning_state = warning;
cleanup_warn = onCleanup(@() warning(warning_state)); %#ok<NASGU>
warning('off');

animal_code        = strsplit(npx_data_dir, '\');
animal_code        = animal_code{end};
animal_code_params = strsplit(animal_code, ' ');
animal_batch       = animal_code_params{1};
repeated_animal    = animal_code_params{3};

%% Parsing parameters
n_rand            = 500;
osc_range         = [Hd_freq.CutoffFrequency1 Hd_freq.CutoffFrequency2];
bin_size          = 0.001;
smoothing_time    = 0.04;
smooth_window     = round(smoothing_time / bin_size);
max_coverage      = 0.9;
rate_wind         = 0.001;
norm_range        = [-3.5 3.5];
conv_function     = 500 * normpdf(linspace(norm_range(1), norm_range(2), 101));
conv_function     = conv_function / sum(normpdf(linspace(norm_range(1), norm_range(2), 101)));
hist_range        = bin_size * ceil(1 / (osc_range(1) * bin_size)) * [-1 1];
hist_edges        = hist_range(1):bin_size:hist_range(2);
hist_edges_centers = hist_edges(2:end) - bin_size / 2;
angle_hist_range  = [-pi pi];
angle_bin_size    = pi / 16;
angle_hist_edges  = angle_hist_range(1):angle_bin_size:angle_hist_range(2);
angle_hist_edges_centers = angle_hist_edges(2:end) - angle_bin_size / 2;
temporal_shift_range  = [-0.250 0.250];
temporal_shift_step   = 0.01;
temporal_shift_values = temporal_shift_range(1):temporal_shift_step:temporal_shift_range(2);
autocorr_range    = [-0.5 0.5];
autocorr_bin_size = 0.01;
autocorr_bin_edges = autocorr_range(1):autocorr_bin_size:autocorr_range(2);
treshold          = 2;
min_run_length    = 0.3;
min_run_spikes    = 4;
run_length_limits = [min_run_length 10];
n_psth_samples    = 50;

%% Determining NPX type and channel map
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
            parsing_struct = [];
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

%% Load spikes
spike_times    = double(readNPY([npx_data_dir, '\spike_times_', area2analyze, '.npy'])) / 30000;
spike_clusters = readNPY([npx_data_dir, '\spike_clusters_', area2analyze, '.npy']);
cluster_info   = readtable([npx_data_dir, '\cluster_info_', area2analyze, '.tsv'], 'FileType', 'text', 'Delimiter', '\t');
cluster_info   = cluster_info(ismember(cluster_info.group, {'mua', 'good'}), :);
these_neurons_areas = areas_by_channel(cluster_info.ch + 1);
cluster_info.area   = these_neurons_areas;
cluster_info.ch     = cluster_info.ch + 1;
channels_with_spikes = unique(cluster_info.ch)';

%% Load LFP
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

config = [];
config.n_rand = n_rand;
config.osc_range = osc_range;
config.bin_size = bin_size;
config.smoothing_time = smoothing_time;
config.smooth_window = smooth_window;
config.max_coverage = max_coverage;
config.rate_wind = rate_wind;
config.norm_range = norm_range;
config.conv_function = conv_function;
config.hist_range = hist_range;
config.hist_edges = hist_edges;
config.hist_edges_centers = hist_edges_centers;
config.angle_hist_range = angle_hist_range;
config.angle_bin_size = angle_bin_size;
config.angle_hist_edges = angle_hist_edges;
config.angle_hist_edges_centers = angle_hist_edges_centers;
config.autocorr_range = autocorr_range;
config.autocorr_bin_size = autocorr_bin_size;
config.autocorr_bin_edges = autocorr_bin_edges;
config.acg_edge_corrected = true;
config.treshold = treshold;
config.min_run_length = min_run_length;
config.min_run_spikes = min_run_spikes;
config.temporal_shift_range = temporal_shift_range;
config.temporal_shift_step = temporal_shift_step;
config.temporal_shift_values = temporal_shift_values;
config.n_psth_samples = n_psth_samples;
config.repeated_animal = repeated_animal;
config.phase_entrainment_fields = {'MVL', 'MVL_p', 'MVL_p_maxshift', 'preferred_phase', 'PPC', 'PPC_p'};

run_events_table_per_cell = cell(size(cluster_info, 1), 1);
assigned = false(size(cluster_info, 1), 1);

%% Per-channel LFP, then per-neuron firing-run parsing
n_channels = numel(channels_with_spikes);
disp(['Processing ', num2str(n_channels), ' channels'])
for ch_i = 1:n_channels
    ch_n = channels_with_spikes(ch_i);
    disp(['Channel ', num2str(ch_i), ' of ', num2str(n_channels), ' (ch #', num2str(ch_n), ')'])

    filtered_osc = filtfilt(Hd_freq.Coefficients, 1, PAG_LFP(ch_n, :));
    hilbert_osc  = hilbert(filtered_osc);
    osc_envelope = abs(hilbert_osc);
    uniform_phases = circular_uniformize(angle(hilbert_osc));

    peak_prominence = 0.1 * std(filtered_osc);
    peak_distance   = sr_LFP / osc_range(2);
    [~, min_locs]   = findpeaks(-filtered_osc, ...
        'MinPeakDistance', peak_distance, ...
        'MinPeakProminence', peak_prominence);

    neurons_this_channel = find(ismember(cluster_info.ch, ch_n))';
    disp([num2str(numel(neurons_this_channel)), ' neuron(s) in this channel'])
    disp([num2str(sum(assigned == 0)), ' neuron(s) to be analyzed'])

    for nn = neurons_this_channel
        this_id     = cluster_info.cluster_id(nn);
        this_spikes = spike_times(spike_clusters == this_id);
        this_spikes = this_spikes(this_spikes >= lfp_time(1) & this_spikes <= lfp_time(end));
        ct          = cluster_group_as_char(cluster_info.group, nn);
        this_area   = cluster_area_as_char(cluster_info.area, nn);

        this_neuron_phases = interp1(lfp_time, uniform_phases, this_spikes);
        valid_spk          = ~isnan(this_neuron_phases);
        this_neuron_phases = this_neuron_phases(valid_spk);
        this_neuron_spikes = this_spikes(valid_spk);
        this_osc_power     = interp1(lfp_time, osc_envelope, this_neuron_spikes);

        last_relevant_event = max([lfp_time(end); this_neuron_spikes]);
        all_rate = zeros(ceil(last_relevant_event / rate_wind), 1);
        spike_index = unique(ceil(this_neuron_spikes / rate_wind));
        spike_index(spike_index < 1 | spike_index > numel(all_rate)) = [];
        all_rate(spike_index) = 1;
        all_rate = conv(all_rate', conv_function);
        all_rate = all_rate(51:end - 50);

        run_events_table = find_runs(all_rate', treshold, rate_wind, this_neuron_spikes, min_run_spikes, run_length_limits, max_coverage);

        if isempty(run_events_table) || all(isnan(run_events_table.RunStartTime))
            run_events_table_per_cell{nn} = run_events_table;
            assigned(nn) = true;
            continue
        end

        n_runs = size(run_events_table, 1);
        polar_histogram               = nan(n_runs, numel(angle_hist_edges_centers));
        mean_osc_alignment            = nan(n_runs, numel(hist_edges_centers));
        real_osc_alignment            = nan(n_runs, numel(hist_edges_centers));
        mean_raw_osc_alignment        = nan(n_runs, numel(hist_edges_centers));
        real_raw_osc_alignment        = nan(n_runs, numel(hist_edges_centers));
        fitted_osc_alignment          = nan(n_runs, numel(hist_edges_centers));
        sin_fit_param                 = nan(n_runs, 6);
        sin_fit_param_no_trend        = nan(n_runs, 6);
        sin_fit_param_no_mean         = nan(n_runs, 6);
        lfp_param                     = nan(n_runs, 3);
        spike_osc_power               = nan(n_runs, 1);
        this_neuron_autoccorrelograms = nan(n_runs, numel(autocorr_bin_edges) - 1);
        phase_entrainment             = nan(n_runs, 6);
        temporal_shifts               = nan(n_runs, numel(temporal_shift_values));

        for j = 1:n_runs
            run_start = run_events_table.RunStartTime(j);
            run_end   = run_events_table.RunEndTime(j);
            if isnan(run_start) || isnan(run_end)
                continue
            end

            this_run_peaks = min_locs(lfp_time(min_locs) >= run_start & lfp_time(min_locs) <= run_end);
            this_run_peaks = lfp_time(this_run_peaks);
            this_run_lfp   = filtered_osc(lfp_time >= run_start & lfp_time <= run_end);
            this_run_env   = osc_envelope(lfp_time >= run_start & lfp_time <= run_end);
            if isempty(this_run_lfp)
                lfp_range = NaN;
                mean_envelope = NaN;
                mean_freq = NaN;
            else
                lfp_range = range(this_run_lfp) / 2;
                mean_envelope = median(this_run_env);
                if numel(this_run_peaks) > 1
                    mean_freq = 1 / mean(diff(this_run_peaks));
                else
                    mean_freq = NaN;
                end
            end
            lfp_param(j, :) = [lfp_range mean_envelope mean_freq];

            this_run_spikes = this_neuron_spikes(this_neuron_spikes >= run_start & this_neuron_spikes <= run_end);
            run_spk         = this_neuron_spikes >= run_start & this_neuron_spikes <= run_end;
            spike_osc_power(j) = mean(this_osc_power(run_spk), 'omitnan');
            all_lags = autocorrelogram(this_run_spikes, autocorr_range);
            acg_counts = histcounts(all_lags, autocorr_bin_edges);
            this_neuron_autoccorrelograms(j, :) = edge_correct_acg( ...
                acg_counts, autocorr_bin_edges, run_end - run_start);

            this_run_psth         = nan(numel(this_run_peaks), numel(hist_edges_centers));
            this_run_psth_samples = min(n_psth_samples, 2^numel(this_run_peaks));
            mean_responses        = nan(this_run_psth_samples, numel(hist_edges_centers));
            raw_responses         = nan(this_run_psth_samples, numel(hist_edges_centers));

            for nps = 1:this_run_psth_samples
                spikes_this_psth = this_neuron_spikes;
                if numel(this_run_peaks) > 0
                    order2collect = randsample(numel(this_run_peaks), numel(this_run_peaks));
                    for peak_n = order2collect'
                        peak_time = this_run_peaks(peak_n);
                        in_win = spikes_this_psth >= peak_time + hist_range(1) & spikes_this_psth <= peak_time + hist_range(2);
                        spikes2plot = spikes_this_psth(in_win) - peak_time;
                        spikes_this_psth(in_win) = [];
                        this_run_psth(peak_n, :) = histcounts(spikes2plot, hist_edges);
                    end
                end
                if size(this_run_psth, 1) > 1
                    mean_response = movmean(mean(this_run_psth), smooth_window) / bin_size;
                    raw_response  = mean(this_run_psth) / bin_size;
                elseif size(this_run_psth, 1) == 1
                    mean_response = movmean(this_run_psth, smooth_window) / bin_size;
                    raw_response  = this_run_psth / bin_size;
                else
                    mean_response = nan(1, numel(hist_edges_centers));
                    raw_response  = nan(1, numel(hist_edges_centers));
                end
                mean_responses(nps, :) = mean_response;
                raw_responses(nps, :)  = raw_response;
            end

            if this_run_psth_samples > 1
                mean_response = mean(mean_responses, 'omitnan');
                raw_response  = mean(raw_responses, 'omitnan');
            else
                mean_response = mean_responses;
                raw_response  = raw_responses;
            end

            real_index = find(~any(isnan(mean_responses), 2), 1, 'first');
            if ~isempty(real_index)
                real_osc_values = mean_responses(real_index, :);
                real_osc_alignment(j, :)     = mean_responses(real_index, :);
                real_raw_osc_alignment(j, :) = raw_responses(real_index, :);
                y = real_osc_values(ceil(0.5 * smooth_window):end - ceil(0.5 * smooth_window));
                rate_range = range(y);
                min_rate = min(y);
                if ~any(isnan(y)) && ~isempty(y) && rate_range > 0
                    y = (y - min_rate) / rate_range;
                    x = hist_edges_centers(ceil(0.5 * smooth_window):end - ceil(0.5 * smooth_window));
                    [~, peaks_mean_positive] = findpeaks(y, 'MinPeakDistance', (1 / osc_range(2)) / bin_size, 'MinPeakProminence', std(y));
                    [~, peaks_mean_negative] = findpeaks(-y, 'MinPeakDistance', (1 / osc_range(2)) / bin_size, 'MinPeakProminence', std(y));
                    consecutive_peaks = sort([peaks_mean_positive, peaks_mean_negative]);
                    per = 2 * mean(diff(consecutive_peaks)) * bin_size;
                    if numel(consecutive_peaks) < 2 || per == 0
                        per = 1 / mean(osc_range);
                    end
                    mdl = fittype('a*sin(b*x + c) + d*x + e', 'indep', 'x');
                    [fittedmdl, gof] = fit(x', y', mdl, 'start', [rand(), 1 / (per / (2 * pi)), rand(), rand(), rand()]);
                    sin_fit_param_no_mean(j, 1) = fittedmdl.a * rate_range;
                    sin_fit_param_no_mean(j, 2) = fittedmdl.b;
                    sin_fit_param_no_mean(j, 3) = fittedmdl.c;
                    sin_fit_param_no_mean(j, 4) = fittedmdl.d * rate_range;
                    sin_fit_param_no_mean(j, 5) = fittedmdl.e * rate_range + min_rate;
                    sin_fit_param_no_mean(j, 6) = gof.rsquare;
                end
            end

            y = mean_response(ceil(0.5 * smooth_window):end - ceil(0.5 * smooth_window));
            rate_range = range(y);
            min_rate = min(y);
            if ~any(isnan(y)) && ~isempty(y) && rate_range > 0
                y = (y - min_rate) / rate_range;
                x = hist_edges_centers(ceil(0.5 * smooth_window):end - ceil(0.5 * smooth_window));
                [~, peaks_mean_positive] = findpeaks(y, 'MinPeakDistance', (1 / osc_range(2)) / bin_size, 'MinPeakProminence', std(y));
                [~, peaks_mean_negative] = findpeaks(-y, 'MinPeakDistance', (1 / osc_range(2)) / bin_size, 'MinPeakProminence', std(y));
                consecutive_peaks = sort([peaks_mean_positive, peaks_mean_negative]);
                per = 2 * mean(diff(consecutive_peaks)) * bin_size;
                if numel(consecutive_peaks) < 2 || per == 0
                    per = 1 / mean(osc_range);
                end
                mdl = fittype('a*sin(b*x + c) + d*x + e', 'indep', 'x');
                [fittedmdl, gof] = fit(x', y', mdl, 'start', [rand(), 1 / (per / (2 * pi)), rand(), rand(), rand()]);
                sin_fit_param(j, 1) = fittedmdl.a * rate_range;
                sin_fit_param(j, 2) = fittedmdl.b;
                sin_fit_param(j, 3) = fittedmdl.c;
                sin_fit_param(j, 4) = fittedmdl.d * rate_range;
                sin_fit_param(j, 5) = fittedmdl.e * rate_range + min_rate;
                sin_fit_param(j, 6) = gof.rsquare;
                fitted_osc_alignment(j, ceil(0.5 * smooth_window):end - ceil(0.5 * smooth_window)) = fittedmdl(x) * rate_range + min_rate;

                [fittedmdl, gof] = fit_sin(x', y', [rand(), 1 / (per / (2 * pi)), rand(), rand()]);
                sin_fit_param_no_trend(j, 1) = fittedmdl.a * rate_range;
                sin_fit_param_no_trend(j, 2) = fittedmdl.b;
                sin_fit_param_no_trend(j, 3) = fittedmdl.c;
                sin_fit_param_no_trend(j, 4) = fittedmdl.d * rate_range + min_rate;
                sin_fit_param_no_trend(j, 5) = gof.rsquare;
            end

            mean_osc_alignment(j, ceil(0.5 * smooth_window):(end - ceil(0.5 * smooth_window))) = ...
                mean_response(ceil(0.5 * smooth_window):end - ceil(0.5 * smooth_window));
            mean_raw_osc_alignment(j, :) = raw_response;

            this_phases = this_neuron_phases(this_neuron_spikes >= run_start & this_neuron_spikes <= run_end);
            time_shift_index = 1;
            for time_shift = temporal_shift_values
                shifted_times = this_run_spikes + time_shift;
                if ~isempty(shifted_times)
                    shifted_phases = interp1(lfp_time, uniform_phases, shifted_times);
                    shifted_phases = shifted_phases(~isnan(shifted_phases));
                    if ~isempty(shifted_phases)
                        temporal_shifts(j, time_shift_index) = circ_r(shifted_phases);
                    end
                end
                time_shift_index = time_shift_index + 1;
            end

            if ~isempty(this_phases)
                r = circ_r(this_phases);
                phase = circ_mean(this_phases);
                this_ppc = compute_ppc_fast(this_phases);
            else
                r = NaN;
                phase = NaN;
                this_ppc = NaN;
            end
            if numel(this_run_spikes) > 1
                shifted_spiketrain = shift_spiketrain(this_run_spikes, run_start, run_end, rand(n_rand, 1) * (run_end - run_start));
            else
                shifted_spiketrain = [];
            end
            if ~isempty(shifted_spiketrain)
                this_neuron_random_phases = interp1(lfp_time, uniform_phases, shifted_spiketrain);
                valid_rand = ~any(isnan(this_neuron_random_phases), 2);
                if any(valid_rand) && ~isnan(r)
                    rand_mvl = circ_r(this_neuron_random_phases(valid_rand, :)');
                    rand_mvl = sort(rand_mvl);
                    [~, loc] = min(abs(r - rand_mvl));
                    r_stat = mean(loc) / n_rand;
                    [~, loc] = min(abs(max(temporal_shift_values) - rand_mvl));
                    r_stat_max = mean(loc) / n_rand;

                    rand_ppc = nan(sum(valid_rand), 1);
                    rand_phase_rows = this_neuron_random_phases(valid_rand, :);
                    for nr = 1:size(rand_phase_rows, 1)
                        rand_ppc(nr) = compute_ppc_fast(rand_phase_rows(nr, :)');
                    end
                    rand_ppc = sort(rand_ppc);
                    [~, loc] = min(abs(this_ppc - rand_ppc));
                    ppc_stat = mean(loc) / n_rand;

                    phase_entrainment(j, :) = [r r_stat r_stat_max phase this_ppc ppc_stat];
                end
            end
            polar_histogram(j, :) = histcounts(this_phases, angle_hist_edges);
        end

        run_events_table.SinFit                 = sin_fit_param;
        run_events_table.SinFit_NOTREND         = sin_fit_param_no_trend;
        run_events_table.SinFit_NOMEAN          = sin_fit_param_no_mean;
        run_events_table.OscPsth                = mean_osc_alignment;
        run_events_table.OscPsthReal            = real_osc_alignment;
        run_events_table.RowOscPsth             = mean_raw_osc_alignment;
        run_events_table.RowOscPsthReal         = real_raw_osc_alignment;
        run_events_table.autoccorrelograms      = this_neuron_autoccorrelograms;
        run_events_table.fitted_osc_alignment   = fitted_osc_alignment;
        run_events_table.lfp_param              = lfp_param;
        run_events_table.spike_osc_power        = spike_osc_power;
        run_events_table.phase_entrainment      = phase_entrainment;
        run_events_table.temporal_shifts        = temporal_shifts;
        run_events_table.polar_histogram        = polar_histogram;
        run_events_table.Id                     = ones(n_runs, 1) * this_id;
        run_events_table.Ch                     = ones(n_runs, 1) * ch_n;
        run_events_table.ct                     = repmat({ct}, n_runs, 1);
        run_events_table.area                   = repmat({this_area}, n_runs, 1);

        run_events_table_per_cell{nn} = run_events_table;
        assigned(nn) = true;
    end
end

if any(~assigned)
    disp('Not all neurons assigned')
end

ALL_TABLES = table();
for j = 1:numel(run_events_table_per_cell)
    trans_table = run_events_table_per_cell{j};
    if isempty(trans_table) || all(isnan(trans_table.RunStartTime))
        continue
    end
    trans_names = trans_table.Properties.VariableNames;
    trans_names(ismember(trans_names, 'RunMaxate')) = {'RunMaxRate'};
    trans_table.Properties.VariableNames = trans_names;
    ALL_TABLES = [ALL_TABLES; trans_table]; %#ok<AGROW>
end

parsing_struct = [];
parsing_struct.cluster_info              = cluster_info;
parsing_struct.run_events_table_per_cell = run_events_table_per_cell;
parsing_struct.ALL_TABLES                = ALL_TABLES;
parsing_struct.config                    = config;
parsing_struct.areas_by_channel          = areas_by_channel;
parsing_struct.channel_map               = channel_map;
parsing_struct.animal_code               = animal_code;
parsing_struct.NPX_Type                  = NPX_Type;
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
    this_area = this_area{1};
end
this_area = char(string(this_area));
end
