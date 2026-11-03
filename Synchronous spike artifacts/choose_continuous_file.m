function [continuous_file, checks] = choose_continuous_file(candidates, ...
    spike_samples, neuron_channel, channel_range, choice)
%CHOOSE_CONTINUOUS_FILE  Pick the continuous AP file that contains the spike.
%
% candidates are paths from NPX_raw_data_folders.xlsx for one session. A
% directory is read as continuous.dat inside it. When several paths match,
% each file is scored with the mean spike on the neuron channel (at most
% n_check spikes, spread through the train). The largest mean is kept.
% A flat mean, below min_spike_scale counts after baseline subtraction,
% means the file is the wrong recording and the session should stop.

if nargin < 5 || isempty(choice)
    choice = struct();
end
choice = fill_choice(choice);
checks = struct('path', {}, 'scale', {}, 'message', {});
paths = candidate_paths(candidates);
if isempty(paths)
    error('choose_continuous_file:none', ...
        'No continuous-file path was provided for this session.');
end

best_scale = -inf;
best_index = 0;
for path_index = 1:numel(paths)
    file_path = paths{path_index};
    checks(end + 1).path = file_path; %#ok<AGROW>
    checks(end).scale = NaN;
    checks(end).message = '';
    try
        scale = neuron_scale(file_path, spike_samples, neuron_channel, ...
            channel_range, choice);
        checks(end).scale = scale;
        if scale > best_scale
            best_scale = scale;
            best_index = path_index;
        end
    catch read_error
        checks(end).message = read_error.message;
        if numel(paths) == 1
            rethrow(read_error);
        end
    end
end

if best_index == 0
    detail = '';
    for check_index = 1:numel(checks)
        detail = sprintf('%s\n  %s  %s', detail, checks(check_index).path, ...
            checks(check_index).message);
    end
    error('choose_continuous_file:read', ...
        'None of the continuous files could be read.%s', detail);
end
if ~(best_scale >= choice.min_scale)
    error('choose_continuous_file:flat', ...
        ['The continuous file does not show a spike on channel %g ' ...
         '(best mean scale %.3g, minimum %.3g).'], ...
        neuron_channel, best_scale, choice.min_scale);
end
continuous_file = paths{best_index};
end

function scale = neuron_scale(file_path, spike_samples, neuron_channel, ...
    channel_range, choice)
picked = spread_samples(spike_samples, choice.n_check);
[waveforms, channels_kept, time_samples] = load_waveforms_from_continuous( ...
    file_path, picked, neuron_channel, channel_range, ...
    'NChannels', choice.n_channels, ...
    'SamplesBefore', choice.samples_before, ...
    'SamplesAfter', choice.samples_after);
neuron_row = find(channels_kept == neuron_channel, 1);
if isempty(neuron_row) || isempty(waveforms)
    scale = 0;
    return
end
centered = center_waveform_snippets(waveforms, time_samples);
traces = reshape(centered(:, neuron_row, :), size(centered, 1), []);
finite_samples = isfinite(traces);
traces(~finite_samples) = 0;
counts = sum(finite_samples, 1);
if ~any(counts)
    scale = 0;
    return
end
mean_trace = sum(traces, 1) ./ max(counts, 1);
mean_trace(counts == 0) = NaN;
scale = max(abs(mean_trace));
if isempty(scale) || ~isfinite(scale)
    scale = 0;
end
end

function picked = spread_samples(spike_samples, n_check)
spike_samples = round(double(spike_samples(:)));
if isempty(spike_samples)
    error('choose_continuous_file:spikes', ...
        'No spike samples were available to check the continuous file.');
end
n_check = min(n_check, numel(spike_samples));
if n_check == numel(spike_samples)
    picked = spike_samples;
    return
end
index = unique(round(linspace(1, numel(spike_samples), n_check)));
picked = spike_samples(index);
end

function paths = candidate_paths(candidates)
if ischar(candidates)
    candidates = {candidates};
elseif exist('isstring', 'builtin') && isstring(candidates)
    candidates = cellstr(candidates);
end
paths = {};
for index = 1:numel(candidates)
    path = candidates{index};
    if exist('isstring', 'builtin') && isstring(path)
        path = char(path);
    end
    if iscell(path)
        continue
    end
    path = strtrim(path);
    if isempty(path)
        continue
    end
    if exist(path, 'dir') == 7
        path = fullfile(path, 'continuous.dat');
    end
    if any(strcmp(paths, path))
        continue
    end
    paths{end + 1} = path; %#ok<AGROW>
end
end

function choice = fill_choice(choice)
defaults = sync_analysis_options(struct());
choice = fill_one(choice, 'min_scale', defaults.min_spike_scale);
choice = fill_one(choice, 'n_channels', defaults.n_ap_channels);
choice = fill_one(choice, 'n_check', defaults.n_check_spikes);
choice = fill_one(choice, 'samples_before', defaults.samples_before);
choice = fill_one(choice, 'samples_after', defaults.samples_after);
end

function choice = fill_one(choice, name, value)
if ~isfield(choice, name) || isempty(choice.(name))
    choice.(name) = value;
end
end
