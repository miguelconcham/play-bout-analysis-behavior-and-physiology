function clusters = good_clusters_from_sorting(cluster_info, spike_times, spike_clusters, source_name)
%GOOD_CLUSTERS_FROM_SORTING  Good Phy/Kilosort units and their spike samples.
%
%   clusters = good_clusters_from_sorting(cluster_info, spike_times, spike_clusters)
%
% cluster_info is a table or struct with id (or cluster_id), ch, and group
% (or KSLabel). ch is 0-based, as in cluster_info.tsv in this repository.
% The returned channel is ch + 1, the row in the continuous int16 file.
% spike_times are the integers stored in spike_times.npy. A zero minimum is
% treated as a 0-based index and shifted by one sample.

if nargin < 4
    source_name = '';
end

ids = column_values(cluster_info, {'id', 'cluster_id'});
channels = column_values(cluster_info, {'ch', 'channel'});
groups = column_values(cluster_info, {'group', 'KSLabel'});
if numel(ids) ~= numel(channels) || numel(ids) ~= numel(groups)
    error('good_clusters_from_sorting:size', ...
        'id, channel, and group columns must have the same number of rows.');
end

spike_samples = spike_times_to_samples(spike_times);
spike_clusters = double(spike_clusters(:));
if numel(spike_samples) ~= numel(spike_clusters)
    error('good_clusters_from_sorting:size', ...
        'spike_times and spike_clusters must have the same number of spikes.');
end

clusters = struct('id', {}, 'channel', {}, 'spike_samples', {}, ...
    'n_spikes', {}, 'source', {});
for unit = 1:numel(ids)
    if ~strcmp(group_label(groups, unit), 'good')
        continue
    end
    unit_id = scalar_number(ids, unit);
    channel = scalar_number(channels, unit) + 1;
    if ~isfinite(unit_id) || ~isfinite(channel) || channel ~= round(channel) || channel < 1
        continue
    end
    samples = spike_samples(spike_clusters == unit_id);
    samples = samples(isfinite(samples) & samples == round(samples) & samples >= 1);
    if isempty(samples)
        continue
    end
    clusters(end + 1).id = unit_id; %#ok<AGROW>
    clusters(end).channel = channel;
    clusters(end).spike_samples = samples(:);
    clusters(end).n_spikes = numel(samples);
    clusters(end).source = source_name;
end
end

function spike_samples = spike_times_to_samples(spike_times)
spike_samples = double(spike_times(:));
if isempty(spike_samples)
    return
end
if any(~isfinite(spike_samples)) || any(spike_samples ~= round(spike_samples))
    error('good_clusters_from_sorting:spikeTime', ...
        'spike_times must be integer sample indices.');
end
if min(spike_samples) == 0
    spike_samples = spike_samples + 1;
elseif any(spike_samples < 0)
    error('good_clusters_from_sorting:spikeTime', ...
        'spike_times must be non-negative sample indices.');
end
spike_samples = round(spike_samples);
end

function values = column_values(cluster_info, names)
for name_index = 1:numel(names)
    name = names{name_index};
    if isstruct(cluster_info) && isfield(cluster_info, name)
        values = cluster_info.(name);
        return
    end
    if exist('istable', 'builtin') && istable(cluster_info) ...
            && any(strcmp(cluster_info.Properties.VariableNames, name))
        values = cluster_info.(name);
        return
    end
end
error('good_clusters_from_sorting:column', ...
    'cluster_info needs a column named %s.', strjoin(names, ' or '));
end

function label = group_label(groups, index)
if iscell(groups)
    value = groups{index};
elseif ischar(groups)
    value = groups(index, :);
elseif exist('isstring', 'builtin') && isstring(groups)
    value = char(groups(index));
elseif exist('iscategorical', 'builtin') && iscategorical(groups)
    value = char(groups(index));
else
    value = groups(index);
end
if exist('isstring', 'builtin') && isstring(value)
    value = char(value);
end
if isnumeric(value)
    label = '';
    return
end
label = lower(strtrim(value));
end

function number = scalar_number(values, index)
if iscell(values)
    number = double(values{index});
else
    number = double(values(index));
end
end
