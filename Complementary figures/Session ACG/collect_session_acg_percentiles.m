function records = collect_session_acg_percentiles(save_root)
%COLLECT_SESSION_ACG_PERCENTILES  Percentile curve saved on each neuron.
%
% Reads neuron_*.mat under save_root, one session folder at a time. The
% file is loaded in full because the null matrices share that variable.
% Only the lock group and the score percentile are kept. A duplicate
% session and cluster id is skipped. Spike times are not in these files.

if exist(save_root, 'dir') ~= 7
    error('session_acg:percentile', ...
        'No session ACG folder at %s.', save_root);
end
files = dir(fullfile(save_root, '*', 'neuron_*.mat'));
files = files(~[files.isdir]);
if isempty(files)
    records = empty_records();
    return
end
paths = fullfile({files.folder}, {files.name});
[~, order] = sort(paths);
files = files(order);

curve_names = session_acg_names('curves');
records = empty_records();
keys = cell(0, 1);
n_files = numel(files);
for f = 1:n_files
    [~, base_name, extension] = fileparts(files(f).name);
    if ~strcmp(extension, '.mat') || isempty(regexp(base_name, '^neuron_\d+$', 'once'))
        continue
    end
    neuron_file = fullfile(files(f).folder, files(f).name);
    try
        loaded = load(neuron_file, 'neuron');
    catch failure
        session_acg_log('', sprintf('Skipping %s: %s', neuron_file, failure.message));
        continue
    end
    if ~isfield(loaded, 'neuron') || ~isstruct(loaded.neuron) || numel(loaded.neuron) ~= 1
        session_acg_log('', sprintf('Skipping %s: no neuron struct.', neuron_file));
        clear loaded
        continue
    end
    neuron = loaded.neuron;
    clear loaded
    if ~isfield(neuron, 'session') || ~isfield(neuron, 'cluster_id') || ~isfield(neuron, 'lock_group')
        session_acg_log('', sprintf('Skipping %s: session, cluster id, or lock group is missing.', neuron_file));
        clear neuron
        continue
    end
    session_name = session_acg_char(neuron.session);
    cluster_id = double(neuron.cluster_id);
    if ~isscalar(cluster_id) || ~isfinite(cluster_id)
        session_acg_log('', sprintf('Skipping %s: cluster id is not one finite number.', neuron_file));
        clear neuron
        continue
    end
    key = sprintf('%s|%g', session_name, cluster_id);
    if any(strcmp(keys, key))
        session_acg_log('', sprintf('Skipping duplicate %s', key));
        clear neuron
        continue
    end
    rec = blank_record(curve_names);
    rec.session = session_name;
    rec.cluster_id = cluster_id;
    rec.lock_group = session_acg_char(neuron.lock_group);
    for c = 1:numel(curve_names)
        rec.curves.(curve_names{c}) = curve_percentile(neuron, curve_names{c});
    end
    keys{end + 1, 1} = key; %#ok<AGROW>
    if isempty(records)
        records = rec;
    else
        records(end + 1, 1) = rec; %#ok<AGROW>
    end
    clear neuron rec
    if mod(f, 25) == 0
        session_acg_log('', sprintf('Read %d/%d neuron files. %d neurons kept.', ...
            f, n_files, numel(records)));
    end
end
session_acg_log('', sprintf('Kept %d neurons from %d files.', numel(records), n_files));
end


function records = empty_records()
records = struct('session', {}, 'cluster_id', {}, 'lock_group', {}, 'curves', {});
end


function rec = blank_record(curve_names)
rec.session = '';
rec.cluster_id = NaN;
rec.lock_group = '';
rec.curves = struct();
for i = 1:numel(curve_names)
    rec.curves.(curve_names{i}) = empty_curve();
end
end


function curve = curve_percentile(neuron, curve_name)
curve = empty_curve();
if ~isfield(neuron, 'scores') || ~isfield(neuron.scores, curve_name)
    return
end
score = neuron.scores.(curve_name);
if isfield(score, 'oscillates') && isscalar(score.oscillates)
    curve.oscillates = logical(score.oscillates);
end
if isfield(score, 'first_lag') && isscalar(score.first_lag)
    curve.first_lag = double(score.first_lag);
end
if isfield(score, 'first_p') && isscalar(score.first_p)
    curve.first_p = double(score.first_p);
end
if ~isfield(score, 'pctl') || ~isfield(score, 'pctl_lag')
    return
end
y = double(score.pctl(:))';
lag = double(score.pctl_lag(:))';
if isempty(y) || numel(y) ~= numel(lag) || ~any(isfinite(y))
    return
end
curve.pctl = y;
curve.pctl_lag = lag;
end


function curve = empty_curve()
curve.pctl = [];
curve.pctl_lag = [];
curve.oscillates = false;
curve.first_lag = NaN;
curve.first_p = NaN;
end
