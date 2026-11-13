function grand = session_acg_grand_average(save_root)
%SESSION_ACG_GRAND_AVERAGE  Mean ACG of each lock group from saved nulls.
%
% Every neuron file keeps its own 500 null autocorrelograms. Null row k of
% the grand mean is the mean of null row k across neurons, after each row
% has been scaled to sum to 1. Trough-locked and peak-locked neurons are
% averaged again together as delta-locked. The delta-peak test is then
% applied to that grand mean. One session folder is read at a time.

param_file = fullfile(save_root, 'parameters.mat');
if exist(param_file, 'file') ~= 2
    error('session_acg:grand', ...
        'No parameters.mat in %s. Run the estimate before the grand average.', save_root);
end
stored = load(param_file, 'settings');
settings = session_acg_settings(stored.settings);
files = dir(fullfile(save_root, '*', 'neuron_*.mat'));
files = files(~[files.isdir]);
if isempty(files)
    error('session_acg:grand', ...
        'No neuron files under %s. The estimate has not saved a session yet.', save_root);
end
assert_session_settings(files, settings);

edges = settings.bin_edges(:)';
bin_w = edges(2) - edges(1);
centers = edges(1:end - 1) + bin_w / 2;
curves = session_acg_names('curves');
groups = session_acg_names('groups');
grand = empty_grand(groups, curves, settings.n_shuffles, numel(centers));
grand.settings = settings;
grand.centers = centers;
grand.keys = cell(0, 1);

for f = 1:numel(files)
    neuron_file = fullfile(files(f).folder, files(f).name);
    loaded = load(neuron_file, 'neuron');
    neuron = loaded.neuron;
    key = sprintf('%s|%g', char(string(neuron.session)), neuron.cluster_id);
    if any(strcmp(grand.keys, key))
        session_acg_log('', sprintf('Skipping duplicate %s', key));
        clear loaded neuron
        continue
    end
    grand = add_neuron(grand, neuron, curves);
    grand.keys{end + 1, 1} = key; %#ok<AGROW>
    clear loaded neuron
end

grand = score_grand(grand, groups, curves, centers, settings.smooth_n, settings.oscillation_alpha);
grand.report = oscillation_report(grand, groups, curves);
report_file = fullfile(save_root, 'grand_average_oscillation.txt');
fid = fopen(report_file, 'w');
if fid >= 0
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s', grand.report);
    clear cleanup
end
fprintf('%s', grand.report);
log_file = fullfile(save_root, 'estimate_log.txt');
if exist(log_file, 'file') == 2
    log_fid = fopen(log_file, 'a');
    if log_fid >= 0
        log_cleanup = onCleanup(@() fclose(log_fid)); %#ok<NASGU>
        fprintf(log_fid, '%s', grand.report);
        clear log_cleanup
    end
end
save(fullfile(save_root, 'grand_average.mat'), 'grand', 'settings', '-v7');
end


function assert_session_settings(files, settings)
folders = unique({files.folder});
for i = 1:numel(folders)
    param_file = fullfile(folders{i}, 'parameters.mat');
    if exist(param_file, 'file') ~= 2
        error('session_acg:parameters', ...
            '%s has neuron files but no parameters.mat.', folders{i});
    end
    saved = load(param_file, 'settings');
    if ~session_acg_parameters_match(saved.settings, settings)
        error('session_acg:parameters', ...
            '%s does not match the settings in the estimate folder.', folders{i});
    end
end
end


function grand = empty_grand(groups, curves, n_shuffles, n_bins)
grand.n_shuffles = n_shuffles;
grand.n_bins = n_bins;
grand.group_order = groups;
grand.curve_order = curves;
grand.groups = struct();
grand.counts = struct();
for g = 1:numel(groups)
    for c = 1:numel(curves)
        grand.groups.(groups{g}).(curves{c}) = empty_block(n_shuffles, n_bins);
        grand.counts.(groups{g}).(curves{c}) = struct('n', 0, 'n_oscillates', 0);
    end
end
end


function block = empty_block(n_shuffles, n_bins)
block.sum_obs = zeros(1, n_bins);
block.n_obs = zeros(1, n_bins);
block.sum_null = zeros(n_shuffles, n_bins);
block.n_null = zeros(n_shuffles, n_bins);
block.n_neurons = 0;
block.mean = nan(1, n_bins);
block.null_mean = nan(n_shuffles, n_bins);
block.result = score_acg_oscillation(nan(1, n_bins), [], [], 1, 1);
end


function grand = add_neuron(grand, neuron, curves)
targets = {char(string(neuron.lock_group))};
if ~any(strcmp(targets{1}, {'trough', 'peak', 'nonentrained'}))
    error('session_acg:group', ...
        'Neuron %g in %s has lock group %s.', neuron.cluster_id, neuron.session, targets{1});
end
if any(strcmp(targets{1}, {'trough', 'peak'}))
    targets{end + 1} = 'delta'; %#ok<AGROW>
end
for c = 1:numel(curves)
    curve = curves{c};
    if ~isfield(neuron.curves, curve)
        continue
    end
    observed = neuron.curves.(curve).observed;
    if ~any(isfinite(observed))
        continue
    end
    null_rows = neuron.curves.(curve).null;
    oscillates = false;
    if isfield(neuron, 'scores') && isfield(neuron.scores, curve) ...
            && isfield(neuron.scores.(curve), 'oscillates')
        oscillates = logical(neuron.scores.(curve).oscillates);
    end
    for g = 1:numel(targets)
        grand.groups.(targets{g}).(curve) = add_curve( ...
            grand.groups.(targets{g}).(curve), observed, null_rows, grand.n_shuffles);
        grand.counts.(targets{g}).(curve).n = grand.counts.(targets{g}).(curve).n + 1;
        if oscillates
            grand.counts.(targets{g}).(curve).n_oscillates = ...
                grand.counts.(targets{g}).(curve).n_oscillates + 1;
        end
    end
end
end


function block = add_curve(block, observed, null_rows, n_shuffles)
observed = double(observed(:))';
finite = isfinite(observed);
block.sum_obs(finite) = block.sum_obs(finite) + observed(finite);
block.n_obs(finite) = block.n_obs(finite) + 1;
block.n_neurons = block.n_neurons + 1;
if isempty(null_rows)
    return
end
null_rows = double(null_rows);
if size(null_rows, 1) ~= n_shuffles || size(null_rows, 2) ~= numel(observed)
    error('session_acg:null', ...
        'Saved nulls must be %d shuffles by %d bins.', n_shuffles, numel(observed));
end
finite_null = isfinite(null_rows);
block.sum_null(finite_null) = block.sum_null(finite_null) + null_rows(finite_null);
block.n_null(finite_null) = block.n_null(finite_null) + 1;
end


function grand = score_grand(grand, groups, curves, centers, smooth_n, alpha)
for g = 1:numel(groups)
    for c = 1:numel(curves)
        block = grand.groups.(groups{g}).(curves{c});
        if block.n_neurons < 1
            grand.groups.(groups{g}).(curves{c}) = block;
            continue
        end
        mean_row = nan(1, grand.n_bins);
        observed_ok = block.n_obs > 0;
        mean_row(observed_ok) = block.sum_obs(observed_ok) ./ block.n_obs(observed_ok);
        null_mean = nan(grand.n_shuffles, grand.n_bins);
        null_ok = block.n_null > 0;
        null_mean(null_ok) = block.sum_null(null_ok) ./ block.n_null(null_ok);
        block.mean = mean_row;
        block.null_mean = null_mean;
        block.result = score_acg_oscillation(mean_row, null_mean, centers, smooth_n, alpha);
        grand.groups.(groups{g}).(curves{c}) = block;
    end
end
end


function report = oscillation_report(grand, groups, curves)
group_title = struct( ...
    'trough', 'Trough-locked', ...
    'peak', 'Peak-locked', ...
    'delta', 'Delta-locked', ...
    'nonentrained', 'Non-entrained');
curve_title = struct( ...
    'session', 'Whole-session spikes', ...
    'fooof_low', 'Low FOOOF-delta windows', ...
    'fooof_high', 'High FOOOF-delta windows', ...
    'hilbert_low', 'Low Hilbert-delta windows', ...
    'hilbert_high', 'High Hilbert-delta windows');
stamp = char(datetime('now', 'TimeZone', 'Europe/Berlin', 'Format', 'yyyy-MM-dd HH:mm:ss'));
settings = grand.settings;
lines = {
    sprintf('Session ACG grand average, %s Europe/Berlin', stamp)
    sprintf('Shuffles %g, windows %g s, minimum spikes %g, tail %g%%, oscillation alpha %g, method %s.', ...
        settings.n_shuffles, settings.window_s, settings.min_spikes, ...
        settings.tail_percent, settings.oscillation_alpha, settings.acg_method)
    sprintf('Areas: %s.', strjoin(settings.areas_to_keep, ', '))
    'A neuron oscillates in delta when its autocorrelogram has a delta peak and that peak prominence is larger than the null.'
    'The grand-mean test scores the mean autocorrelogram against the mean of each saved shuffle row.'
    ''
    };
for c = 1:numel(curves)
    for g = 1:numel(groups)
        counts = grand.counts.(groups{g}).(curves{c});
        result = grand.groups.(groups{g}).(curves{c}).result;
        lines{end + 1} = sprintf('%s, %s', curve_title.(curves{c}), group_title.(groups{g})); %#ok<AGROW>
        lines{end + 1} = sprintf('  neurons in the mean: %d', counts.n); %#ok<AGROW>
        if counts.n > 0
            lines{end + 1} = sprintf('  neurons with a delta oscillation: %d/%d', ...
                counts.n_oscillates, counts.n); %#ok<AGROW>
        end
        lines{end + 1} = grand_peak_line(result); %#ok<AGROW>
        lines{end + 1} = ''; %#ok<AGROW>
    end
end
report = strjoin(lines, newline);
report = [report, newline];
end


function line = grand_peak_line(result)
if ~isfinite(result.first_p)
    line = '  grand mean: no delta-peak test (no neurons, or the null was not saved).';
    return
end
if result.oscillates
    call = 'oscillates';
else
    call = 'does not clear the delta-peak test';
end
lag_text = 'no delta peak';
if result.first_present && isfinite(result.first_lag)
    lag_text = sprintf('peak %.3f s', result.first_lag);
end
percentile = percentile_at_peak(result);
if isfinite(percentile)
    percentile_text = sprintf(', percentile at that lag %.1f', percentile);
else
    percentile_text = '';
end
line = sprintf('  grand mean: %s, prominence %.3g, shuffle p %.3g, %s%s.', ...
    lag_text, result.first_prominence, result.first_p, call, percentile_text);
end


function percentile = percentile_at_peak(result)
percentile = NaN;
if ~result.first_present || isempty(result.pctl) || ~isfinite(result.first_lag)
    return
end
[~, index] = min(abs(result.pctl_lag - result.first_lag));
percentile = result.pctl(index);
end
