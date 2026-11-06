%% Analyze_locking_period_acg_shuffle
% How exceptional the delta-range autocorrelogram peaks are.
%
% This is separate from Analyze_locking_period_acg because each neuron is
% repeated many times. Spikes inside a firing run are replaced by the same
% number of independent uniform times in that run. Intervals are not kept,
% so the autocorrelogram of a rhythmic train is destroyed. The low 25% and
% high 25% runs (PPC, and delta power) are then pooled and scored with the
% same first-peak and second-cycle measures as the main script.
%
% A peak is exceptional when its prominence is larger than the prominences
% from those shuffles. The one-sided p-value is
% (1 + number of shuffles at least as large) / (number of shuffles + 1).
% The second cycle is compared only when the observed curve could measure
% it. A shuffle that finds no second cycle counts as prominence 0.
%
% The result file is rewritten after every neuron. Stop the script and run
% it again: finished neurons are skipped. Lower n_shuffles to shorten a
% first pass. Only runs that fall in a plotted tail are shuffled.
% Only areas_to_keep are analyzed. SupCol is the superior-colliculus
% label in the depth tables and is kept as SC. VLPAG and DLPAG are not
% treated as LPAG.
%
% Run this file, not the unit test, to produce the figures:
%   Analyze_locking_period_acg_shuffle

%% Paths and parameters
data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
run('\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
phase_folder   = fullfile(data_root, 'Analysis results', 'phase locking data');
npx_root       = fullfile(data_root, 'NPX data', 'NPX raw data');

n_shuffles       = 200;
resume           = true;
alpha            = 0.01;
non_entrained_p  = 0.1;
low_percent      = 25;
high_percent     = 75;
acg_lag          = 1;
smooth_half      = 0.04;
acg_refractory_s = 0.005;
areas_to_keep = {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'};

group_order = {'trough', 'peak', 'nonentrained'};
group_title = {'Trough-locked', 'Peak-locked', 'Non-entrained'};
checkpoint = fullfile(parsing_folder, ...
    sprintf('locking_period_acg_shuffle_n%d.mat', n_shuffles));

parameters = struct( ...
    'n_shuffles', n_shuffles, ...
    'alpha', alpha, ...
    'non_entrained_p', non_entrained_p, ...
    'low_percent', low_percent, ...
    'high_percent', high_percent, ...
    'acg_lag', acg_lag, ...
    'smooth_half', smooth_half, ...
    'acg_refractory_s', acg_refractory_s, ...
    'areas_to_keep', {areas_to_keep}, ...
    'null', 'uniform spike times inside each firing run');

records = {};
if resume && exist(checkpoint, 'file') == 2
    stored = load(checkpoint, 'shuffle_records', 'parameters');
    if ~shuffle_checkpoint_matches(stored, parameters)
        error('Analyze_locking_period_acg_shuffle:checkpoint', ...
            ['%s was built with different parameters. Delete it or set ' ...
            'resume = false.'], checkpoint);
    end
    records = stored.shuffle_records;
    if ~iscell(records)
        records = num2cell(records(:));
    end
    records = records(:);
    records = records_in_areas(records, areas_to_keep);
    disp(['Resuming from ', num2str(numel(records)), ' neurons in ', checkpoint])
else
    disp(['Checkpoint: ', checkpoint])
end
done_keys = neuron_keys(records);

%% Delta-lock labels
load(fullfile(phase_folder, 'delta_all_neurons_v2.mat'), 'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) = {'isRt'};
id_all = double(cluster_id_column(all_neurons));
n_cells = height(all_neurons);
lock_group = cell(n_cells, 1);
for i = 1:n_cells
    lock_group{i} = delta_lock_group( ...
        all_neurons.EntireSession.PPCPval(i), ...
        all_neurons.EntireSession.PreferedAngle(i), ...
        alpha, non_entrained_p);
end

%% One parsing file at a time
parsing_files = list_parsing_mat_files(parsing_folder);
if isempty(parsing_files)
    error('Analyze_locking_period_acg_shuffle:empty', ...
        'No spike-train parsing files in %s.', parsing_folder);
end

loaded_codes = {};
edges = [];
centers = [];
n_bins = [];
bin_w = [];
smooth_n = [];

for f = 1:numel(parsing_files)
    disp(['Loading ', parsing_files(f).name, ' ...'])
    S = load(parsing_files(f).path, 'parsing_struct');
    if ~isfield(S, 'parsing_struct')
        warning('Analyze_locking_period_acg_shuffle:contents', ...
            'No parsing_struct in %s', parsing_files(f).name);
        clear S
        continue
    end
    while ~isempty(S.parsing_struct)
        ps = S.parsing_struct(1);
        S.parsing_struct(1) = [];
        if isfield(ps, 'ALL_TABLES')
            ps = rmfield(ps, 'ALL_TABLES');
        end
        animal_code = char(string(ps.animal_code));
        if ismember(animal_code, loaded_codes)
            disp(['  already loaded, skip ', animal_code])
            clear ps
            continue
        end
        saved_edges = ps.config.autocorr_bin_edges(:)';
        if isempty(edges)
            if isfield(ps.config, 'autocorr_bin_size') && isfinite(ps.config.autocorr_bin_size)
                bin_w = ps.config.autocorr_bin_size;
            else
                bin_w = mean(diff(saved_edges));
            end
            edges = (-acg_lag):bin_w:acg_lag;
            centers = edges(1:end - 1) + bin_w / 2;
            n_bins = numel(centers);
            smooth_n = max(1, round(smooth_half / bin_w));
        end
        saved_covers = numel(saved_edges) == numel(edges) ...
            && max(abs(saved_edges - edges)) < 1e-6;

        ci = ps.cluster_info;
        id_col = double(cluster_id_column(ci));
        run_cell = ps.run_events_table_per_cell;
        ppc_col = ppc_column_index(ps);
        already_corrected = isfield(ps, 'config') && isfield(ps.config, 'acg_edge_corrected') ...
            && ps.config.acg_edge_corrected;
        if isempty(ppc_col)
            warning('Analyze_locking_period_acg_shuffle:ppc', ...
                'No PPC column in %s; session skipped.', animal_code);
            clear ps
            continue
        end

        [session_times, session_clusters] = load_pag_spikes(npx_root, animal_code);
        if isempty(session_times)
            warning('Analyze_locking_period_acg_shuffle:spikes', ...
                'Spike times for %s were not found. Session skipped.', animal_code);
            clear ps
            continue
        end

        for nn = 1:height(ci)
            this_id = id_col(nn);
            this_key = neuron_key(animal_code, this_id);
            if ismember(this_key, done_keys)
                continue
            end
            hit = find(session_matches(all_neurons.session, animal_code) & id_all == this_id, 1);
            if isempty(hit)
                continue
            end
            this_area = cluster_area_char(ci, nn);
            if isempty(strtrim(this_area))
                this_area = cluster_area_char(all_neurons, hit);
            end
            this_area = canonical_kept_area(this_area, areas_to_keep);
            if isempty(this_area)
                continue
            end
            this_group = lock_group{hit};
            if ~any(strcmp(group_order, this_group))
                continue
            end
            rt = run_cell{nn};
            if isempty(rt) || ~istable(rt) || ~ismember('autoccorrelograms', rt.Properties.VariableNames)
                continue
            end
            if size(rt.phase_entrainment, 2) < ppc_col
                continue
            end
            if ~ismember('RunStartTime', rt.Properties.VariableNames) ...
                    || ~ismember('RunEndTime', rt.Properties.VariableNames)
                continue
            end
            if saved_covers
                acg = rt.autoccorrelograms;
                this_corrected = already_corrected;
            else
                acg = run_window_acgs(session_times(session_clusters == this_id), ...
                    rt.RunStartTime, rt.RunEndTime, edges);
                this_corrected = true;
            end
            if size(acg, 2) ~= n_bins || size(acg, 1) ~= height(rt)
                warning('Analyze_locking_period_acg_shuffle:width', ...
                    'ACG size mismatch in %s id %g; cell skipped.', animal_code, this_id);
                continue
            end

            run_length = acg_run_lengths(rt);
            ppc = rt.phase_entrainment(:, ppc_col);
            power = run_delta_power(rt);
            pooled = pool_locking_period_acgs(acg, edges, run_length, ppc, ...
                this_corrected, low_percent, high_percent);
            pooled_power = pool_locking_period_acgs(acg, edges, run_length, power, ...
                this_corrected, low_percent, high_percent);
            lock_ok = split_is_usable(pooled);
            power_ok = split_is_usable(pooled_power);
            if ~lock_ok && ~power_ok
                continue
            end

            use = false(height(rt), 1);
            if lock_ok
                use = use | pooled.low_index(:) | pooled.high_index(:);
            end
            if power_ok
                use = use | pooled_power.low_index(:) | pooled_power.high_index(:);
            end
            spikes = session_times(session_clusters == this_id);
            disp(sprintf('  %s id %g (%s): %d runs, %d shuffles', ...
                animal_code, this_id, this_group, sum(use), n_shuffles))
            rng(neuron_shuffle_seed(animal_code, this_id), 'twister');
            nulls = empty_nulls(n_shuffles);
            t0 = tic;
            for s = 1:n_shuffles
                shuffled = shuffle_spikes_within_runs(spikes, ...
                    rt.RunStartTime(use), rt.RunEndTime(use));
                acg_s = nan(height(rt), n_bins);
                acg_s(use, :) = run_window_acgs(shuffled, ...
                    rt.RunStartTime(use), rt.RunEndTime(use), edges);
                if lock_ok
                    sh_lock = pool_locking_period_acgs(acg_s, edges, run_length, ppc, ...
                        true, low_percent, high_percent);
                    nulls = store_null_tail(nulls, 'locking', 'low', s, ...
                        score_pooled_row(sh_lock.low, centers, smooth_n, acg_refractory_s));
                    nulls = store_null_tail(nulls, 'locking', 'high', s, ...
                        score_pooled_row(sh_lock.high, centers, smooth_n, acg_refractory_s));
                end
                if power_ok
                    sh_power = pool_locking_period_acgs(acg_s, edges, run_length, power, ...
                        true, low_percent, high_percent);
                    nulls = store_null_tail(nulls, 'power', 'low', s, ...
                        score_pooled_row(sh_power.low, centers, smooth_n, acg_refractory_s));
                    nulls = store_null_tail(nulls, 'power', 'high', s, ...
                        score_pooled_row(sh_power.high, centers, smooth_n, acg_refractory_s));
                end
                if mod(s, 25) == 0
                    disp(sprintf('    shuffle %d/%d (%.0f s)', s, n_shuffles, toc(t0)))
                end
            end

            rec = blank_record(n_shuffles);
            rec.session = animal_code;
            rec.cluster_id = this_id;
            rec.area = this_area;
            rec.lock_group = this_group;
            if lock_ok
                rec.locking = pack_split(pooled, centers, smooth_n, acg_refractory_s, ...
                    nulls.locking);
            end
            if power_ok
                rec.power = pack_split(pooled_power, centers, smooth_n, acg_refractory_s, ...
                    nulls.power);
            end
            records{end + 1, 1} = rec; %#ok<SAGROW>
            done_keys{end + 1, 1} = this_key; %#ok<SAGROW>
            shuffle_records = records; %#ok<NASGU>
            save(checkpoint, 'shuffle_records', 'parameters');
            disp(['    saved ', num2str(numel(records)), ' neurons'])
        end

        loaded_codes{end + 1} = animal_code; %#ok<AGROW>
        clear ps ci run_cell session_times session_clusters
    end
    clear S
end

disp(['Neurons with a shuffle null: ', num2str(numel(records))])
plot_shuffle_exception(records, group_order, group_title, 'locking', ...
    'Mean PPC of the runs');
plot_shuffle_exception(records, group_order, group_title, 'power', ...
    'Mean delta power of the runs');


function rec = blank_record(n_shuffles)
rec.session = '';
rec.cluster_id = NaN;
rec.area = '';
rec.lock_group = '';
rec.locking = blank_split(n_shuffles);
rec.power = blank_split(n_shuffles);
end


function split = blank_split(n_shuffles)
split.usable = false;
split.score_low = NaN;
split.score_high = NaN;
split.low = blank_tail(n_shuffles);
split.high = blank_tail(n_shuffles);
end


function tail = blank_tail(n_shuffles)
tail.first_prominence = NaN;
tail.first_lag = NaN;
tail.first_present = false;
tail.first_p = NaN;
tail.second_prominence = NaN;
tail.second_lag = NaN;
tail.second_present = false;
tail.second_measurable = false;
tail.second_p = NaN;
tail.null_first = nan(n_shuffles, 1);
tail.null_second = nan(n_shuffles, 1);
end


function nulls = empty_nulls(n_shuffles)
nulls.locking.low.first = nan(n_shuffles, 1);
nulls.locking.low.second = nan(n_shuffles, 1);
nulls.locking.high.first = nan(n_shuffles, 1);
nulls.locking.high.second = nan(n_shuffles, 1);
nulls.power.low.first = nan(n_shuffles, 1);
nulls.power.low.second = nan(n_shuffles, 1);
nulls.power.high.first = nan(n_shuffles, 1);
nulls.power.high.second = nan(n_shuffles, 1);
end


function nulls = store_null_tail(nulls, split_name, tail_name, sample, peaks)
nulls.(split_name).(tail_name).first(sample) = null_first_value(peaks);
nulls.(split_name).(tail_name).second(sample) = null_second_value(peaks);
end


function v = null_first_value(peaks)
if ~peaks.first_measurable
    v = NaN;
else
    v = peaks.first_prominence;
end
end


function v = null_second_value(peaks)
if ~peaks.first_measurable
    v = NaN;
elseif peaks.second_measurable
    v = peaks.second_prominence;
else
    v = 0;
end
end


function split = pack_split(pooled, centers, smooth_n, refractory, null_split)
split = blank_split(numel(null_split.low.first));
split.usable = true;
split.score_low = pooled.score_low;
split.score_high = pooled.score_high;
split.low = pack_tail(pooled.low, centers, smooth_n, refractory, ...
    null_split.low.first, null_split.low.second);
split.high = pack_tail(pooled.high, centers, smooth_n, refractory, ...
    null_split.high.first, null_split.high.second);
end


function tail = pack_tail(acg_row, centers, smooth_n, refractory, null_first, null_second)
peaks = score_pooled_row(acg_row, centers, smooth_n, refractory);
tail = blank_tail(numel(null_first));
tail.first_prominence = peaks.first_prominence;
tail.first_lag = peaks.first_lag;
tail.first_present = peaks.first_present;
tail.second_prominence = peaks.second_prominence;
tail.second_lag = peaks.second_lag;
tail.second_present = peaks.second_present;
tail.second_measurable = peaks.second_measurable;
tail.null_first = null_first(:);
tail.null_second = null_second(:);
tail.first_p = exception_p(peaks.first_prominence, null_first);
if peaks.second_measurable
    tail.second_p = exception_p(peaks.second_prominence, null_second);
end
end


function peaks = score_pooled_row(acg_row, centers, smooth_n, refractory)
if isempty(acg_row) || ~any(isfinite(acg_row))
    peaks = acg_delta_peaks(nan(1, numel(centers)), centers, smooth_n);
    return
end
row = blank_acg_zero_lag(acg_row(:)', centers, refractory);
row = normalize_acg_rows(row);
peaks = acg_delta_peaks(row, centers, smooth_n);
end


function p = exception_p(observed, nulls)
nulls = nulls(isfinite(nulls));
if ~isfinite(observed) || isempty(nulls)
    p = NaN;
    return
end
p = (1 + sum(nulls >= observed)) / (numel(nulls) + 1);
end


function ok = split_is_usable(pooled)
ok = pooled.n_low >= 1 && pooled.n_high >= 1 ...
    && any(isfinite(pooled.low)) && any(isfinite(pooled.high));
end


function match = shuffle_checkpoint_matches(stored, parameters)
match = isfield(stored, 'parameters') && isfield(stored, 'shuffle_records') ...
    && stored.parameters.n_shuffles == parameters.n_shuffles ...
    && stored.parameters.alpha == parameters.alpha ...
    && stored.parameters.non_entrained_p == parameters.non_entrained_p ...
    && stored.parameters.low_percent == parameters.low_percent ...
    && stored.parameters.high_percent == parameters.high_percent ...
    && stored.parameters.acg_lag == parameters.acg_lag;
end


function keys = neuron_keys(records)
keys = cell(numel(records), 1);
for i = 1:numel(records)
    keys{i} = neuron_key(records{i}.session, records{i}.cluster_id);
end
end


function key = neuron_key(session, cluster_id)
key = [char(string(session)), '|', sprintf('%.0f', double(cluster_id))];
end


function seed = neuron_shuffle_seed(session, cluster_id)
letters = double(char(string(session)));
seed = mod(sum(letters .* (1:numel(letters))) + round(double(cluster_id)) * 997, 2^31 - 2);
seed = seed + 1;
end


function plot_shuffle_exception(records, group_order, group_title, split_name, score_label)
if isempty(records)
    return
end
if strcmp(split_name, 'locking')
    split_label = 'locking';
else
    split_label = 'delta power';
end
measures = {'first_p', 'second_p'};
measure_labels = {'First-peak exceptionality', 'Second-cycle exceptionality'};
colors = [0.20 0.45 0.80; 0.85 0.25 0.20; 0.45 0.45 0.45];
figure('Name', ['ACG peak shuffle, ', split_label], 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.08 0.08 0.84 0.8])
for m = 1:numel(measures)
    [p_low, p_high, ok, x_low, x_high, groups] = paired_shuffle_values( ...
        records, split_name, measures{m});
    y_low = -log10(p_low);
    y_high = -log10(p_high);
    subplot(2, 2, 2 * m - 1)
    plot_paired_exception(groups(ok), y_low(ok), y_high(ok), ...
        group_order, colors, measure_labels{m})
    subplot(2, 2, 2 * m)
    plot_exception_against_score(groups(ok), x_low(ok), x_high(ok), ...
        y_low(ok), y_high(ok), group_order, colors, score_label, measure_labels{m})
    report_exception_fraction(groups, p_low, p_high, ok, group_order, group_title, ...
        split_label, measure_labels{m});
end
hold on
handles = gobjects(1, numel(group_order));
for g = 1:numel(group_order)
    handles(g) = plot(nan, nan, '-', 'Color', colors(g, :), 'LineWidth', 2.5);
end
legend(handles, group_title, 'Location', 'eastoutside')
end


function [p_low, p_high, ok, x_low, x_high, groups] = paired_shuffle_values(records, split_name, measure_name)
n = numel(records);
p_low = nan(n, 1);
p_high = nan(n, 1);
x_low = nan(n, 1);
x_high = nan(n, 1);
groups = cell(n, 1);
for i = 1:n
    groups{i} = records{i}.lock_group;
    block = records{i}.(split_name);
    x_low(i) = block.score_low;
    x_high(i) = block.score_high;
    p_low(i) = block.low.(measure_name);
    p_high(i) = block.high.(measure_name);
    if strcmp(measure_name, 'second_p')
        if ~block.low.second_measurable || ~block.high.second_measurable
            p_low(i) = NaN;
            p_high(i) = NaN;
        end
    end
end
ok = isfinite(p_low) & isfinite(p_high) & p_low > 0 & p_high > 0;
end


function plot_paired_exception(lock_group, y_low, y_high, group_order, colors, y_label)
hold on
threshold = -log10(0.05);
plot([0.6 2.4], [threshold threshold], '--', 'Color', [0.45 0.45 0.45], ...
    'HandleVisibility', 'off')
for g = 1:numel(group_order)
    idx = strcmp(lock_group, group_order{g});
    if ~any(idx)
        continue
    end
    plot([ones(sum(idx), 1), 2 * ones(sum(idx), 1)]', [y_low(idx), y_high(idx)]', '-', ...
        'Color', [colors(g, :), 0.25], 'HandleVisibility', 'off')
    plot([1 2], [mean(y_low(idx), 'omitnan'), mean(y_high(idx), 'omitnan')], '-', ...
        'Color', colors(g, :), 'LineWidth', 2.5, 'HandleVisibility', 'off')
end
xlim([0.6 2.4])
set(gca, 'XTick', [1 2], 'XTickLabel', {'Low 25%', 'High 25%'})
ylabel([y_label, ', -log10 p'])
title(sprintf('Paired, n = %d. Dashed line is p = 0.05', numel(y_low)))
box off
end


function plot_exception_against_score(lock_group, x_low, x_high, y_low, y_high, group_order, colors, x_label, y_label)
hold on
span = [x_low(:); x_high(:)];
span = span(isfinite(span));
if numel(span) >= 2
    plot([min(span) max(span)], [-log10(0.05) -log10(0.05)], '--', ...
        'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
end
for g = 1:numel(group_order)
    idx = strcmp(lock_group, group_order{g});
    if ~any(idx)
        continue
    end
    plot([x_low(idx), x_high(idx)]', [y_low(idx), y_high(idx)]', '-', ...
        'Color', [colors(g, :), 0.3], 'HandleVisibility', 'off')
end
xlabel(x_label)
ylabel([y_label, ', -log10 p'])
title('Each line is one neuron')
box off
end


function report_exception_fraction(groups, p_low, p_high, ok, group_order, group_title, split_label, measure_label)
for g = 1:numel(group_order)
    idx = ok & strcmp(groups, group_order{g});
    n = sum(idx);
    if n < 1
        continue
    end
    disp(sprintf('  %s, %s, %s: low p<0.05 = %d/%d, high p<0.05 = %d/%d', ...
        group_title{g}, split_label, measure_label, ...
        sum(p_low(idx) < 0.05), n, sum(p_high(idx) < 0.05), n))
end
end


function power = run_delta_power(rt)
if ~ismember('lfp_param', rt.Properties.VariableNames) || size(rt.lfp_param, 2) < 2
    power = nan(height(rt), 1);
    return
end
power = double(rt.lfp_param(:, 2)) .^ 2;
end


function T = acg_run_lengths(rt)
if ismember('RunLength', rt.Properties.VariableNames)
    T = double(rt.RunLength);
elseif ismember('RunStartTime', rt.Properties.VariableNames) ...
        && ismember('RunEndTime', rt.Properties.VariableNames)
    T = double(rt.RunEndTime) - double(rt.RunStartTime);
else
    T = nan(height(rt), 1);
end
T = T(:);
end


function files = list_parsing_mat_files(parsing_folder)
files = struct('path', {}, 'name', {}, 'kind', {});
combined = fullfile(parsing_folder, 'spike_train_parsing_structure.mat');
if exist(combined, 'file')
    files(end + 1) = struct('path', combined, ...
        'name', 'spike_train_parsing_structure.mat', 'kind', 'pooled'); %#ok<AGROW>
end
flist = dir(fullfile(parsing_folder, '*spike_train_parsing*.mat'));
skip = {'spike_train_parsing_structure.mat', 'spike_train_parsing_animal_names.mat', ...
    'spike_train_parsing_animal_names_v1.mat'};
for f = 1:numel(flist)
    if ismember(flist(f).name, skip)
        continue
    end
    files(end + 1) = struct( ...
        'path', fullfile(flist(f).folder, flist(f).name), ...
        'name', flist(f).name, 'kind', 'session'); %#ok<AGROW>
end
end


function id_col = cluster_id_column(ci)
if ismember('cluster_id', ci.Properties.VariableNames)
    id_col = ci.cluster_id;
elseif ismember('id', ci.Properties.VariableNames)
    id_col = ci.id;
else
    error('Analyze_locking_period_acg_shuffle:id', ...
        'cluster table has neither cluster_id nor id');
end
end


function col = ppc_column_index(ps)
col = [];
if isfield(ps, 'config') && isfield(ps.config, 'phase_entrainment_fields')
    hit = find(strcmp(ps.config.phase_entrainment_fields, 'PPC'), 1);
    if ~isempty(hit)
        col = hit;
        return
    end
end
for j = 1:numel(ps.run_events_table_per_cell)
    rt = ps.run_events_table_per_cell{j};
    if ~isempty(rt) && istable(rt) && ismember('phase_entrainment', rt.Properties.VariableNames)
        if size(rt.phase_entrainment, 2) >= 5
            col = 5;
        end
        return
    end
end
end


function tf = session_matches(session_col, animal_code)
if iscell(session_col)
    tf = strcmp(session_col, animal_code);
    empty = cellfun(@isempty, session_col);
    tf(empty) = false;
    nested = cellfun(@iscell, session_col);
    if any(nested)
        tf(nested) = cellfun(@(x) strcmp(x{1}, animal_code), session_col(nested));
    end
else
    tf = strcmp(string(session_col), animal_code);
end
end


function a = cluster_area_char(ci, nn)
if ~ismember('area', ci.Properties.VariableNames)
    a = '';
    return
end
a = ci.area(nn, :);
if iscell(a)
    a = a{1};
    if iscell(a)
        a = a{1};
    end
end
a = char(string(a));
end


function name = canonical_kept_area(label, areas_to_keep)
% SupCol and SC are the superior colliculus. VLPAG does not match LPAG.
name = '';
text = strtrim(char(string(label)));
if isempty(text)
    return
end
tokens = regexp(text, '[,;/]+', 'split');
for k = 1:numel(tokens)
    token = upper(regexprep(strtrim(tokens{k}), '\s+', ''));
    if any(strcmp(token, {'SC', 'SUPCOL', 'SUPERIORCOLLICULUS'}))
        token = 'SC';
    end
    hit = strcmpi(areas_to_keep, token);
    if any(hit)
        name = areas_to_keep{find(hit, 1)};
        return
    end
end
end


function records = records_in_areas(records, areas_to_keep)
keep = false(numel(records), 1);
for i = 1:numel(records)
    if ~isstruct(records{i}) || ~isfield(records{i}, 'area')
        continue
    end
    keep(i) = ~isempty(canonical_kept_area(records{i}.area, areas_to_keep));
end
n_drop = sum(~keep);
records = records(keep);
if n_drop > 0
    disp(['Dropped ', num2str(n_drop), ' checkpoint neurons outside the selected areas'])
end
end


function Y = blank_acg_zero_lag(X, centers, flank_s)
Y = X;
if isempty(X) || isempty(centers)
    return
end
t = centers(:)';
dt = median(diff(t));
if ~isfinite(dt) || dt <= 0
    dt = flank_s;
end
half = abs(dt) / 2;
lo = t - half;
hi = t + half;
Y(:, lo < flank_s & hi > -flank_s) = NaN;
end


function [spike_times, spike_clusters] = load_pag_spikes(npx_root, animal_code)
spike_times = [];
spike_clusters = [];
folder = fullfile(npx_root, animal_code);
times_file = fullfile(folder, 'spike_times_PAG.npy');
clusters_file = fullfile(folder, 'spike_clusters_PAG.npy');
if exist(times_file, 'file') ~= 2 || exist(clusters_file, 'file') ~= 2
    return
end
spike_times = double(readNPY(times_file)) / 30000;
spike_clusters = double(readNPY(clusters_file));
end
