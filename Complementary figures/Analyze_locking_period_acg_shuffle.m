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
% Bin edges and the saved analysis parameters are loaded from
% locking_period_acg_summary.mat. Smoothing width, the refractory blank,
% and any other setting that summary does not store are read from the
% assignments at the top of Analyze_locking_period_acg.m.
%
% Those same shuffles are averaged across neurons. The mean sum-to-1
% autocorrelogram of each lock group, including delta-locked neurons
% (peak-locked or trough-locked), is smoothed with the same Gaussian and
% scored for a first peak and a second cycle. Shuffle index s of every
% neuron is averaged the same way, which is one draw of the null mean.
% With 10 ms bins and a +/- 1 s axis the running sum is about 10 MB.
% The autocorrelograms are not computed again. A checkpoint that stored
% only prominences has no shuffled curves, so it cannot be resumed.
%
% Run this file, not the unit test, to produce the figures:
%   Analyze_locking_period_acg_shuffle

%% Paths and parameters
data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
run('\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
phase_folder   = fullfile(data_root, 'Analysis results', 'phase locking data');
npx_root       = fullfile(data_root, 'NPX data', 'NPX raw data');

n_shuffles = 200;
resume     = true;

summary_path = fullfile(parsing_folder, 'locking_period_acg_summary.mat');
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = fileparts(which('Analyze_locking_period_acg_shuffle'));
end
script_path = fullfile(script_dir, 'Analyze_locking_period_acg.m');
[analysis, edges, centers] = load_locking_acg_parameters(summary_path, script_path);
alpha = analysis.alpha;
non_entrained_p = analysis.non_entrained_p;
low_percent = analysis.low_percent;
high_percent = analysis.high_percent;
acg_lag = analysis.acg_lag;
smooth_half = analysis.smooth_half;
skip_abs_lag = analysis.skip_abs_lag;
acg_refractory_s = analysis.acg_refractory_s;
areas_to_keep = analysis.areas_to_keep;
bin_w = median(diff(centers));
n_bins = numel(centers);
smooth_n = max(1, round(smooth_half / bin_w));
if abs(skip_abs_lag - 0.02) > 1e-9
    warning('Analyze_locking_period_acg_shuffle:smoothGap', ...
        ['Peak scoring leaves lags within 20 ms unsmoothed. The loaded ' ...
        'skip_abs_lag is %g s.'], skip_abs_lag);
end
disp(['ACG axis from ', summary_path, ': +/- ', num2str(max(abs(centers))), ...
    ' s, bin ', num2str(bin_w), ' s, ', num2str(n_bins), ' bins'])
disp(['Smoothing from ', analysis.smooth_source, ': smooth_half = ', ...
    num2str(smooth_half), ' s, refractory blank = ', num2str(acg_refractory_s), ' s'])

group_order = {'trough', 'peak', 'nonentrained'};
group_title = {'Trough-locked', 'Peak-locked', 'Non-entrained'};
grand_order = {'trough', 'peak', 'delta', 'nonentrained'};
grand_title = {'Trough-locked', 'Peak-locked', 'Delta-locked', 'Non-entrained'};
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
    'skip_abs_lag', skip_abs_lag, ...
    'acg_refractory_s', acg_refractory_s, ...
    'areas_to_keep', {areas_to_keep}, ...
    'grand_average', true, ...
    'null', 'uniform spike times inside each firing run');

grand_bytes = numel(grand_order) * 2 * 2 * n_shuffles * n_bins * 8 * 2;
disp(sprintf(['Grand-average accumulator is %.1f MB (%d groups, %d shuffles, ' ...
    '%d bins). Each neuron is still shuffled once.'], ...
    grand_bytes / 1e6, numel(grand_order), n_shuffles, n_bins))

records = {};
grand = [];
if resume && exist(checkpoint, 'file') == 2
    stored = load(checkpoint, 'shuffle_records', 'parameters', 'grand_average');
    if ~isfield(stored, 'grand_average')
        error('Analyze_locking_period_acg_shuffle:checkpoint', ...
            ['%s stored peak prominences only. The grand-average null is ' ...
            'accumulated while shuffling, so this file cannot be resumed. ' ...
            'Delete it or set resume = false.'], checkpoint);
    end
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
    grand = stored.grand_average;
    dropped = setdiff(grand.included_keys, neuron_keys(records));
    if ~isempty(dropped)
        error('Analyze_locking_period_acg_shuffle:checkpoint', ...
            ['The area filter removed neurons already added to the ' ...
            'grand average in %s. Delete it or set resume = false.'], checkpoint);
    end
    if grand.n_shuffles ~= n_shuffles || grand.n_bins ~= n_bins
        error('Analyze_locking_period_acg_shuffle:checkpoint', ...
            ['%s grand average does not match %d shuffles and %d bins. ' ...
            'Delete it or set resume = false.'], checkpoint, n_shuffles, n_bins);
    end
    disp(['Resuming from ', num2str(numel(records)), ' neurons in ', checkpoint])
else
    disp(['Checkpoint: ', checkpoint])
end
if isempty(grand)
    grand = empty_grand(grand_order, n_shuffles, n_bins);
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
            null_curves = blank_null_curves(n_shuffles, n_bins);
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
                    [peaks, row] = score_pooled_row(sh_lock.low, centers, smooth_n, acg_refractory_s);
                    nulls = store_null_tail(nulls, 'locking', 'low', s, peaks);
                    null_curves.locking.low(s, :) = row;
                    [peaks, row] = score_pooled_row(sh_lock.high, centers, smooth_n, acg_refractory_s);
                    nulls = store_null_tail(nulls, 'locking', 'high', s, peaks);
                    null_curves.locking.high(s, :) = row;
                end
                if power_ok
                    sh_power = pool_locking_period_acgs(acg_s, edges, run_length, power, ...
                        true, low_percent, high_percent);
                    [peaks, row] = score_pooled_row(sh_power.low, centers, smooth_n, acg_refractory_s);
                    nulls = store_null_tail(nulls, 'power', 'low', s, peaks);
                    null_curves.power.low(s, :) = row;
                    [peaks, row] = score_pooled_row(sh_power.high, centers, smooth_n, acg_refractory_s);
                    nulls = store_null_tail(nulls, 'power', 'high', s, peaks);
                    null_curves.power.high(s, :) = row;
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
            obs_curves = blank_obs_curves(n_bins);
            if lock_ok
                [~, obs_curves.locking.low] = score_pooled_row( ...
                    pooled.low, centers, smooth_n, acg_refractory_s);
                [~, obs_curves.locking.high] = score_pooled_row( ...
                    pooled.high, centers, smooth_n, acg_refractory_s);
            end
            if power_ok
                [~, obs_curves.power.low] = score_pooled_row( ...
                    pooled_power.low, centers, smooth_n, acg_refractory_s);
                [~, obs_curves.power.high] = score_pooled_row( ...
                    pooled_power.high, centers, smooth_n, acg_refractory_s);
            end
            grand = add_neuron_to_grand(grand, this_group, this_key, obs_curves, null_curves);
            records{end + 1, 1} = rec; %#ok<SAGROW>
            done_keys{end + 1, 1} = this_key; %#ok<SAGROW>
            shuffle_records = records; %#ok<NASGU>
            grand_average = grand; %#ok<NASGU>
            save(checkpoint, 'shuffle_records', 'parameters', 'grand_average');
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
report_grand_shuffle(grand, centers, smooth_n, grand_order, grand_title);


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


function [peaks, row] = score_pooled_row(acg_row, centers, smooth_n, refractory)
row = normalized_pooled_row(acg_row, centers, refractory);
peaks = acg_delta_peaks(row, centers, smooth_n);
end


function row = normalized_pooled_row(acg_row, centers, refractory)
row = nan(1, numel(centers));
if isempty(acg_row) || ~any(isfinite(acg_row))
    return
end
row = blank_acg_zero_lag(acg_row(:)', centers, refractory);
row = normalize_acg_rows(row);
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
    && isfield(stored.parameters, 'grand_average') && stored.parameters.grand_average ...
    && stored.parameters.n_shuffles == parameters.n_shuffles ...
    && stored.parameters.alpha == parameters.alpha ...
    && stored.parameters.non_entrained_p == parameters.non_entrained_p ...
    && stored.parameters.low_percent == parameters.low_percent ...
    && stored.parameters.high_percent == parameters.high_percent ...
    && stored.parameters.acg_lag == parameters.acg_lag ...
    && stored.parameters.smooth_half == parameters.smooth_half ...
    && stored.parameters.acg_refractory_s == parameters.acg_refractory_s;
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


function [analysis, edges, centers] = load_locking_acg_parameters(summary_path, script_path)
if exist(summary_path, 'file') ~= 2
    error('Analyze_locking_period_acg_shuffle:nosummary', ...
        ['No autocorrelogram summary at %s. Run Analyze_locking_period_acg ' ...
        'first so the shuffle uses that axis and those parameters.'], summary_path);
end
if exist(script_path, 'file') ~= 2
    error('Analyze_locking_period_acg_shuffle:noscript', ...
        'Analyze_locking_period_acg.m was not found next to this script.');
end
loaded = load(summary_path, 'locking_acg');
if ~isfield(loaded, 'locking_acg') || ~isfield(loaded.locking_acg, 'centers') ...
        || ~isfield(loaded.locking_acg, 'bin_edges')
    error('Analyze_locking_period_acg_shuffle:nosummary', ...
        ['%s has no saved bin edges. Run Analyze_locking_period_acg again ' ...
        'so the shuffle can use the same axis.'], summary_path);
end
centers = loaded.locking_acg.centers(:).';
edges = loaded.locking_acg.bin_edges(:).';
if isfield(loaded.locking_acg, 'parameters')
    analysis = loaded.locking_acg.parameters;
else
    analysis = struct();
end
needed = {'alpha', 'non_entrained_p', 'low_percent', 'high_percent', 'acg_lag', ...
    'smooth_half', 'skip_abs_lag', 'acg_refractory_s', 'delta_lag'};
from_script = read_script_assignments(script_path, needed);
had_smooth = isfield(analysis, 'smooth_half') && ~isempty(analysis.smooth_half);
for i = 1:numel(needed)
    name = needed{i};
    if ~isfield(analysis, name) || isempty(analysis.(name))
        analysis.(name) = from_script.(name);
    end
end
if ~isfield(analysis, 'areas_to_keep') || isempty(analysis.areas_to_keep)
    analysis.areas_to_keep = read_script_areas(script_path);
end
if had_smooth
    analysis.smooth_source = summary_path;
else
    analysis.smooth_source = script_path;
end
end


function values = read_script_assignments(script_path, names)
text = fileread(script_path);
values = struct();
for i = 1:numel(names)
    name = names{i};
    token = regexp(text, ['(?m)^', name, '\s*=\s*([^;%]+)'], 'tokens', 'once');
    if isempty(token)
        error('Analyze_locking_period_acg_shuffle:script', ...
            'Analyze_locking_period_acg.m has no assignment for %s.', name);
    end
    values.(name) = str2num(token{1}); %#ok<ST2NM>
    if isempty(values.(name))
        error('Analyze_locking_period_acg_shuffle:script', ...
            'Could not read %s from Analyze_locking_period_acg.m.', name);
    end
end
end


function areas = read_script_areas(script_path)
text = fileread(script_path);
token = regexp(text, '(?m)^areas_to_keep\s*=\s*(\{[^;]+\})', 'tokens', 'once');
if isempty(token)
    error('Analyze_locking_period_acg_shuffle:script', ...
        'Analyze_locking_period_acg.m has no areas_to_keep assignment.');
end
quoted = regexp(token{1}, '''([^'']+)''', 'tokens');
areas = cell(1, numel(quoted));
for i = 1:numel(quoted)
    areas{i} = quoted{i}{1};
end
end


function grand = empty_grand(group_order, n_shuffles, n_bins)
grand.n_shuffles = n_shuffles;
grand.n_bins = n_bins;
grand.included_keys = cell(0, 1);
grand.group_order = group_order;
splits = {'locking', 'power'};
tails = {'low', 'high'};
for s = 1:numel(splits)
    for t = 1:numel(tails)
        for g = 1:numel(group_order)
            block.sum_obs = zeros(1, n_bins);
            block.n_obs = zeros(1, n_bins);
            block.sum_null = zeros(n_shuffles, n_bins);
            block.n_null = zeros(n_shuffles, n_bins);
            block.n_neurons = 0;
            grand.(splits{s}).(tails{t}).(group_order{g}) = block;
        end
    end
end
end


function curves = blank_null_curves(n_shuffles, n_bins)
curves.locking.low = nan(n_shuffles, n_bins);
curves.locking.high = nan(n_shuffles, n_bins);
curves.power.low = nan(n_shuffles, n_bins);
curves.power.high = nan(n_shuffles, n_bins);
end


function curves = blank_obs_curves(n_bins)
curves.locking.low = nan(1, n_bins);
curves.locking.high = nan(1, n_bins);
curves.power.low = nan(1, n_bins);
curves.power.high = nan(1, n_bins);
end


function grand = add_neuron_to_grand(grand, lock_group, neuron_key, obs_curves, null_curves)
targets = {lock_group};
if any(strcmp(lock_group, {'trough', 'peak'}))
    targets{end + 1} = 'delta'; %#ok<AGROW>
end
splits = {'locking', 'power'};
tails = {'low', 'high'};
for g = 1:numel(targets)
    for s = 1:numel(splits)
        for t = 1:numel(tails)
            row = obs_curves.(splits{s}).(tails{t});
            if ~any(isfinite(row))
                continue
            end
            grand.(splits{s}).(tails{t}).(targets{g}) = add_grand_curve( ...
                grand.(splits{s}).(tails{t}).(targets{g}), row, ...
                null_curves.(splits{s}).(tails{t}));
        end
    end
end
grand.included_keys{end + 1, 1} = neuron_key; %#ok<AGROW>
end


function block = add_grand_curve(block, row, null_rows)
finite = isfinite(row);
block.sum_obs(finite) = block.sum_obs(finite) + row(finite);
block.n_obs(finite) = block.n_obs(finite) + 1;
block.n_neurons = block.n_neurons + 1;
for s = 1:size(null_rows, 1)
    sample = null_rows(s, :);
    finite = isfinite(sample);
    block.sum_null(s, finite) = block.sum_null(s, finite) + sample(finite);
    block.n_null(s, finite) = block.n_null(s, finite) + 1;
end
end


function report_grand_shuffle(grand, centers, smooth_n, group_order, group_title)
splits = {'locking', 'power'};
split_labels = {'locking', 'delta power'};
tails = {'low', 'high'};
tail_labels = {'Low 25%', 'High 25%'};
disp(' ')
disp('Grand-average shuffle. Peaks are scored on the Gaussian-smoothed mean sum-to-1 autocorrelogram.')
disp('The all-runs curve is not tested: runs between the two tails are not shuffled.')
for s = 1:numel(splits)
    figure('Name', ['Grand average shuffle, ', split_labels{s}], 'Color', 'w', ...
        'Units', 'normalized', 'OuterPosition', [0.04 0.08 0.92 0.8])
    for t = 1:numel(tails)
        for g = 1:numel(group_order)
            block = grand.(splits{s}).(tails{t}).(group_order{g});
            subplot(numel(tails), numel(group_order), (t - 1) * numel(group_order) + g)
            result = grand_peak_result(block, centers, smooth_n);
            plot_grand_panel(result, [group_title{g}, ', ', tail_labels{t}])
            disp(sprintf(['  %s, %s, %s: n = %d, first prom = %.3g, shuffle p = %.3g, ' ...
                'second prom = %.3g, shuffle p = %.3g'], ...
                group_title{g}, split_labels{s}, tail_labels{t}, result.n, ...
                result.first_prominence, result.first_p, ...
                result.second_prominence, result.second_p))
        end
    end
    sgtitle({['Grand mean versus shuffle, ', split_labels{s}], ...
        'Black is the observed mean. Grey is the central 95% of shuffled means.'})
end
end


function result = grand_peak_result(block, centers, smooth_n)
result.n = block.n_neurons;
result.first_prominence = NaN;
result.first_p = NaN;
result.second_prominence = NaN;
result.second_p = NaN;
result.first_lag = NaN;
result.second_lag = NaN;
result.fold_t = centers(centers > 0);
result.observed = nan(size(result.fold_t));
result.null_lo = result.observed;
result.null_hi = result.observed;
result.null_median = result.observed;
if block.n_neurons < 1
    return
end
mu = block.sum_obs ./ block.n_obs;
mu(block.n_obs == 0) = NaN;
obs = acg_delta_peaks(mu, centers, smooth_n);
result.fold_t = obs.fold_t;
result.observed = obs.fold_y;
result.first_lag = obs.first_lag;
result.second_lag = obs.second_lag;
result.first_prominence = obs.first_prominence;
result.second_prominence = obs.second_prominence;
n_s = size(block.sum_null, 1);
null_fold = nan(n_s, numel(obs.fold_t));
null_first = nan(n_s, 1);
null_second = nan(n_s, 1);
for s = 1:n_s
    sample = block.sum_null(s, :) ./ block.n_null(s, :);
    sample(block.n_null(s, :) == 0) = NaN;
    pk = acg_delta_peaks(sample, centers, smooth_n);
    null_first(s) = null_first_value(pk);
    null_second(s) = null_second_value(pk);
    if numel(pk.fold_y) == numel(obs.fold_t)
        null_fold(s, :) = pk.fold_y;
    end
end
result.first_p = exception_p(obs.first_prominence, null_first);
if obs.second_measurable
    result.second_p = exception_p(obs.second_prominence, null_second);
end
result.null_lo = prctile(null_fold, 2.5, 1);
result.null_hi = prctile(null_fold, 97.5, 1);
result.null_median = median(null_fold, 1, 'omitnan');
end


function plot_grand_panel(result, panel_title)
if result.n < 1
    title({panel_title, 'n = 0'})
    axis off
    return
end
t = result.fold_t;
hold on
fill([t, fliplr(t)], [result.null_lo, fliplr(result.null_hi)], ...
    [0.75 0.75 0.75], 'EdgeColor', 'none', 'HandleVisibility', 'off')
plot(t, result.null_median, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1)
plot(t, result.observed, 'k', 'LineWidth', 1.6)
mark_grand_peak(result.first_lag, result.observed, result.fold_t)
mark_grand_peak(result.second_lag, result.observed, result.fold_t)
xlabel('Lag (s)')
ylabel('Smoothed mean ACG')
title({panel_title, sprintf('n = %d\nfirst p = %.3g, second p = %.3g', ...
    result.n, result.first_p, result.second_p)})
box off
end


function mark_grand_peak(lag, y, t)
if ~isfinite(lag) || isempty(t)
    return
end
[~, idx] = min(abs(t - lag));
plot(t(idx), y(idx), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 5, 'HandleVisibility', 'off')
end
