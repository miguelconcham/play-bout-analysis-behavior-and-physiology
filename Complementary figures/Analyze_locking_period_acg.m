%% Analyze_locking_period_acg
% Mean autocorrelogram of spikes in low and high periods, and in every
% firing run. Reads the delta-lock table and the firing-run tables written
% by Estimate_spike_train_parsing. Does not recompute phase locking.
%
% Each neuron is labeled from EntireSession in delta_all_neurons_v2.mat:
% trough-locked, peak-locked (PPC p < alpha, split at +/- pi/2), or
% non-entrained (PPC p >= 0.1). For each neuron the saved runs are split
% twice, each time at that neuron's own values: bottom 25% and top 25% of
% PPC, and bottom 25% and top 25% of delta power. A neuron is plotted for
% a split only when at least four of its runs have that score, so the two
% tails are separated. The third curve pools every firing run. Delta power
% is the square of the median delta-band envelope stored in lfp_param
% column 2. Squaring does not change which runs fall in the outer
% quartiles. The estimate stores autocorrelograms inside firing runs, so
% spikes that fall in the gaps between runs are not in these curves.
%
% Run autocorrelograms are overlap-weighted (pair counts / (T-|lag|)), the
% lag-0 bins are blanked, and each neuron is scaled to sum to 1 so a cell
% with more spikes does not dominate the mean. One figure per lock group
% for the locking split, and one figure per lock group for the power split.
% Each split also has a figure of all plotted neurons together.
% The lag axis runs to +/- 1 s. Saved histograms that stop at +/- 0.5 s are
% rebuilt from the session spike times and the saved run intervals.
% Shading marks lags of 0.2-0.5 s (2-5 Hz). A 1 Hz side peak sits at the
% edge of this window.
%
% Each plotted autocorrelogram is also scored neuron by neuron. The burst
% is the mean level from the first positive bin through 50 ms. Oscillatory
% peaks are local maxima from 100 ms onward; height is findpeaks
% prominence, the height above the higher neighboring valley. The first
% peak is the earliest one. A second peak counts only near twice that lag,
% and only when that window was observed. Low and high are the same
% neuron, so the plots connect the two conditions with a line. A second
% cycle that the runs are too short to measure is left out of the pair.
% For each measure the script tests whether the high tail is larger than
% the low tail in the same neuron, and whether a bigger change in PPC or
% delta power goes with a bigger change in that measure.
% A reference figure for each split shows the strongest t1-only neuron and
% the strongest t1-and-t2 neuron in each lock group. The printed burst,
% prominences, lags, and mean score are the values on the scatter plots.
% Example neurons are the higher-firing-rate cells among those with that
% peak pattern. In the paired plots, higher-rate neurons are drawn in front.
% The mean curve still gives every neuron the same weight.

%% Paths and parameters
data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
run('\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
phase_folder   = fullfile(data_root, 'Analysis results', 'phase locking data');
npx_root       = fullfile(data_root, 'NPX data', 'NPX raw data');

alpha            = 0.01;
non_entrained_p  = 0.1;
low_percent      = 25;
high_percent     = 75;
delta_lag        = [0.2 0.5];   % s, 5 Hz to 2 Hz
acg_lag          = 1;       % s, half-width of the autocorrelogram axis
smooth_half      = 0.04;        % s, full smooth width; Gaussian sigma is a quarter of this
skip_abs_lag     = 0.02;        % s, do not smooth across the lag-0 gap
acg_refractory_s = 0.005;

group_order = {'trough', 'peak', 'nonentrained'};
group_title = {'Trough-locked', 'Peak-locked', 'Non-entrained'};

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
for g = 1:numel(group_order)
    disp([group_title{g}, ': ', num2str(sum(strcmp(lock_group, group_order{g}))), ...
        ' / ', num2str(n_cells)])
end

%% One parsing file at a time
parsing_files = list_parsing_mat_files(parsing_folder);
if isempty(parsing_files)
    error('Analyze_locking_period_acg:empty', ...
        'No spike-train parsing files in %s. Run Estimate_spike_train_parsing first.', ...
        parsing_folder);
end

loaded_codes = {};
edges = [];
centers = [];
n_bins = [];
bin_w = [];
acg_raw = empty_group_store(group_order, 0);
acg_power = empty_group_store(group_order, 0);
n_in_parsing = zeros(1, numel(group_order));
n_kept = zeros(1, numel(group_order));
n_unsplit = zeros(1, numel(group_order));
n_kept_power = zeros(1, numel(group_order));
n_unsplit_power = zeros(1, numel(group_order));
n_missing = 0;

for f = 1:numel(parsing_files)
    disp(['Loading ', parsing_files(f).name, ' ...'])
    S = load(parsing_files(f).path, 'parsing_struct');
    if ~isfield(S, 'parsing_struct')
        warning('Analyze_locking_period_acg:contents', ...
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
            acg_raw = empty_group_store(group_order, n_bins);
            acg_power = empty_group_store(group_order, n_bins);
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
            warning('Analyze_locking_period_acg:ppc', ...
                'No PPC column in %s; session skipped.', animal_code);
            clear ps
            continue
        end

        session_times = [];
        session_clusters = [];
        if ~saved_covers
            [session_times, session_clusters] = load_pag_spikes(npx_root, animal_code);
            if isempty(session_times)
                warning('Analyze_locking_period_acg:spikes', ...
                    ['Spike times for %s were not found, so the autocorrelogram ' ...
                    'could not be extended to +/- %g s. Session skipped.'], ...
                    animal_code, acg_lag);
                clear ps
                continue
            end
            disp(['  recomputing autocorrelograms to +/- ', num2str(acg_lag), ' s'])
        end

        for nn = 1:height(ci)
            this_id = id_col(nn);
            hit = find(session_matches(all_neurons.session, animal_code) & id_all == this_id, 1);
            if isempty(hit)
                n_missing = n_missing + 1;
                continue
            end
            this_group = lock_group{hit};
            g = find(strcmp(group_order, this_group), 1);
            if isempty(g)
                continue
            end
            n_in_parsing(g) = n_in_parsing(g) + 1;

            rt = run_cell{nn};
            if isempty(rt) || ~istable(rt) || ~ismember('autoccorrelograms', rt.Properties.VariableNames)
                n_unsplit(g) = n_unsplit(g) + 1;
                n_unsplit_power(g) = n_unsplit_power(g) + 1;
                continue
            end
            if size(rt.phase_entrainment, 2) < ppc_col
                n_unsplit(g) = n_unsplit(g) + 1;
                n_unsplit_power(g) = n_unsplit_power(g) + 1;
                continue
            end
            if saved_covers
                acg = rt.autoccorrelograms;
                this_corrected = already_corrected;
            elseif ~ismember('RunStartTime', rt.Properties.VariableNames) ...
                    || ~ismember('RunEndTime', rt.Properties.VariableNames)
                acg = [];
                this_corrected = true;
            else
                acg = run_window_acgs(session_times(session_clusters == this_id), ...
                    rt.RunStartTime, rt.RunEndTime, edges);
                this_corrected = true;
            end
            if size(acg, 2) ~= n_bins || size(acg, 1) ~= height(rt)
                warning('Analyze_locking_period_acg:width', ...
                    'ACG size mismatch in %s id %g; cell skipped.', animal_code, this_id);
                n_unsplit(g) = n_unsplit(g) + 1;
                n_unsplit_power(g) = n_unsplit_power(g) + 1;
                continue
            end
            run_length = acg_run_lengths(rt);
            this_area = cluster_area_char(ci, nn);
            this_rate = neuron_run_rate(rt);
            pooled = pool_locking_period_acgs(acg, edges, run_length, ...
                rt.phase_entrainment(:, ppc_col), this_corrected, low_percent, high_percent);
            if split_is_usable(pooled)
                acg_raw.(group_order{g}) = append_pooled( ...
                    acg_raw.(group_order{g}), pooled, animal_code, this_id, this_area, this_rate);
                n_kept(g) = n_kept(g) + 1;
            else
                n_unsplit(g) = n_unsplit(g) + 1;
            end

            pooled_power = pool_locking_period_acgs(acg, edges, run_length, ...
                run_delta_power(rt), this_corrected, low_percent, high_percent);
            if split_is_usable(pooled_power)
                acg_power.(group_order{g}) = append_pooled( ...
                    acg_power.(group_order{g}), pooled_power, animal_code, this_id, this_area, this_rate);
                n_kept_power(g) = n_kept_power(g) + 1;
            else
                n_unsplit_power(g) = n_unsplit_power(g) + 1;
            end
        end

        loaded_codes{end + 1} = animal_code; %#ok<AGROW>
        disp(['  processed ', animal_code])
        clear ps ci run_cell
    end
    clear S
end

disp(['Parsing sessions: ', num2str(numel(loaded_codes))])
disp(['Cells in the parsing files but not in the delta table: ', num2str(n_missing)])
if isempty(edges)
    error('Analyze_locking_period_acg:nobins', 'No parsing struct contained ACG bin edges.');
end

%% Normalize and plot
smooth_n = max(1, round(smooth_half / bin_w));
curve_colors = [0.20 0.45 0.75; 0.85 0.25 0.20; 0.15 0.15 0.15];
parameters = struct( ...
    'alpha', alpha, ...
    'non_entrained_p', non_entrained_p, ...
    'low_percent', low_percent, ...
    'high_percent', high_percent, ...
    'delta_lag', delta_lag, ...
    'acg_lag', acg_lag, ...
    'delta_power', 'square of median delta-band envelope (lfp_param column 2)');

locking_labels = { ...
    sprintf('Low %g%% locking', low_percent), ...
    sprintf('High %g%% locking', 100 - high_percent), ...
    'All firing runs'};
power_labels = { ...
    sprintf('Low %g%% delta power', low_percent), ...
    sprintf('High %g%% delta power', 100 - high_percent), ...
    'All firing runs'};

locking_acg = summarize_and_plot(acg_raw, group_order, group_title, ...
    locking_labels, 'locking', centers, delta_lag, skip_abs_lag, smooth_n, ...
    acg_refractory_s, curve_colors, n_in_parsing, n_unsplit);
power_acg = summarize_and_plot(acg_power, group_order, group_title, ...
    power_labels, 'delta power', centers, delta_lag, skip_abs_lag, smooth_n, ...
    acg_refractory_s, curve_colors, n_in_parsing, n_unsplit_power);
locking_acg.centers = centers;
locking_acg.bin_edges = edges;
locking_acg.parameters = parameters;
power_acg.centers = centers;
power_acg.bin_edges = edges;
power_acg.parameters = parameters;
locking_acg = add_peak_scores(locking_acg, centers, smooth_n, group_order);
power_acg = add_peak_scores(power_acg, centers, smooth_n, group_order);
plot_paired_peak_scores(locking_acg, group_order, group_title, ...
    'Mean PPC of the runs', 'locking');
plot_paired_peak_scores(power_acg, group_order, group_title, ...
    'Mean delta power of the runs', 'delta power');
plot_peak_examples(locking_acg, centers, smooth_n, group_order, group_title, ...
    'PPC', 'locking');
plot_peak_examples(power_acg, centers, smooth_n, group_order, group_title, ...
    'delta power', 'delta power');

summary_path = fullfile(parsing_folder, 'locking_period_acg_summary.mat');
save(summary_path, 'locking_acg', 'power_acg');
disp(['Saved ', summary_path])


function store = empty_group_store(group_order, n_bins)
store = struct();
for g = 1:numel(group_order)
    block = struct();
    block.n = 0;
    block.low = nan(0, n_bins);
    block.high = nan(0, n_bins);
    block.all = nan(0, n_bins);
    block.n_low = zeros(0, 1);
    block.n_high = zeros(0, 1);
    block.n_all = zeros(0, 1);
    block.score_low = zeros(0, 1);
    block.score_high = zeros(0, 1);
    block.score_all = zeros(0, 1);
    block.firing_rate = zeros(0, 1);
    block.session = {};
    block.cluster_id = zeros(0, 1);
    block.area = {};
    store.(group_order{g}) = block;
end
end


function ok = split_is_usable(pooled)
ok = pooled.n_low >= 1 && pooled.n_high >= 1 ...
    && any(isfinite(pooled.low)) && any(isfinite(pooled.high)) ...
    && any(isfinite(pooled.all));
end


function block = append_pooled(block, pooled, animal_code, this_id, this_area, firing_rate)
row = block.n + 1;
block.n = row;
block.low(row, :) = pooled.low;
block.high(row, :) = pooled.high;
block.all(row, :) = pooled.all;
block.n_low(row, 1) = pooled.n_low;
block.n_high(row, 1) = pooled.n_high;
block.n_all(row, 1) = pooled.n_all;
block.score_low(row, 1) = pooled.score_low;
block.score_high(row, 1) = pooled.score_high;
block.score_all(row, 1) = pooled.score_all;
block.firing_rate(row, 1) = firing_rate;
block.session{row, 1} = animal_code;
block.cluster_id(row, 1) = this_id;
block.area{row, 1} = this_area;
end


function power = run_delta_power(rt)
% Square of the median delta-band envelope. lfp_param is
% [half range, median envelope, mean frequency].
if ~ismember('lfp_param', rt.Properties.VariableNames) || size(rt.lfp_param, 2) < 2
    power = nan(height(rt), 1);
    return
end
power = double(rt.lfp_param(:, 2)) .^ 2;
end


function summary = summarize_and_plot(acg_raw, group_order, group_title, curve_labels, split_name, centers, delta_lag, skip_abs_lag, smooth_n, acg_refractory_s, curve_colors, n_in_parsing, n_unsplit)
curve_names = {'low', 'high', 'all'};
summary = struct();
for g = 1:numel(group_order)
    name = group_order{g};
    n = acg_raw.(name).n;
    disp([group_title{g}, ', ', split_name, ': ', num2str(n_in_parsing(g)), ...
        ' in parsing, ', num2str(n), ' plotted, ', num2str(n_unsplit(g)), ...
        ' skipped (no separated 25/75 split or no spikes)'])
    if n < 1
        warning('Analyze_locking_period_acg:group', ...
            'No %s neurons with both %s tails.', group_title{g}, split_name);
        continue
    end
    keep = 1:n;
    for c = 1:numel(curve_names)
        raw = blank_acg_zero_lag(acg_raw.(name).(curve_names{c})(keep, :), centers, acg_refractory_s);
        summary.(name).(curve_names{c}) = normalize_acg_rows(raw);
    end
    summary.(name).session = acg_raw.(name).session(keep);
    summary.(name).cluster_id = acg_raw.(name).cluster_id(keep);
    summary.(name).area = acg_raw.(name).area(keep);
    summary.(name).n_low = acg_raw.(name).n_low(keep);
    summary.(name).n_high = acg_raw.(name).n_high(keep);
    summary.(name).n_all = acg_raw.(name).n_all(keep);
    summary.(name).score_low = acg_raw.(name).score_low(keep);
    summary.(name).score_high = acg_raw.(name).score_high(keep);
    summary.(name).score_all = acg_raw.(name).score_all(keep);
    summary.(name).firing_rate = acg_raw.(name).firing_rate(keep);
    summary.(name) = plot_acg_curves(summary.(name), curve_names, curve_labels, ...
        [group_title{g}, ', ', split_name], centers, delta_lag, skip_abs_lag, ...
        smooth_n, curve_colors, split_name);
end

summary.all_neurons = stack_plotted_neurons(summary, group_order);
if summary.all_neurons.n < 1
    warning('Analyze_locking_period_acg:all', ...
        'No neurons to pool for the %s figure.', split_name);
    return
end
n_all = summary.all_neurons.n;
disp(['All neurons, ', split_name, ': ', num2str(n_all), ' plotted'])
summary.all_neurons = plot_acg_curves(summary.all_neurons, curve_names, curve_labels, ...
    ['All neurons, ', split_name], centers, delta_lag, skip_abs_lag, ...
    smooth_n, curve_colors, split_name);
end


function combined = stack_plotted_neurons(summary, group_order)
curve_names = {'low', 'high', 'all'};
combined = struct('n', 0);
for c = 1:numel(curve_names)
    combined.(curve_names{c}) = [];
end
combined.session = cell(0, 1);
combined.cluster_id = zeros(0, 1);
combined.area = cell(0, 1);
combined.lock_group = cell(0, 1);
combined.n_low = zeros(0, 1);
combined.n_high = zeros(0, 1);
combined.n_all = zeros(0, 1);
combined.score_low = zeros(0, 1);
combined.score_high = zeros(0, 1);
combined.score_all = zeros(0, 1);
combined.firing_rate = zeros(0, 1);
for g = 1:numel(group_order)
    name = group_order{g};
    if ~isfield(summary, name) || ~isfield(summary.(name), 'low')
        continue
    end
    n = size(summary.(name).low, 1);
    for c = 1:numel(curve_names)
        combined.(curve_names{c}) = [combined.(curve_names{c}); summary.(name).(curve_names{c})]; %#ok<AGROW>
    end
    combined.session = [combined.session; summary.(name).session]; %#ok<AGROW>
    combined.cluster_id = [combined.cluster_id; summary.(name).cluster_id]; %#ok<AGROW>
    combined.area = [combined.area; summary.(name).area]; %#ok<AGROW>
    combined.lock_group = [combined.lock_group; repmat(group_order(g), n, 1)]; %#ok<AGROW>
    combined.n_low = [combined.n_low; summary.(name).n_low]; %#ok<AGROW>
    combined.n_high = [combined.n_high; summary.(name).n_high]; %#ok<AGROW>
    combined.n_all = [combined.n_all; summary.(name).n_all]; %#ok<AGROW>
    combined.score_low = [combined.score_low; summary.(name).score_low]; %#ok<AGROW>
    combined.score_high = [combined.score_high; summary.(name).score_high]; %#ok<AGROW>
    combined.score_all = [combined.score_all; summary.(name).score_all]; %#ok<AGROW>
    combined.firing_rate = [combined.firing_rate; summary.(name).firing_rate]; %#ok<AGROW>
    combined.n = combined.n + n;
end
end


function group_acg = plot_acg_curves(group_acg, curve_names, curve_labels, figure_label, centers, delta_lag, skip_abs_lag, smooth_n, curve_colors, split_name)
figure('Name', figure_label, 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.15 0.15 0.55 0.55])
hold on
shade = shade_delta_lags(delta_lag);
plotted = gobjects(1, numel(curve_names));
for c = 1:numel(curve_names)
    shown = smooth_acg_away_from_zero( ...
        group_acg.(curve_names{c}), centers, skip_abs_lag, smooth_n);
    plotted(c) = plot_mean_sem(centers, shown, curve_colors(c, :));
    mass = delta_side_mass(group_acg.(curve_names{c}), centers, delta_lag);
    disp(sprintf('  %s  %s   mean ACG in [%g %g] s = %.4g   n = %d', ...
        figure_label, curve_labels{c}, delta_lag(1), delta_lag(2), ...
        mean(mass, 'omitnan'), sum(isfinite(mass))))
    group_acg.([curve_names{c}, '_delta_mass']) = mass;
end
report_split_contrast(group_acg, figure_label, split_name);
xline(0, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off')
half_bin = median(diff(centers)) / 2;
xlim([centers(1) - half_bin, centers(end) + half_bin])
xlabel('Lag (s)')
ylabel('Mean normalized autocorrelogram')
title({figure_label, sprintf('n = %d    shaded: 2-5 Hz', size(group_acg.low, 1))})
legend(plotted, curve_labels, 'Location', 'best')
place_delta_shade(shade)
box off
end


function summary = add_peak_scores(summary, centers, smooth_n, group_order)
names = [group_order, {'all_neurons'}];
for k = 1:numel(names)
    if ~isfield(summary, names{k}) || ~isfield(summary.(names{k}), 'low')
        continue
    end
    summary.(names{k}) = score_condition_acgs(summary.(names{k}), centers, smooth_n);
end
end


function group = score_condition_acgs(group, centers, smooth_n)
curve_names = {'low', 'high', 'all'};
for c = 1:numel(curve_names)
    group.([curve_names{c}, '_peaks']) = score_acg_matrix( ...
        group.(curve_names{c}), centers, smooth_n);
end
end


function scores = score_acg_matrix(acg, centers, smooth_n)
n = size(acg, 1);
scores.burst = nan(n, 1);
scores.first_measurable = false(n, 1);
scores.first_present = false(n, 1);
scores.first_lag = nan(n, 1);
scores.first_prominence = nan(n, 1);
scores.second_measurable = false(n, 1);
scores.second_present = false(n, 1);
scores.second_lag = nan(n, 1);
scores.second_prominence = nan(n, 1);
scores.n_cycles = nan(n, 1);
for row = 1:n
    peaks = acg_delta_peaks(acg(row, :), centers, smooth_n);
    scores.burst(row) = peaks.burst;
    scores.first_measurable(row) = peaks.first_measurable;
    scores.first_present(row) = peaks.first_present;
    scores.first_lag(row) = peaks.first_lag;
    scores.first_prominence(row) = peaks.first_prominence;
    scores.second_measurable(row) = peaks.second_measurable;
    scores.second_present(row) = peaks.second_present;
    scores.second_lag(row) = peaks.second_lag;
    scores.second_prominence(row) = peaks.second_prominence;
    scores.n_cycles(row) = peaks.n_cycles;
end
end


function plot_paired_peak_scores(summary, group_order, group_title, score_label, split_name)
if ~isfield(summary, 'all_neurons') || ~isfield(summary.all_neurons, 'low_peaks') ...
        || summary.all_neurons.n < 1
    return
end
neurons = summary.all_neurons;
colors = lock_group_colors();
measures = {'burst', 'first_prominence', 'second_prominence'};
measure_labels = { ...
    'Burst level, 0-50 ms', ...
    'First-peak prominence', ...
    'Second-cycle prominence'};
figure('Name', ['Paired ACG peaks, ', split_name], 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.05 0.05 0.9 0.9])
for m = 1:numel(measures)
    [y_low, y_high, ok] = paired_peak_values(neurons.low_peaks, neurons.high_peaks, measures{m});
    subplot(3, 2, 2 * m - 1)
    plot_paired_categories(neurons.lock_group(ok), y_low(ok), y_high(ok), ...
        neurons.firing_rate(ok), group_order, colors, measure_labels{m})
    change = paired_change_stats(y_low(ok), y_high(ok), ...
        neurons.score_low(ok), neurons.score_high(ok));
    title(sprintf('Paired, n = %d\nmedian high-low %.3g, increase p = %.3g', ...
        change.n, change.median_delta, change.p_increase))
    subplot(3, 2, 2 * m)
    plot_paired_against_score(neurons.lock_group(ok), ...
        neurons.score_low(ok), neurons.score_high(ok), y_low(ok), y_high(ok), ...
        neurons.firing_rate(ok), group_order, colors, score_label, measure_labels{m})
    report_paired_increase(summary, neurons, group_order, group_title, ...
        measures{m}, measure_labels{m}, split_name, score_label);
end
legend_group_colors(group_order, group_title, colors)

figure('Name', ['ACG cycle counts, ', split_name], 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.2 0.2 0.55 0.45])
plot_paired_cycle_counts(neurons, group_order, group_title, colors, split_name)
end


function [y_low, y_high, ok] = paired_peak_values(low_peaks, high_peaks, measure_name)
y_low = low_peaks.(measure_name);
y_high = high_peaks.(measure_name);
ok = isfinite(y_low) & isfinite(y_high);
if strcmp(measure_name, 'second_prominence')
    ok = ok & low_peaks.second_measurable & high_peaks.second_measurable;
elseif strcmp(measure_name, 'first_prominence')
    ok = ok & low_peaks.first_measurable & high_peaks.first_measurable;
end
end


function plot_paired_categories(lock_group, y_low, y_high, firing_rate, group_order, colors, y_label)
hold on
for g = 1:numel(group_order)
    idx = strcmp(lock_group, group_order{g});
    if ~any(idx)
        continue
    end
    order = idx_high_rate_last(idx, firing_rate);
    x = [ones(numel(order), 1), 2 * ones(numel(order), 1)];
    plot(x', [y_low(order), y_high(order)]', '-', ...
        'Color', [colors(g, :), 0.25], 'HandleVisibility', 'off')
    plot([1 2], [mean(y_low(order), 'omitnan'), mean(y_high(order), 'omitnan')], ...
        '-', 'Color', colors(g, :), 'LineWidth', 2.5, 'HandleVisibility', 'off')
end
xlim([0.6 2.4])
set(gca, 'XTick', [1 2], 'XTickLabel', {'Low 25%', 'High 25%'})
ylabel(y_label)
title(sprintf('Paired, n = %d', numel(y_low)))
box off
end


function plot_paired_against_score(lock_group, x_low, x_high, y_low, y_high, firing_rate, group_order, colors, x_label, y_label)
hold on
for g = 1:numel(group_order)
    idx = strcmp(lock_group, group_order{g});
    if ~any(idx)
        continue
    end
    order = idx_high_rate_last(idx, firing_rate);
    plot([x_low(order), x_high(order)]', [y_low(order), y_high(order)]', '-', ...
        'Color', [colors(g, :), 0.3], 'HandleVisibility', 'off')
    scatter(x_low(order), y_low(order), 18, colors(g, :), 'filled', ...
        'MarkerFaceAlpha', 0.7, 'HandleVisibility', 'off')
    scatter(x_high(order), y_high(order), 28, colors(g, :), 'filled', ...
        'MarkerFaceAlpha', 0.9, 'HandleVisibility', 'off')
end
xlabel(x_label)
ylabel(y_label)
title('Each line is one neuron')
box off
end


function report_paired_increase(summary, neurons, group_order, group_title, measure_name, measure_label, split_name, score_label)
disp(['Paired increase, ', split_name, ', ', measure_label])
[y_low, y_high, ok] = paired_peak_values(neurons.low_peaks, neurons.high_peaks, measure_name);
print_paired_increase('All neurons', y_low(ok), y_high(ok), ...
    neurons.score_low(ok), neurons.score_high(ok), score_label);
for g = 1:numel(group_order)
    if ~isfield(summary, group_order{g}) || ~isfield(summary.(group_order{g}), 'low_peaks')
        continue
    end
    block = summary.(group_order{g});
    [g_low, g_high, g_ok] = paired_peak_values(block.low_peaks, block.high_peaks, measure_name);
    print_paired_increase(group_title{g}, g_low(g_ok), g_high(g_ok), ...
        block.score_low(g_ok), block.score_high(g_ok), score_label);
end
end


function print_paired_increase(label, y_low, y_high, x_low, x_high, score_label)
stats = paired_change_stats(y_low, y_high, x_low, x_high);
disp(sprintf(['  %s: median high-low = %.3g, increased %d/%d, ' ...
    'signrank p = %.3g, increase p = %.3g, Spearman with %s change rho = %.3g p = %.3g'], ...
    label, stats.median_delta, stats.n_increased, stats.n, ...
    stats.p_two, stats.p_increase, score_label, stats.rho, stats.p_rho))
end


function stats = paired_change_stats(y_low, y_high, x_low, x_high)
delta = y_high(:) - y_low(:);
score_change = x_high(:) - x_low(:);
stats.n = numel(delta);
stats.n_increased = sum(delta > 0);
if isempty(delta)
    stats.median_delta = NaN;
else
    stats.median_delta = median(delta, 'omitnan');
end
stats.p_two = paired_signrank(y_high, y_low, 'both');
stats.p_increase = paired_signrank(y_high, y_low, 'right');
[stats.rho, stats.p_rho] = spearman_change(score_change, delta);
end


function [rho, p] = spearman_change(dx, dy)
ok = isfinite(dx) & isfinite(dy);
dx = dx(ok);
dy = dy(ok);
if numel(dx) < 3 || numel(unique(dx)) < 2 || numel(unique(dy)) < 2
    rho = NaN;
    p = NaN;
    return
end
[rho, p] = corr(dx(:), dy(:), 'Type', 'Spearman');
end


function plot_paired_cycle_counts(neurons, group_order, group_title, ~, split_name)
classifiable_low = neurons.low_peaks.n_cycles == 0 | neurons.low_peaks.second_measurable;
classifiable_high = neurons.high_peaks.n_cycles == 0 | neurons.high_peaks.second_measurable;
ok = classifiable_low & classifiable_high;
fractions = zeros(numel(group_order), 3, 2);
ns = zeros(numel(group_order), 1);
for g = 1:numel(group_order)
    idx = ok & strcmp(neurons.lock_group, group_order{g});
    ns(g) = sum(idx);
    if ~any(idx)
        continue
    end
    fractions(g, :, 1) = cycle_fraction(neurons.low_peaks.n_cycles(idx));
    fractions(g, :, 2) = cycle_fraction(neurons.high_peaks.n_cycles(idx));
end
low_bars = squeeze(fractions(:, :, 1));
high_bars = squeeze(fractions(:, :, 2));
if numel(group_order) == 1
    low_bars = low_bars(:)';
    high_bars = high_bars(:)';
end
bar_x = (1:numel(group_order))';
cycle_colors = [0.75 0.75 0.75; 0.30 0.55 0.85; 0.80 0.30 0.25];
hold on
low_bars_handle = bar(bar_x - 0.16, low_bars, 0.3, 'stacked');
high_bars_handle = bar(bar_x + 0.16, high_bars, 0.3, 'stacked');
for s = 1:3
    low_bars_handle(s).FaceColor = cycle_colors(s, :);
    low_bars_handle(s).FaceAlpha = 0.45;
    low_bars_handle(s).EdgeColor = 'none';
    high_bars_handle(s).FaceColor = cycle_colors(s, :);
    high_bars_handle(s).EdgeColor = 'none';
end
legend(high_bars_handle, {'0 cycles', '1 cycle', '2 cycles'}, 'Location', 'eastoutside')
xlim([0.4, numel(group_order) + 0.6])
set(gca, 'XTick', 1:numel(group_order), 'XTickLabel', group_title)
ylabel('Fraction of paired neurons')
title({['Cycles in the autocorrelogram, ', split_name], ...
    'Left bar low 25%, right bar high 25%. Stack: 0, 1, 2 cycles'})
for g = 1:numel(group_order)
    text(bar_x(g), 1.02, sprintf('n=%d', ns(g)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8)
end
ylim([0 1.12])
box off
end


function fraction = cycle_fraction(n_cycles)
fraction = [mean(n_cycles == 0), mean(n_cycles == 1), mean(n_cycles == 2)];
end


function colors = lock_group_colors()
colors = [0.20 0.45 0.80; 0.85 0.25 0.20; 0.45 0.45 0.45];
end


function legend_group_colors(group_order, group_title, colors)
hold on
handles = gobjects(1, numel(group_order));
for g = 1:numel(group_order)
    handles(g) = plot(nan, nan, '-', 'Color', colors(g, :), 'LineWidth', 2.5);
end
legend(handles, group_title, 'Location', 'eastoutside')
end


function plot_peak_examples(summary, centers, smooth_n, group_order, group_title, score_label, split_name)
if ~isfield(summary, 'all_neurons') || ~isfield(summary.all_neurons, 'low') ...
        || summary.all_neurons.n < 1
    return
end
neurons = summary.all_neurons;
t1_rows = select_strong_examples(neurons, group_order, 't1');
both_rows = select_strong_examples(neurons, group_order, 'both');
if isempty(t1_rows) && isempty(both_rows)
    disp(['No peak examples for ', split_name])
    return
end
n_col = max(numel(t1_rows), numel(both_rows));
figure('Name', ['ACG peak examples, ', split_name], 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.04 0.05 0.92 0.88])
row_rows = {t1_rows, both_rows};
row_names = {'Strong t1, no t2', 'Strong t1 and t2'};
for r = 1:2
    rows = row_rows{r};
    for c = 1:n_col
        subplot(2, n_col, (r - 1) * n_col + c)
        if c > numel(rows)
            axis off
            continue
        end
        draw_peak_example(neurons, rows(c), centers, smooth_n, ...
            group_order, group_title, score_label, row_names{r}, split_name)
    end
end
sgtitle({['Reference for the peak scatters, ', split_name], ...
    'Higher firing rate is preferred. Grey band is the 0-50 ms burst. The stem is the prominence on the scatter plots.'})
end


function order = idx_high_rate_last(idx, firing_rate)
positions = find(idx);
rate = firing_rate(positions);
rate(~isfinite(rate)) = -inf;
[~, sort_idx] = sort(rate, 'ascend');
order = positions(sort_idx);
end


function rows = select_strong_examples(neurons, group_order, kind)
strength = example_strength(neurons, kind);
if isfield(neurons, 'firing_rate')
    rate = neurons.firing_rate(:);
else
    rate = nan(numel(strength), 1);
end
rows = zeros(0, 1);
for g = 1:numel(group_order)
    idx = find(strcmp(neurons.lock_group, group_order{g}) & strength > 0);
    if isempty(idx)
        continue
    end
    rows(end + 1, 1) = pick_high_rate_example(idx, strength, rate); %#ok<AGROW>
end
end


function pick = pick_high_rate_example(idx, strength, rate)
peak = strength(idx);
strong = peak >= median(peak);
if ~any(strong)
    strong = true(size(peak));
end
candidates = idx(strong);
candidate_rate = rate(candidates);
if ~any(isfinite(candidate_rate))
    [~, local] = max(strength(candidates));
    pick = candidates(local);
    return
end
candidate_rate(~isfinite(candidate_rate)) = -inf;
[~, local] = max(candidate_rate);
pick = candidates(local);
end


function strength = example_strength(neurons, kind)
n = neurons.n;
strength = -inf(n, 1);
for i = 1:n
    if strcmp(kind, 't1') && ...
            (neurons.low_peaks.second_present(i) || neurons.high_peaks.second_present(i))
        continue
    end
    strength(i) = max( ...
        side_strength(neurons.low_peaks, i, kind), ...
        side_strength(neurons.high_peaks, i, kind));
end
end


function s = side_strength(peaks, i, kind)
s = -inf;
if ~peaks.first_present(i) || ~isfinite(peaks.first_prominence(i)) ...
        || peaks.first_prominence(i) <= 0
    return
end
if strcmp(kind, 't1')
    if peaks.second_measurable(i) && ~peaks.second_present(i)
        s = peaks.first_prominence(i);
    end
    return
end
if peaks.second_present(i) && isfinite(peaks.second_prominence(i)) ...
        && peaks.second_prominence(i) > 0
    s = min(peaks.first_prominence(i), peaks.second_prominence(i));
end
end


function draw_peak_example(neurons, row, centers, smooth_n, group_order, group_title, score_label, kind_name, split_name)
low_color = [0.20 0.45 0.75];
high_color = [0.85 0.25 0.20];
low_pk = acg_delta_peaks(neurons.low(row, :), centers, smooth_n);
high_pk = acg_delta_peaks(neurons.high(row, :), centers, smooth_n);
hold on
low_line = plot(low_pk.fold_t, low_pk.fold_y, 'Color', low_color, 'LineWidth', 1.6);
high_line = plot(high_pk.fold_t, high_pk.fold_y, 'Color', high_color, 'LineWidth', 1.6);
mark_scored_peak(low_pk, 'first', low_color);
mark_scored_peak(low_pk, 'second', low_color);
mark_scored_peak(high_pk, 'first', high_color);
mark_scored_peak(high_pk, 'second', high_color);
xlim([0 1])
xlabel('Lag (s)')
ylabel('Folded normalized ACG')
group_name = neurons.lock_group{row};
hit = find(strcmp(group_order, group_name), 1);
if ~isempty(hit)
    group_name = group_title{hit};
end
title({kind_name, sprintf('%s, %s, id %g, %s, %s', ...
    group_name, neurons.session{row}, neurons.cluster_id(row), neurons.area{row}, ...
    rate_label(neurons, row))}, ...
    'Interpreter', 'none')
yl = ylim;
burst = patch([0 0.05 0.05 0], [yl(1) yl(1) yl(2) yl(2)], [0.75 0.75 0.75], ...
    'FaceAlpha', 0.28, 'EdgeColor', 'none', 'HandleVisibility', 'off');
ylim(yl)
uistack(burst, 'bottom')
if ~isempty(low_line) && ~isempty(high_line)
    legend([low_line, high_line], {'Low 25%', 'High 25%'}, 'Location', 'southeast')
end
text(0.98, 0.98, example_parameter_text(neurons, row, score_label), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 3, 'Interpreter', 'none')
box off
disp(sprintf('  example %s, %s: %s id %g\n%s', split_name, kind_name, ...
    neurons.session{row}, neurons.cluster_id(row), ...
    example_parameter_text(neurons, row, score_label)))
end


function mark_scored_peak(peaks, which, color)
if strcmp(which, 'first')
    present = peaks.first_present;
    lag = peaks.first_lag;
    prom = peaks.first_prominence;
    label = 't1';
else
    present = peaks.second_present;
    lag = peaks.second_lag;
    prom = peaks.second_prominence;
    label = 't2';
end
if ~present || ~isfinite(lag) || ~isfinite(prom) || isempty(peaks.fold_t)
    return
end
y = interp1(peaks.fold_t, peaks.fold_y, lag, 'nearest');
if ~isfinite(y)
    return
end
plot([lag lag], [max(y - prom, 0), y], '-', 'Color', color, ...
    'LineWidth', 1.3, 'HandleVisibility', 'off')
plot(lag, y, 'o', 'Color', color, 'MarkerFaceColor', color, ...
    'MarkerSize', 5, 'HandleVisibility', 'off')
text(lag, y, ['  ' label], 'Color', color, 'FontSize', 8, ...
    'VerticalAlignment', 'bottom', 'Interpreter', 'none', 'HandleVisibility', 'off')
end


function txt = example_parameter_text(neurons, row, score_label)
txt = sprintf('%s\n%s', ...
    tail_parameter_line('Low', neurons.low_peaks, row, neurons.score_low(row), score_label), ...
    tail_parameter_line('High', neurons.high_peaks, row, neurons.score_high(row), score_label));
end


function line = tail_parameter_line(name, peaks, row, score, score_label)
line = sprintf('%s  burst %.3g   t1 %s   t2 %s   %s %.3g', ...
    name, peaks.burst(row), ...
    peak_phrase(peaks.first_present(row), peaks.first_measurable(row), ...
        peaks.first_lag(row), peaks.first_prominence(row)), ...
    peak_phrase(peaks.second_present(row), peaks.second_measurable(row), ...
        peaks.second_lag(row), peaks.second_prominence(row)), ...
    score_label, score);
end


function phrase = peak_phrase(present, measurable, lag, prom)
if present && isfinite(lag) && isfinite(prom)
    phrase = sprintf('%d ms, prom %.3g', round(1000 * lag), prom);
    return
end
if measurable
    phrase = 'none, prom 0';
    return
end
phrase = 'n/a';
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


function rate = neuron_run_rate(rt)
% Spikes per second across firing runs. RunNumSpikes is the count stored
% with each run.
if ~ismember('RunNumSpikes', rt.Properties.VariableNames)
    rate = NaN;
    return
end
n_spikes = double(rt.RunNumSpikes(:));
duration = acg_run_lengths(rt);
ok = isfinite(n_spikes) & n_spikes >= 0 & isfinite(duration) & duration > 0;
if ~any(ok)
    rate = NaN;
    return
end
rate = sum(n_spikes(ok)) / sum(duration(ok));
end


function label = rate_label(neurons, row)
if ~isfield(neurons, 'firing_rate') || ~isfinite(neurons.firing_rate(row))
    label = 'rate n/a';
    return
end
label = sprintf('%.3g sp/s', neurons.firing_rate(row));
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
    error('Analyze_locking_period_acg:id', 'cluster table has neither cluster_id nor id');
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


function Y = smooth_acg_away_from_zero(X, centers, skip_abs_lag, smooth_n)
Y = X;
if isempty(X)
    return
end
pos = centers > skip_abs_lag;
neg = centers < -skip_abs_lag;
for j = 1:size(X, 1)
    if any(pos)
        Y(j, pos) = gaussian_smooth_acg(X(j, pos), smooth_n / 4);
    end
    if any(neg)
        Y(j, neg) = gaussian_smooth_acg(X(j, neg), smooth_n / 4);
    end
end
end


function mass = delta_side_mass(X, centers, lag_win)
lag_idx = abs(centers) >= lag_win(1) & abs(centers) <= lag_win(2);
mass = mean(X(:, lag_idx), 2, 'omitnan');
end


function report_split_contrast(group_acg, group_label, split_name)
p_low = paired_signrank(group_acg.high_delta_mass, group_acg.low_delta_mass);
p_all = paired_signrank(group_acg.high_delta_mass, group_acg.all_delta_mass);
disp(sprintf('  %s signrank, high vs low 25%% %s: p = %.3g', group_label, split_name, p_low))
disp(sprintf('  %s signrank, high vs all runs: p = %.3g', group_label, p_all))
end


function p = paired_signrank(a, b, tail)
if nargin < 3 || isempty(tail)
    tail = 'both';
end
ok = isfinite(a) & isfinite(b);
a = a(ok);
b = b(ok);
if numel(a) < 2 || all(a == b)
    p = NaN;
    return
end
p = signrank(a, b, 'tail', tail);
end


function h = plot_mean_sem(t, M, col)
mu = mean(M, 1, 'omitnan');
n = sum(isfinite(M), 1);
sem = std(M, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
fill_sem_finite(t, mu, sem, col)
h = plot(t, mu, 'Color', col, 'LineWidth', 2);
end


function fill_sem_finite(t, mu, sem, col)
lo = mu - sem;
hi = mu + sem;
ok = isfinite(mu) & isfinite(lo) & isfinite(hi);
d = diff([false, ok, false]);
i1 = find(d == 1);
i2 = find(d == -1) - 1;
for k = 1:numel(i1)
    sl = i1(k):i2(k);
    if numel(sl) < 2
        continue
    end
    fill([t(sl), fliplr(t(sl))], [lo(sl), fliplr(hi(sl))], col, ...
        'FaceAlpha', 0.22, 'EdgeColor', 'none', 'HandleVisibility', 'off')
end
end


function shade = shade_delta_lags(lag_win)
shade = gobjects(1, 2);
for k = 1:2
    sign = 2 * (k - 1) - 1;
    x1 = sign * lag_win(1);
    x2 = sign * lag_win(2);
    shade(k) = patch([x1 x2 x2 x1], [0 0 1 1], [0.95 0.85 0.40], ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end


function place_delta_shade(shade)
yl = ylim;
for k = 1:numel(shade)
    x = get(shade(k), 'XData');
    set(shade(k), 'YData', [yl(1) yl(1) yl(2) yl(2)], 'XData', x)
end
uistack(shade, 'bottom')
end
