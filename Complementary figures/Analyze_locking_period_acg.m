%% Analyze_locking_period_acg
% Mean autocorrelogram of spikes in low and high periods, and in every
% firing run. Reads the delta-lock table and the firing-run tables written
% by Estimate_spike_train_parsing. Does not recompute phase locking.
%
% Each neuron is labeled from EntireSession in delta_all_neurons_v2.mat:
% trough-locked, peak-locked (PPC p < alpha, split at +/- pi/2), or
% non-entrained (PPC p >= 0.1). Delta-locked is the peak-locked and
% trough-locked neurons together. For each neuron the saved runs are split
% twice, each time at that neuron's own values: bottom 25% and top 25% of
% PPC, and bottom 25% and top 25% of delta power. A neuron is plotted for
% a split only when at least four of its runs have that score, so the two
% tails are separated. The third curve pools every firing run. Delta power
% is the square of the median delta-band envelope stored in lfp_param
% column 2. Squaring does not change which runs fall in the outer
% quartiles. The estimate stores autocorrelograms inside firing runs, so
% spikes that fall in the gaps between runs are not in these curves.
%
% Only areas_to_keep are analyzed. The depth tables label superior
% colliculus as SupCol; that label is kept as SC. Matching ignores case
% and spaces. VLPAG and DLPAG are not treated as LPAG, and DMPAG, InfCol,
% isRt, and PAG are left out.
%
% Run autocorrelograms are overlap-weighted (pair counts / (T-|lag|)), the
% lag-0 bins are blanked, and each neuron is scaled to sum to 1 so a cell
% with more spikes does not dominate the mean. One figure per lock group
% for the locking split, and one figure per lock group for the power split.
% Each split also has a figure of all plotted neurons together.
% The lag axis runs to +/- 1 s. Saved histograms that stop at +/- 0.5 s are
% rebuilt from the session spike times and the saved run intervals.
% Shading marks lags of 0.2-0.5 s (2-5 Hz). A 1 Hz side peak sits at the
% edge of this window. The shading does not set the vertical scale. Every
% mean autocorrelogram uses one y limit, taken from the highest and lowest
% mean plus or minus its standard error, and rounded onto a short decimal.
% On each mean autocorrelogram, grey bars in the top tenth of the axis mark
% clusters that survive the sign-flip test. The bars are mirrored because
% the autocorrelogram is. The statistics report lists every cluster, with
% its positive-lag range and cluster p-value, including those that do not
% survive. A cluster survives when that p-value is at or below acg_cluster_alpha.
%
% Each plotted autocorrelogram is also scored neuron by neuron. The burst
% is the mean level from the first positive bin through 50 ms. Oscillatory
% peaks are local maxima from 100 ms onward; height is findpeaks
% prominence, the height above the higher neighboring valley. The first
% peak is the earliest one. A second peak counts only near twice that lag,
% and only when that window was observed. Low and high are the same
% neuron. The low-versus-high panel is a scatter of the low tail against
% the high tail, with the identity line as reference. The mean is a cross,
% so an increase sits above that line. The panel against PPC or delta power
% still connects the two tails of each neuron. A second
% cycle that the runs are too short to measure is left out of the pair.
% The sum-to-1 scaling is not a division by firing rate. Set
% plot_rate_scaled_peaks to also score the same peaks on
% (pair counts / overlap) / (firing rate of those runs ^ 2 * bin width).
% Each curve uses the firing rate of its own tail. A flat Poisson
% autocorrelogram then sits at 1, so a change in rate does not lift the
% curve. The same switch also redraws the mean autocorrelograms in that
% ratio, with their own y limit and the same low-versus-high cluster test.
% The same all-runs autocorrelogram is then used in a linear model of
% mean delta power. Predictors are the mean firing rate of those runs,
% their mean PPC, the 0-50 ms burst, and the first-peak and second-cycle
% prominences.
% Relative weights give each predictor's share of the explained variance.
% For each measure the script tests whether the high tail is larger than
% the low tail in the same neuron, and whether a bigger change in PPC or
% delta power goes with a bigger change in that measure. The same paired
% plots add a third column for the latency of the first peak and of the
% second cycle, and test whether that lag shifts between the two tails.
% A separate figure shows the distribution of those first-peak latencies
% in the low tail and the high tail, using the same paired neurons.
% A reference figure for each split shows a t1-only neuron and a
% t1-and-t2 neuron in each lock group. Among cells with that peak pattern,
% the higher firing rate is chosen. The printed burst, prominences, lags,
% and mean score are the values on the scatter plots.
% Example neurons are the higher-firing-rate cells among those with that
% peak pattern. In the paired plots, higher-rate neurons are drawn in front.
% The mean curve still gives every neuron the same weight.
%
% Analyze_locking_period_acg_shuffle is a separate, slower script. It
% rebuilds the same peaks after placing spikes uniformly inside each run,
% and it can be stopped and continued.
% figure_folder, at the top of the script, is where figures are written
% as SVG and JPG and where statistics.txt collects the command-window
% tests. Leave it empty to save nothing. Set reprint_stats_only to reload
% the saved summary and print those tests again without reloading spikes
% or opening figures. Set redraw_figures_only to reload that summary and
% draw the figures, including cluster bars computed from the saved curves.
% Spikes are not reloaded.

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
acg_cluster_alpha = 0.05;  % uncorrected paired threshold, and the cluster cutoff
acg_cluster_perm  = 1000;  % sign flips for the cluster p-value
plot_rate_scaled_peaks = true;  % false: sum-to-1 autocorrelograms only
                                 % true: also the rate-ratio autocorrelograms (Poisson = 1), their cluster tests, and the peak figures
areas_to_keep = {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'};
% Leave empty to keep figures on screen only. Otherwise every figure is
% written there as SVG and JPG, and statistics.txt gets the command window.
figure_folder = '';
% true reloads locking_period_acg_summary.mat and reprints the tests.
% Spikes are not reloaded and no figures are opened. A non-empty
% figure_folder still receives a new statistics.txt.
reprint_stats_only = false;
% true reloads the same summary and draws the figures. Spike files are
% not read. Cluster bars are computed from the saved curves when the
% summary does not already contain them.
redraw_figures_only = false;

source_order = {'trough', 'peak', 'nonentrained'};
source_title = {'Trough-locked', 'Peak-locked', 'Non-entrained'};
group_order = {'trough', 'peak', 'delta', 'nonentrained'};
group_title = {'Trough-locked', 'Peak-locked', 'Delta-locked', 'Non-entrained'};

if exist('diary_cleanup', 'var')
    delete(diary_cleanup)
    clear diary_cleanup
end
start_statistics_diary(figure_folder)

if reprint_stats_only
    summary_path = fullfile(parsing_folder, 'locking_period_acg_summary.mat');
    if exist(summary_path, 'file') ~= 2
        error('Analyze_locking_period_acg:nosummary', ...
            'No saved summary at %s. Run once with reprint_stats_only = false.', summary_path);
    end
    loaded = load(summary_path, 'locking_acg', 'power_acg');
    smooth_n = max(1, round(smooth_half / median(diff(loaded.locking_acg.centers))));
    reprint_saved_statistics(loaded.locking_acg, loaded.power_acg, ...
        group_order, group_title, plot_rate_scaled_peaks, smooth_n, ...
        skip_abs_lag, acg_cluster_alpha, acg_cluster_perm);
    disp(['Reprinted statistics from ', summary_path])
    stop_statistics_diary()
    return
end

if redraw_figures_only
    summary_path = fullfile(parsing_folder, 'locking_period_acg_summary.mat');
    if exist(summary_path, 'file') ~= 2
        error('Analyze_locking_period_acg:nosummary', ...
            'No saved summary at %s. Run once with redraw_figures_only = false.', summary_path);
    end
    loaded = load(summary_path, 'locking_acg', 'power_acg');
    smooth_n = max(1, round(smooth_half / median(diff(loaded.locking_acg.centers))));
    [locking_acg, power_acg] = redraw_saved_figures(loaded.locking_acg, loaded.power_acg, ...
        group_order, group_title, plot_rate_scaled_peaks, smooth_n, ...
        skip_abs_lag, delta_lag, acg_cluster_alpha, acg_cluster_perm);
    save(summary_path, 'locking_acg', 'power_acg');
    disp(['Drew figures from ', summary_path])
    if ~isempty(figure_folder)
        save_open_figures(figure_folder)
    end
    stop_statistics_diary()
    return
end

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
in_area = false(n_cells, 1);
for i = 1:n_cells
    in_area(i) = ~isempty(canonical_kept_area(cluster_area_char(all_neurons, i), areas_to_keep));
end
disp_block('Lock groups')
for g = 1:numel(source_order)
    disp([source_title{g}, ': ', num2str(sum(strcmp(lock_group, source_order{g}) & in_area)), ...
        ' / ', num2str(sum(in_area)), ' in ', strjoin(areas_to_keep, ', ')])
end
disp(['Delta-locked: ', num2str(sum(ismember(lock_group, {'trough', 'peak'}) & in_area)), ...
    ' / ', num2str(sum(in_area)), ' in ', strjoin(areas_to_keep, ', '), ...
    ' (peak-locked or trough-locked)'])

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
acg_raw = empty_group_store(source_order, 0);
acg_power = empty_group_store(source_order, 0);
n_in_parsing = zeros(1, numel(source_order));
n_kept = zeros(1, numel(source_order));
n_unsplit = zeros(1, numel(source_order));
n_kept_power = zeros(1, numel(source_order));
n_unsplit_power = zeros(1, numel(source_order));
n_missing = 0;
n_outside = 0;

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
            acg_raw = empty_group_store(source_order, n_bins);
            acg_power = empty_group_store(source_order, n_bins);
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
            this_area = cluster_area_char(ci, nn);
            if isempty(strtrim(this_area))
                this_area = cluster_area_char(all_neurons, hit);
            end
            this_area = canonical_kept_area(this_area, areas_to_keep);
            if isempty(this_area)
                n_outside = n_outside + 1;
                continue
            end
            this_group = lock_group{hit};
            g = find(strcmp(source_order, this_group), 1);
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
            this_rate = neuron_run_rate(rt);
            this_ppc = rt.phase_entrainment(:, ppc_col);
            this_power = run_delta_power(rt);
            run_means = all_run_means(rt, acg, run_length, this_ppc, this_power);
            pooled = pool_locking_period_acgs(acg, edges, run_length, ...
                this_ppc, this_corrected, low_percent, high_percent);
            if split_is_usable(pooled)
                acg_raw.(source_order{g}) = append_pooled( ...
                    acg_raw.(source_order{g}), pooled, animal_code, this_id, this_area, ...
                    this_rate, tail_rate_struct(rt, pooled), run_means);
                n_kept(g) = n_kept(g) + 1;
            else
                n_unsplit(g) = n_unsplit(g) + 1;
            end

            pooled_power = pool_locking_period_acgs(acg, edges, run_length, ...
                this_power, this_corrected, low_percent, high_percent);
            if split_is_usable(pooled_power)
                acg_power.(source_order{g}) = append_pooled( ...
                    acg_power.(source_order{g}), pooled_power, animal_code, this_id, this_area, ...
                    this_rate, tail_rate_struct(rt, pooled_power), run_means);
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

disp_block('Parsing')
disp(['Parsing sessions: ', num2str(numel(loaded_codes))])
disp(['Cells in the parsing files but not in the delta table: ', num2str(n_missing)])
disp(['Cells outside ', strjoin(areas_to_keep, ', '), ': ', num2str(n_outside)])
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
    'plot_rate_scaled_peaks', plot_rate_scaled_peaks, ...
    'areas_to_keep', {areas_to_keep}, ...
    'figure_folder', figure_folder, ...
    'acg_cluster_alpha', acg_cluster_alpha, ...
    'acg_cluster_perm', acg_cluster_perm, ...
    'delta_power', 'square of median delta-band envelope (lfp_param column 2)');

locking_labels = { ...
    sprintf('Low %g%% locking', low_percent), ...
    sprintf('High %g%% locking', 100 - high_percent), ...
    'All firing runs'};
power_labels = { ...
    sprintf('Low %g%% delta power', low_percent), ...
    sprintf('High %g%% delta power', 100 - high_percent), ...
    'All firing runs'};

locking_acg = summarize_and_plot(acg_raw, source_order, source_title, ...
    locking_labels, 'locking', centers, delta_lag, skip_abs_lag, smooth_n, ...
    acg_refractory_s, curve_colors, n_in_parsing, n_unsplit, ...
    acg_cluster_alpha, acg_cluster_perm);
power_acg = summarize_and_plot(acg_power, source_order, source_title, ...
    power_labels, 'delta power', centers, delta_lag, skip_abs_lag, smooth_n, ...
    acg_refractory_s, curve_colors, n_in_parsing, n_unsplit_power, ...
    acg_cluster_alpha, acg_cluster_perm);
acg_ylim = shared_mean_acg_ylim(locking_acg, power_acg, group_order, ...
    centers, skip_abs_lag, smooth_n);
parameters.acg_ylim = acg_ylim;
locking_acg.centers = centers;
locking_acg.bin_edges = edges;
locking_acg.parameters = parameters;
power_acg.centers = centers;
power_acg.bin_edges = edges;
power_acg.parameters = parameters;
locking_acg = plot_summary_acgs(locking_acg, group_order, group_title, ...
    locking_labels, 'locking', centers, delta_lag, skip_abs_lag, smooth_n, ...
    curve_colors, acg_ylim);
power_acg = plot_summary_acgs(power_acg, group_order, group_title, ...
    power_labels, 'delta power', centers, delta_lag, skip_abs_lag, smooth_n, ...
    curve_colors, acg_ylim);
if plot_rate_scaled_peaks
    plot_rate_scaled_mean_acgs(locking_acg, power_acg, bin_w, centers, ...
        group_order, group_title, locking_labels, power_labels, delta_lag, ...
        skip_abs_lag, smooth_n, curve_colors, acg_cluster_alpha, acg_cluster_perm);
end
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
if plot_rate_scaled_peaks
    plot_rate_scaled_peaks_for(locking_acg, bin_w, centers, smooth_n, ...
        group_order, group_title, 'Mean PPC of the runs', 'locking');
    plot_rate_scaled_peaks_for(power_acg, bin_w, centers, smooth_n, ...
        group_order, group_title, 'Mean delta power of the runs', 'delta power');
end
locking_acg.delta_power_model = fit_delta_power_contributions( ...
    locking_acg, group_order, group_title);

summary_path = fullfile(parsing_folder, 'locking_period_acg_summary.mat');
save(summary_path, 'locking_acg', 'power_acg');
disp(['Saved ', summary_path])
if ~isempty(figure_folder)
    save_open_figures(figure_folder)
end
stop_statistics_diary()


function start_statistics_diary(folder)
% A script keeps variables in the base workspace. An onCleanup left there
% runs when the next run replaces it, and that call turns diary off before
% any line is printed, so statistics.txt stays empty.
if strcmp(get(0, 'Diary'), 'on')
    diary('off')
end
if isempty(folder)
    return
end
if exist(folder, 'dir') ~= 7
    mkdir(folder)
end
stats_path = fullfile(folder, 'statistics.txt');
if exist(stats_path, 'file') == 2
    delete(stats_path)
end
diary(stats_path)
end


function stop_statistics_diary()
if strcmp(get(0, 'Diary'), 'on')
    diary('off')
end
end


function save_open_figures(folder)
figs = findall(0, 'Type', 'figure');
if isempty(figs)
    disp(['No figures to save in ', folder])
    return
end
[~, order] = sort([figs.Number]);
figs = figs(order);
for i = 1:numel(figs)
    name = strtrim(figs(i).Name);
    if isempty(name)
        name = sprintf('figure_%d', figs(i).Number);
    end
    name = regexprep(name, '[^\w]+', '_');
    name = regexprep(name, '^_|_$', '');
    if isempty(name)
        name = sprintf('figure_%d', figs(i).Number);
    end
    base = name;
    suffix = 1;
    while exist(fullfile(folder, [base '.svg']), 'file') == 2
        suffix = suffix + 1;
        base = sprintf('%s_%d', name, suffix);
    end
    set(figs(i), 'PaperPositionMode', 'auto')
    svg_path = fullfile(folder, [base '.svg']);
    jpg_path = fullfile(folder, [base '.jpg']);
    print(figs(i), svg_path, '-dsvg')
    print(figs(i), jpg_path, '-djpeg', '-r200')
    disp(['Saved ', svg_path])
    disp(['Saved ', jpg_path])
end
end


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
    block.rate_low = zeros(0, 1);
    block.rate_high = zeros(0, 1);
    block.rate_all = zeros(0, 1);
    block.model_fr = zeros(0, 1);
    block.model_ppc = zeros(0, 1);
    block.model_power = zeros(0, 1);
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


function block = append_pooled(block, pooled, animal_code, this_id, this_area, firing_rate, tail_rates, run_means)
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
block.rate_low(row, 1) = tail_rates.low;
block.rate_high(row, 1) = tail_rates.high;
block.rate_all(row, 1) = tail_rates.all;
block.model_fr(row, 1) = run_means.fr;
block.model_ppc(row, 1) = run_means.ppc;
block.model_power(row, 1) = run_means.power;
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


function summary = summarize_and_plot(acg_raw, group_order, group_title, curve_labels, split_name, centers, delta_lag, skip_abs_lag, smooth_n, acg_refractory_s, curve_colors, n_in_parsing, n_unsplit, cluster_alpha, cluster_perm)
curve_names = {'low', 'high', 'all'};
summary = struct();
disp_block(['Mean autocorrelogram, ', split_name])
for g = 1:numel(group_order)
    name = group_order{g};
    n = acg_raw.(name).n;
    disp_rule()
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
        summary.(name).(['mass_', curve_names{c}]) = sum(raw, 2, 'omitnan');
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
    summary.(name).rate_low = acg_raw.(name).rate_low(keep);
    summary.(name).rate_high = acg_raw.(name).rate_high(keep);
    summary.(name).rate_all = acg_raw.(name).rate_all(keep);
    summary.(name).model_fr = acg_raw.(name).model_fr(keep);
    summary.(name).model_ppc = acg_raw.(name).model_ppc(keep);
    summary.(name).model_power = acg_raw.(name).model_power(keep);
    summary.(name) = record_delta_mass(summary.(name), curve_names, curve_labels, ...
        [group_title{g}, ', ', split_name], centers, delta_lag, split_name, ...
        skip_abs_lag, smooth_n, cluster_alpha, cluster_perm);
end

summary = add_delta_locked_group(summary);
if isfield(summary, 'delta') && isfield(summary.delta, 'low') && ~isempty(summary.delta.low)
    disp_rule()
    disp(['Delta-locked, ', split_name, ': ', num2str(size(summary.delta.low, 1)), ...
        ' plotted (peak-locked or trough-locked)'])
    summary.delta = record_delta_mass(summary.delta, curve_names, curve_labels, ...
        ['Delta-locked, ', split_name], centers, delta_lag, split_name, ...
        skip_abs_lag, smooth_n, cluster_alpha, cluster_perm);
end

summary.all_neurons = stack_plotted_neurons(summary, group_order);
if summary.all_neurons.n < 1
    warning('Analyze_locking_period_acg:all', ...
        'No neurons to pool for the %s figure.', split_name);
    return
end
n_all = summary.all_neurons.n;
disp_rule()
disp(['All neurons, ', split_name, ': ', num2str(n_all), ' plotted'])
summary.all_neurons = record_delta_mass(summary.all_neurons, curve_names, curve_labels, ...
    ['All neurons, ', split_name], centers, delta_lag, split_name, ...
    skip_abs_lag, smooth_n, cluster_alpha, cluster_perm);
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
combined.rate_low = zeros(0, 1);
combined.rate_high = zeros(0, 1);
combined.rate_all = zeros(0, 1);
combined.model_fr = zeros(0, 1);
combined.model_ppc = zeros(0, 1);
combined.model_power = zeros(0, 1);
combined.mass_low = zeros(0, 1);
combined.mass_high = zeros(0, 1);
combined.mass_all = zeros(0, 1);
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
    combined.rate_low = [combined.rate_low; summary.(name).rate_low]; %#ok<AGROW>
    combined.rate_high = [combined.rate_high; summary.(name).rate_high]; %#ok<AGROW>
    combined.rate_all = [combined.rate_all; summary.(name).rate_all]; %#ok<AGROW>
    combined.model_fr = [combined.model_fr; summary.(name).model_fr]; %#ok<AGROW>
    combined.model_ppc = [combined.model_ppc; summary.(name).model_ppc]; %#ok<AGROW>
    combined.model_power = [combined.model_power; summary.(name).model_power]; %#ok<AGROW>
    combined.mass_low = [combined.mass_low; summary.(name).mass_low]; %#ok<AGROW>
    combined.mass_high = [combined.mass_high; summary.(name).mass_high]; %#ok<AGROW>
    combined.mass_all = [combined.mass_all; summary.(name).mass_all]; %#ok<AGROW>
    combined.n = combined.n + n;
end
end


function summary = add_delta_locked_group(summary)
% Peak-locked and trough-locked neurons as one extra condition.
has_trough = isfield(summary, 'trough') && isfield(summary.trough, 'low') && ~isempty(summary.trough.low);
has_peak = isfield(summary, 'peak') && isfield(summary.peak, 'low') && ~isempty(summary.peak.low);
if has_trough && has_peak
    summary.delta = concat_neuron_blocks(summary.trough, summary.peak);
elseif has_trough
    summary.delta = summary.trough;
elseif has_peak
    summary.delta = summary.peak;
else
    return
end
drop = {'low_high_clusters', 'low_delta_mass', 'high_delta_mass', 'all_delta_mass'};
for k = 1:numel(drop)
    if isfield(summary.delta, drop{k})
        summary.delta = rmfield(summary.delta, drop{k});
    end
end
end


function group = concat_neuron_blocks(a, b)
fn = fieldnames(a);
group = struct();
for i = 1:numel(fn)
    name = fn{i};
    if ~isfield(b, name)
        continue
    end
    if strcmp(name, 'low_high_clusters') || endsWith(name, '_delta_mass')
        continue
    end
    va = a.(name);
    vb = b.(name);
    if isstruct(va)
        group.(name) = concat_struct_fields(va, vb);
    elseif isnumeric(va) || islogical(va) || iscell(va)
        group.(name) = [va; vb]; %#ok<AGROW>
    end
end
end


function out = concat_struct_fields(a, b)
out = struct();
fn = fieldnames(a);
for i = 1:numel(fn)
    name = fn{i};
    if ~isfield(b, name)
        continue
    end
    va = a.(name);
    vb = b.(name);
    if isnumeric(va) || islogical(va)
        out.(name) = [va; vb]; %#ok<AGROW>
    end
end
end


function summary = plot_summary_acgs(summary, group_order, group_title, curve_labels, split_name, centers, delta_lag, skip_abs_lag, smooth_n, curve_colors, acg_ylim, y_label, name_suffix)
if nargin < 12 || isempty(y_label)
    y_label = 'Mean normalized autocorrelogram';
end
if nargin < 13
    name_suffix = '';
end
curve_names = {'low', 'high', 'all'};
for g = 1:numel(group_order)
    name = group_order{g};
    if ~isfield(summary, name) || ~isfield(summary.(name), 'low')
        continue
    end
    summary.(name) = plot_acg_curves(summary.(name), curve_names, curve_labels, ...
        [group_title{g}, ', ', split_name, name_suffix], centers, delta_lag, skip_abs_lag, ...
        smooth_n, curve_colors, acg_ylim, y_label);
end
if ~isfield(summary, 'all_neurons') || ~isfield(summary.all_neurons, 'low') ...
        || summary.all_neurons.n < 1
    return
end
summary.all_neurons = plot_acg_curves(summary.all_neurons, curve_names, curve_labels, ...
    ['All neurons, ', split_name, name_suffix], centers, delta_lag, skip_abs_lag, ...
    smooth_n, curve_colors, acg_ylim, y_label);
end


function group_acg = record_delta_mass(group_acg, curve_names, curve_labels, figure_label, centers, delta_lag, split_name, skip_abs_lag, smooth_n, cluster_alpha, cluster_perm)
for c = 1:numel(curve_names)
    mass = delta_side_mass(group_acg.(curve_names{c}), centers, delta_lag);
    disp(sprintf('  %s  %s   mean ACG in [%g %g] s = %.4g   n = %d', ...
        figure_label, curve_labels{c}, delta_lag(1), delta_lag(2), ...
        mean(mass, 'omitnan'), sum(isfinite(mass))))
    group_acg.([curve_names{c}, '_delta_mass']) = mass;
end
report_split_contrast(group_acg, figure_label, split_name)
low_s = smooth_acg_away_from_zero(group_acg.low, centers, skip_abs_lag, smooth_n);
high_s = smooth_acg_away_from_zero(group_acg.high, centers, skip_abs_lag, smooth_n);
group_acg.low_high_clusters = paired_acg_clusters( ...
    low_s, high_s, centers, cluster_alpha, cluster_perm);
report_acg_clusters(group_acg.low_high_clusters, figure_label, split_name)
end


function group_acg = plot_acg_curves(group_acg, curve_names, curve_labels, figure_label, centers, delta_lag, skip_abs_lag, smooth_n, curve_colors, acg_ylim, y_label)
if nargin < 11 || isempty(y_label)
    y_label = 'Mean normalized autocorrelogram';
end
figure('Name', figure_label, 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.15 0.15 0.55 0.55])
hold on
shade = shade_delta_lags(delta_lag);
plotted = gobjects(1, numel(curve_names));
for c = 1:numel(curve_names)
    shown = smooth_acg_away_from_zero( ...
        group_acg.(curve_names{c}), centers, skip_abs_lag, smooth_n);
    plotted(c) = plot_mean_sem(centers, shown, curve_colors(c, :));
end
xline(0, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off')
half_bin = median(diff(centers)) / 2;
xlim([centers(1) - half_bin, centers(end) + half_bin])
ylim(acg_ylim)
yticks(readable_acg_ticks(acg_ylim))
ytickformat('%g')
xlabel('Lag (s)')
ylabel(y_label)
title({figure_label, sprintf('n = %d    band: 2-5 Hz    top bars: clusters that survive', size(group_acg.low, 1))})
legend(plotted, curve_labels, 'Location', 'best')
place_delta_shade(shade)
if isfield(group_acg, 'low_high_clusters')
    shade_acg_clusters(group_acg.low_high_clusters, acg_ylim)
end
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


function plot_paired_peak_scores(summary, group_order, group_title, score_label, split_name, unit_label)
if nargin < 6 || isempty(unit_label)
    unit_label = '';
end
if ~isfield(summary, 'all_neurons') || ~isfield(summary.all_neurons, 'low_peaks') ...
        || summary.all_neurons.n < 1
    return
end
neurons = summary.all_neurons;
colors = lock_group_colors(group_order);
point_keep = ~strcmp(group_order, 'delta');
point_order = group_order(point_keep);
point_title = group_title(point_keep);
point_colors = colors(point_keep, :);
measures = {'burst', 'first_prominence', 'second_prominence'};
measure_labels = { ...
    'Burst level, 0-50 ms', ...
    'First-peak prominence', ...
    'Second-cycle prominence'};
if isempty(unit_label)
    scale_label = 'sum-to-1 normalized autocorrelogram';
    figure_name = ['Paired ACG peaks, ', split_name, ', sum-to-1'];
else
    scale_label = 'rate ratio (Poisson = 1)';
    for m = 1:numel(measure_labels)
        measure_labels{m} = [measure_labels{m}, ' (ratio)'];
    end
    figure_name = ['Paired ACG peaks, ', split_name, ', rate ratio'];
end
report_peak_tests(summary, group_order, group_title, score_label, split_name, unit_label);
n_col = 2 + isempty(unit_label);
figure('Name', figure_name, 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.02 0.05 0.96 0.9])
sgtitle(scale_label)
lag_measures = {'', 'first_lag', 'second_lag'};
lag_labels = {'', 'First-peak latency (s)', 'Second-cycle latency (s)'};
for m = 1:numel(measures)
    [y_low, y_high, ok] = paired_peak_values(neurons.low_peaks, neurons.high_peaks, measures{m});
    subplot(3, n_col, n_col * (m - 1) + 1)
    plot_paired_categories(neurons.lock_group(ok), y_low(ok), y_high(ok), ...
        neurons.firing_rate(ok), point_order, point_colors, measure_labels{m})
    change = paired_change_stats(y_low(ok), y_high(ok), ...
        neurons.score_low(ok), neurons.score_high(ok));
    title(sprintf(['All neurons, n = %d\nmedian high-low %.3g\n' ...
        'two-sided signrank p = %.3g\none-sided increase p = %.3g'], ...
        change.n, change.median_delta, change.p_two, change.p_increase))
    subplot(3, n_col, n_col * (m - 1) + 2)
    plot_paired_against_score(neurons.lock_group(ok), ...
        neurons.score_low(ok), neurons.score_high(ok), y_low(ok), y_high(ok), ...
        neurons.firing_rate(ok), point_order, point_colors, score_label, measure_labels{m})
    title(sprintf('All neurons\nSpearman rho = %.3g, two-sided p = %.3g', ...
        change.rho, change.p_rho))
    if n_col == 3 && m > 1
        [lag_low, lag_high, lag_ok] = paired_peak_values( ...
            neurons.low_peaks, neurons.high_peaks, lag_measures{m});
        subplot(3, n_col, n_col * m)
        plot_paired_categories(neurons.lock_group(lag_ok), lag_low(lag_ok), lag_high(lag_ok), ...
            neurons.firing_rate(lag_ok), point_order, point_colors, lag_labels{m})
        lag_change = paired_change_stats(lag_low(lag_ok), lag_high(lag_ok), ...
            neurons.score_low(lag_ok), neurons.score_high(lag_ok));
        title(sprintf(['All neurons, n = %d\nmedian high-low %.3g s\n' ...
            'two-sided signrank p = %.3g\none-sided later p = %.3g'], ...
            lag_change.n, lag_change.median_delta, lag_change.p_two, lag_change.p_increase))
    end
end
legend_group_colors(point_order, point_title, point_colors)
if isempty(unit_label)
    scale_tag = 'sum-to-1';
else
    scale_tag = 'rate ratio';
end
plot_paired_by_lock_group(neurons, group_order, group_title, colors, split_name, scale_label, scale_tag)

if isempty(unit_label)
    figure('Name', ['ACG cycle counts, ', split_name], 'Color', 'w', ...
        'Units', 'normalized', 'OuterPosition', [0.2 0.2 0.55 0.45])
    plot_paired_cycle_counts(neurons, group_order, group_title, colors, split_name)
    plot_first_latency_distributions(neurons, group_order, group_title, colors, split_name)
end
end


function plot_paired_by_lock_group(neurons, group_order, group_title, colors, split_name, scale_label, scale_tag)
row_measures = {'burst', 'first_prominence', 'second_prominence', 'first_lag', 'second_lag'};
row_labels = { ...
    'Burst level, 0-50 ms', ...
    'First-peak prominence', ...
    'Second-cycle prominence', ...
    'First-peak latency (s)', ...
    'Second-cycle latency (s)'};
side_names = { ...
    'one-sided increase p', 'one-sided increase p', 'one-sided increase p', ...
    'one-sided later p', 'one-sided later p'};
units = {'', '', '', ' s', ' s'};
panel_names = [group_order, {'all'}];
panel_titles = [group_title, {'All neurons'}];
for c = 1:numel(panel_names)
    figure('Name', ['Paired, ', panel_titles{c}, ', ', split_name, ', ', scale_tag], ...
        'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.18 0.04 0.42 0.92])
    if strcmp(panel_names{c}, 'all')
        in_panel = true(numel(neurons.lock_group), 1);
    else
        in_panel = in_lock_group(neurons.lock_group, panel_names{c});
    end
    for r = 1:numel(row_measures)
        [y_low, y_high, ok] = paired_peak_values( ...
            neurons.low_peaks, neurons.high_peaks, row_measures{r});
        in = ok & in_panel;
        subplot(numel(row_measures), 1, r)
        lim = paired_cloud_limits(y_low(in), y_high(in));
        plot_one_paired_cloud(neurons.lock_group(in), y_low(in), y_high(in), ...
            neurons.firing_rate(in), group_order, colors, lim, true)
        ylabel({'High 25%', row_labels{r}})
        change = paired_change_stats(y_low(in), y_high(in), ...
            neurons.score_low(in), neurons.score_high(in));
        title(sprintf(['n = %d, median high-low %.3g%s\n' ...
            'two-sided signrank p = %.3g, %s = %.3g'], ...
            change.n, change.median_delta, units{r}, ...
            change.p_two, side_names{r}, change.p_increase))
    end
    sgtitle({[panel_titles{c}, ', ', split_name], scale_label})
end
end


function lim = paired_cloud_limits(y_low, y_high)
vals = [y_low(:); y_high(:)];
vals = vals(isfinite(vals));
if isempty(vals)
    lim = [0 1];
    return
end
lo = min(vals);
hi = max(vals);
span = hi - lo;
if span <= 0
    span = max(abs(lo), 1);
end
lim = [lo - 0.08 * span, hi + 0.08 * span];
end


function plot_one_paired_cloud(lock_group, y_low, y_high, firing_rate, group_order, colors, lim, show_xlabel)
hold on
plot(lim, lim, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1, 'HandleVisibility', 'off')
for g = 1:numel(group_order)
    idx = strcmp(lock_group, group_order{g});
    if ~any(idx)
        continue
    end
    order = idx_high_rate_last(idx, firing_rate);
    scatter(y_low(order), y_high(order), 36, colors(g, :), 'filled', ...
        'MarkerFaceAlpha', 0.75, 'HandleVisibility', 'off')
end
if ~isempty(y_low)
    plot(mean(y_low, 'omitnan'), mean(y_high, 'omitnan'), '+', ...
        'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2, 'HandleVisibility', 'off')
end
xlim(lim)
ylim(lim)
axis square
if show_xlabel
    xlabel('Low 25%')
end
box off
end


function [y_low, y_high, ok] = paired_peak_values(low_peaks, high_peaks, measure_name)
y_low = low_peaks.(measure_name);
y_high = high_peaks.(measure_name);
ok = isfinite(y_low) & isfinite(y_high);
if strcmp(measure_name, 'second_prominence')
    ok = ok & low_peaks.second_measurable & high_peaks.second_measurable;
elseif strcmp(measure_name, 'first_prominence')
    ok = ok & low_peaks.first_measurable & high_peaks.first_measurable;
elseif strcmp(measure_name, 'first_lag')
    ok = ok & low_peaks.first_present & high_peaks.first_present;
elseif strcmp(measure_name, 'second_lag')
    ok = ok & low_peaks.second_present & high_peaks.second_present ...
        & low_peaks.second_measurable & high_peaks.second_measurable;
end
end


function plot_paired_categories(lock_group, y_low, y_high, firing_rate, group_order, colors, y_label)
% Low tail on x, high tail on y. The dashed line is equality. A cross above
% that line is a larger value in the high tail.
hold on
vals = [y_low(:); y_high(:)];
vals = vals(isfinite(vals));
if isempty(vals)
    axis off
    return
end
lo = min(vals);
hi = max(vals);
span = hi - lo;
if span <= 0
    span = max(abs(lo), 1);
end
pad = 0.08 * span;
lim = [lo - pad, hi + pad];
plot(lim, lim, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1, 'HandleVisibility', 'off')
for g = 1:numel(group_order)
    idx = strcmp(lock_group, group_order{g});
    if ~any(idx)
        continue
    end
    order = idx_high_rate_last(idx, firing_rate);
    scatter(y_low(order), y_high(order), 22, colors(g, :), 'filled', ...
        'MarkerFaceAlpha', 0.75, 'HandleVisibility', 'off')
    plot(mean(y_low(order), 'omitnan'), mean(y_high(order), 'omitnan'), '+', ...
        'Color', colors(g, :), 'MarkerSize', 12, 'LineWidth', 1.6, 'HandleVisibility', 'off')
end
plot(mean(y_low, 'omitnan'), mean(y_high, 'omitnan'), '+', ...
    'Color', 'k', 'MarkerSize', 16, 'LineWidth', 2.4, 'HandleVisibility', 'off')
xlim(lim)
ylim(lim)
axis square
xlabel({'Low 25%', y_label})
ylabel({'High 25%', y_label})
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
box off
end


function report_peak_tests(summary, group_order, group_title, score_label, split_name, unit_label)
if nargin < 6 || isempty(unit_label)
    unit_label = '';
    scale_label = 'sum-to-1 normalized autocorrelogram';
else
    scale_label = 'rate ratio (Poisson = 1)';
end
if ~isfield(summary, 'all_neurons') || ~isfield(summary.all_neurons, 'low_peaks')
    return
end
neurons = summary.all_neurons;
measures = {'burst', 'first_prominence', 'second_prominence'};
measure_labels = { ...
    'Burst level, 0-50 ms', ...
    'First-peak prominence', ...
    'Second-cycle prominence'};
lag_measures = {'first_lag', 'second_lag'};
lag_labels = {'First-peak latency (s)', 'Second-cycle latency (s)'};
if ~isempty(unit_label)
    for m = 1:numel(measure_labels)
        measure_labels{m} = [measure_labels{m}, ' (ratio)'];
    end
end
disp_block(['Paired tests, ', split_name, ': ', scale_label])
for m = 1:numel(measures)
    report_paired_increase(summary, neurons, group_order, group_title, ...
        measures{m}, measure_labels{m}, split_name, score_label, scale_label);
    if isempty(unit_label) && m > 1
        report_paired_latency(summary, neurons, group_order, group_title, ...
            lag_measures{m - 1}, lag_labels{m - 1}, split_name, score_label, scale_label);
    end
end
end


function [locking_acg, power_acg] = redraw_saved_figures(locking_acg, power_acg, group_order, group_title, plot_rate_scaled_peaks, smooth_n, skip_abs_lag, delta_lag, cluster_alpha, cluster_perm)
centers = locking_acg.centers(:).';
bin_w = median(diff(centers));
if isfield(locking_acg, 'parameters') && isfield(locking_acg.parameters, 'low_percent')
    low_percent = locking_acg.parameters.low_percent;
    high_percent = locking_acg.parameters.high_percent;
else
    low_percent = 25;
    high_percent = 75;
end
locking_labels = { ...
    sprintf('Low %g%% locking', low_percent), ...
    sprintf('High %g%% locking', 100 - high_percent), ...
    'All firing runs'};
power_labels = { ...
    sprintf('Low %g%% delta power', low_percent), ...
    sprintf('High %g%% delta power', 100 - high_percent), ...
    'All firing runs'};
curve_colors = [0.20 0.45 0.75; 0.85 0.25 0.20; 0.15 0.15 0.15];
locking_acg = add_delta_locked_group(locking_acg);
power_acg = add_delta_locked_group(power_acg);
locking_acg = ensure_acg_clusters(locking_acg, group_order, group_title, ...
    'locking', centers, skip_abs_lag, smooth_n, cluster_alpha, cluster_perm);
power_acg = ensure_acg_clusters(power_acg, group_order, group_title, ...
    'delta power', centers, skip_abs_lag, smooth_n, cluster_alpha, cluster_perm);
if isfield(locking_acg, 'parameters') && isfield(locking_acg.parameters, 'acg_ylim')
    acg_ylim = locking_acg.parameters.acg_ylim;
else
    acg_ylim = shared_mean_acg_ylim(locking_acg, power_acg, group_order, ...
        centers, skip_abs_lag, smooth_n);
end
locking_acg = plot_summary_acgs(locking_acg, group_order, group_title, ...
    locking_labels, 'locking', centers, delta_lag, skip_abs_lag, smooth_n, ...
    curve_colors, acg_ylim);
power_acg = plot_summary_acgs(power_acg, group_order, group_title, ...
    power_labels, 'delta power', centers, delta_lag, skip_abs_lag, smooth_n, ...
    curve_colors, acg_ylim);
if plot_rate_scaled_peaks
    plot_rate_scaled_mean_acgs(locking_acg, power_acg, bin_w, centers, ...
        group_order, group_title, locking_labels, power_labels, delta_lag, ...
        skip_abs_lag, smooth_n, curve_colors, cluster_alpha, cluster_perm);
end
if ~isfield(locking_acg, 'all_neurons') || ~isfield(locking_acg.all_neurons, 'low_peaks')
    locking_acg = add_peak_scores(locking_acg, centers, smooth_n, group_order);
    power_acg = add_peak_scores(power_acg, centers, smooth_n, group_order);
elseif isfield(locking_acg, 'delta') && ~isfield(locking_acg.delta, 'low_peaks')
    locking_acg.delta = score_condition_acgs(locking_acg.delta, centers, smooth_n);
    if isfield(power_acg, 'delta')
        power_acg.delta = score_condition_acgs(power_acg.delta, centers, smooth_n);
    end
end
plot_paired_peak_scores(locking_acg, group_order, group_title, ...
    'Mean PPC of the runs', 'locking');
plot_paired_peak_scores(power_acg, group_order, group_title, ...
    'Mean delta power of the runs', 'delta power');
plot_peak_examples(locking_acg, centers, smooth_n, group_order, group_title, ...
    'PPC', 'locking');
plot_peak_examples(power_acg, centers, smooth_n, group_order, group_title, ...
    'delta power', 'delta power');
if plot_rate_scaled_peaks
    plot_rate_scaled_peaks_for(locking_acg, bin_w, centers, smooth_n, ...
        group_order, group_title, 'Mean PPC of the runs', 'locking');
    plot_rate_scaled_peaks_for(power_acg, bin_w, centers, smooth_n, ...
        group_order, group_title, 'Mean delta power of the runs', 'delta power');
end
model_has_delta = isfield(locking_acg, 'delta_power_model') ...
    && isfield(locking_acg.delta_power_model, 'groups') ...
    && isfield(locking_acg.delta_power_model.groups, 'delta');
if isfield(locking_acg, 'delta_power_model') && isfield(locking_acg.delta_power_model, 'all') ...
        && locking_acg.delta_power_model.all.ok && model_has_delta
    plot_power_contributions(locking_acg.delta_power_model, group_order, group_title, ...
        locking_acg.delta_power_model.all.names);
else
    locking_acg.delta_power_model = fit_delta_power_contributions( ...
        locking_acg, group_order, group_title);
end
end


function summary = ensure_acg_clusters(summary, group_order, group_title, split_name, centers, skip_abs_lag, smooth_n, cluster_alpha, cluster_perm)
names = [group_order, {'all_neurons'}];
titles = [group_title, {'All neurons'}];
disp_block(['ACG clusters, ', split_name])
for g = 1:numel(names)
    if ~isfield(summary, names{g}) || ~isfield(summary.(names{g}), 'low')
        continue
    end
    group = summary.(names{g});
    if ~isfield(group, 'low_high_clusters')
        low_s = smooth_acg_away_from_zero(group.low, centers, skip_abs_lag, smooth_n);
        high_s = smooth_acg_away_from_zero(group.high, centers, skip_abs_lag, smooth_n);
        group.low_high_clusters = paired_acg_clusters( ...
            low_s, high_s, centers, cluster_alpha, cluster_perm);
    end
    disp_rule()
    report_acg_clusters(group.low_high_clusters, [titles{g}, ', ', split_name], split_name)
    summary.(names{g}) = group;
end
end


function reprint_saved_statistics(locking_acg, power_acg, group_order, group_title, plot_rate_scaled_peaks, smooth_n, skip_abs_lag, cluster_alpha, cluster_perm)
locking_acg = add_delta_locked_group(locking_acg);
power_acg = add_delta_locked_group(power_acg);
params = locking_acg.parameters;
centers = locking_acg.centers(:).';
bin_w = median(diff(centers));
locking_labels = { ...
    sprintf('Low %g%% locking', params.low_percent), ...
    sprintf('High %g%% locking', 100 - params.high_percent), ...
    'All firing runs'};
power_labels = { ...
    sprintf('Low %g%% delta power', params.low_percent), ...
    sprintf('High %g%% delta power', 100 - params.high_percent), ...
    'All firing runs'};
reprint_split_mass(locking_acg, group_order, group_title, locking_labels, ...
    'locking', centers, params.delta_lag, skip_abs_lag, smooth_n, cluster_alpha, cluster_perm);
reprint_split_mass(power_acg, group_order, group_title, power_labels, ...
    'delta power', centers, params.delta_lag, skip_abs_lag, smooth_n, cluster_alpha, cluster_perm);
if plot_rate_scaled_peaks
    rate_scaled_acg_stats(locking_acg, bin_w, group_order, group_title, ...
        locking_labels, 'locking', centers, params.delta_lag, skip_abs_lag, ...
        smooth_n, cluster_alpha, cluster_perm);
    rate_scaled_acg_stats(power_acg, bin_w, group_order, group_title, ...
        power_labels, 'delta power', centers, params.delta_lag, skip_abs_lag, ...
        smooth_n, cluster_alpha, cluster_perm);
end
report_peak_tests(locking_acg, group_order, group_title, ...
    'Mean PPC of the runs', 'locking', '');
report_peak_tests(power_acg, group_order, group_title, ...
    'Mean delta power of the runs', 'delta power', '');
if plot_rate_scaled_peaks
    locking_hz = add_peak_scores(rate_scale_summary(locking_acg, bin_w, group_order), ...
        centers, smooth_n, group_order);
    power_hz = add_peak_scores(rate_scale_summary(power_acg, bin_w, group_order), ...
        centers, smooth_n, group_order);
    report_peak_tests(locking_hz, group_order, group_title, ...
        'Mean PPC of the runs', 'locking', 'ratio');
    report_peak_tests(power_hz, group_order, group_title, ...
        'Mean delta power of the runs', 'delta power', 'ratio');
end
if isfield(locking_acg, 'delta_power_model')
    model = locking_acg.delta_power_model;
    if ~isfield(model, 'groups') || ~isfield(model.groups, 'delta')
        model = fit_delta_power_contributions(locking_acg, group_order, group_title, false);
    end
    if model.all.ok
        report_power_model(model.all);
    end
    for g = 1:numel(group_order)
        if isfield(model, 'groups') && isfield(model.groups, group_order{g}) ...
                && model.groups.(group_order{g}).ok
            report_power_model(model.groups.(group_order{g}));
        end
    end
end
end


function reprint_split_mass(summary, group_order, group_title, curve_labels, split_name, centers, delta_lag, skip_abs_lag, smooth_n, cluster_alpha, cluster_perm)
curve_names = {'low', 'high', 'all'};
names = [group_order, {'all_neurons'}];
titles = [group_title, {'All neurons'}];
disp_block(['Mean autocorrelogram, ', split_name])
for g = 1:numel(names)
    if ~isfield(summary, names{g}) || ~isfield(summary.(names{g}), 'low')
        continue
    end
    disp_rule()
    record_delta_mass(summary.(names{g}), curve_names, curve_labels, ...
        [titles{g}, ', ', split_name], centers, delta_lag, split_name, ...
        skip_abs_lag, smooth_n, cluster_alpha, cluster_perm);
end
end


function disp_block(title)
disp_rule()
disp(title)
end


function disp_rule()
disp(repmat('-', 1, 72))
end


function report_paired_increase(summary, neurons, group_order, group_title, measure_name, measure_label, split_name, score_label, scale_label)
disp_rule()
disp(['Paired increase, ', scale_label, ', ', split_name, ', ', measure_label])
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


function report_paired_latency(summary, neurons, group_order, group_title, measure_name, measure_label, split_name, score_label, scale_label)
disp_rule()
disp(['Paired latency, ', scale_label, ', ', split_name, ', ', measure_label])
[y_low, y_high, ok] = paired_peak_values(neurons.low_peaks, neurons.high_peaks, measure_name);
print_paired_shift('All neurons', y_low(ok), y_high(ok), ...
    neurons.score_low(ok), neurons.score_high(ok), score_label);
for g = 1:numel(group_order)
    if ~isfield(summary, group_order{g}) || ~isfield(summary.(group_order{g}), 'low_peaks')
        continue
    end
    block = summary.(group_order{g});
    [g_low, g_high, g_ok] = paired_peak_values(block.low_peaks, block.high_peaks, measure_name);
    print_paired_shift(group_title{g}, g_low(g_ok), g_high(g_ok), ...
        block.score_low(g_ok), block.score_high(g_ok), score_label);
end
end


function print_paired_shift(label, y_low, y_high, x_low, x_high, score_label)
stats = paired_change_stats(y_low, y_high, x_low, x_high);
disp(sprintf(['  %s: median high-low = %.3g s, later %d/%d, ' ...
    'two-sided signrank p = %.3g, one-sided later p = %.3g, ' ...
    'Spearman with %s change rho = %.3g, two-sided p = %.3g'], ...
    label, stats.median_delta, stats.n_increased, stats.n, ...
    stats.p_two, stats.p_increase, score_label, stats.rho, stats.p_rho))
end


function print_paired_increase(label, y_low, y_high, x_low, x_high, score_label)
stats = paired_change_stats(y_low, y_high, x_low, x_high);
disp(sprintf(['  %s: median high-low = %.3g, increased %d/%d, ' ...
    'two-sided signrank p = %.3g, one-sided increase p = %.3g, ' ...
    'Spearman with %s change rho = %.3g, two-sided p = %.3g'], ...
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
    idx = ok & in_lock_group(neurons.lock_group, group_order{g});
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


function plot_first_latency_distributions(neurons, group_order, group_title, colors, split_name)
% Same neurons as the paired first-peak latency test: a first peak in both tails.
[lag_low, lag_high, ok] = paired_peak_values(neurons.low_peaks, neurons.high_peaks, 'first_lag');
if ~any(ok)
    return
end
figure('Name', ['First-peak latency, ', split_name], 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.06 0.08 0.88 0.78])
edges = 0.10:0.05:1.00;
panel_idx = cell(1, 1 + numel(group_order));
panel_idx{1} = ok;
for g = 1:numel(group_order)
    panel_idx{g + 1} = ok & in_lock_group(neurons.lock_group, group_order{g});
end
panel_title = [{'All neurons'}, group_title];
n_col = min(3, numel(panel_idx));
n_row = ceil(numel(panel_idx) / n_col);
for p = 1:numel(panel_idx)
    subplot(n_row, n_col, p)
    idx = panel_idx{p};
    low = lag_low(idx);
    high = lag_high(idx);
    if p == 1
        low_color = [0.75 0.75 0.75];
        high_color = [0.15 0.15 0.15];
    else
        low_color = 0.55 * colors(p - 1, :) + 0.45;
        high_color = colors(p - 1, :);
    end
    hold on
    if ~isempty(low)
        histogram(low, edges, 'Normalization', 'probability', ...
            'FaceColor', low_color, 'FaceAlpha', 0.45, 'EdgeColor', 'none')
        histogram(high, edges, 'Normalization', 'probability', ...
            'FaceColor', high_color, 'FaceAlpha', 0.45, 'EdgeColor', 'none')
        mark_latency_median(median(low, 'omitnan'), low_color)
        mark_latency_median(median(high, 'omitnan'), high_color)
        change = paired_change_stats(low, high, zeros(size(low)), zeros(size(high)));
        title(sprintf('%s, n = %d\nmedian high-low %.3g s, p = %.3g', ...
            panel_title{p}, change.n, change.median_delta, change.p_two))
    else
        title(sprintf('%s, n = 0', panel_title{p}))
    end
    xlim(edges([1 end]))
    xlabel('First-peak latency (s)')
    ylabel('Fraction of neurons')
    box off
    if p == 1 && ~isempty(low)
        legend({'Low 25%', 'High 25%'}, 'Location', 'northeast')
    end
end
end


function mark_latency_median(value, color)
if ~isfinite(value)
    return
end
plot([value value], ylim, '-', 'Color', color, 'LineWidth', 1.5, 'HandleVisibility', 'off')
end


function colors = lock_group_colors(group_order)
palette.trough = [0.20 0.45 0.80];
palette.peak = [0.85 0.25 0.20];
palette.nonentrained = [0.45 0.45 0.45];
palette.delta = [0.10 0.50 0.35];
if nargin < 1 || isempty(group_order)
    group_order = {'trough', 'peak', 'nonentrained'};
end
colors = zeros(numel(group_order), 3);
for g = 1:numel(group_order)
    colors(g, :) = palette.(group_order{g});
end
end


function tf = in_lock_group(lock_group, name)
if strcmp(name, 'delta')
    tf = ismember(lock_group, {'trough', 'peak'});
else
    tf = strcmp(lock_group, name);
end
end


function legend_group_colors(group_order, group_title, colors)
hold on
handles = gobjects(1, numel(group_order));
for g = 1:numel(group_order)
    handles(g) = plot(nan, nan, 'o', 'Color', colors(g, :), ...
        'MarkerFaceColor', colors(g, :), 'MarkerSize', 6, 'LineStyle', 'none');
end
legend(handles, group_title, 'Location', 'eastoutside')
end


function plot_peak_examples(summary, centers, smooth_n, group_order, group_title, score_label, split_name)
if ~isfield(summary, 'all_neurons') || ~isfield(summary.all_neurons, 'low') ...
        || summary.all_neurons.n < 1
    return
end
neurons = summary.all_neurons;
disp_block(['Peak examples, ', split_name])
t1_rows = select_strong_examples(neurons, group_order, 't1');
both_rows = select_strong_examples(neurons, group_order, 'both');
if ~any(t1_rows) && ~any(both_rows)
    disp(['No peak examples for ', split_name])
    return
end
n_col = numel(group_order);
figure('Name', ['ACG peak examples, ', split_name], 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.04 0.05 0.92 0.88])
row_rows = {t1_rows, both_rows};
row_names = {'Strong t1, no t2', 'Strong t1 and t2'};
for r = 1:2
    rows = row_rows{r};
    for c = 1:n_col
        subplot(2, n_col, (r - 1) * n_col + c)
        if c > numel(rows) || rows(c) < 1
            axis off
            continue
        end
        draw_peak_example(neurons, rows(c), centers, smooth_n, ...
            group_order, group_title, score_label, row_names{r}, split_name, group_title{c})
    end
end
sgtitle({['Reference for the peak scatters, ', split_name], ...
    'Line is the Gaussian used to detect the peaks. Grey band is the 0-50 ms burst.'})
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
rows = zeros(numel(group_order), 1);
for g = 1:numel(group_order)
    idx = find(in_lock_group(neurons.lock_group, group_order{g}) & strength > 0);
    if isempty(idx)
        continue
    end
    rows(g) = pick_high_rate_example(idx, strength, rate);
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


function draw_peak_example(neurons, row, centers, smooth_n, group_order, group_title, score_label, kind_name, split_name, panel_group)
if nargin < 10
    panel_group = '';
end
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
if ~isempty(panel_group)
    group_name = panel_group;
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


function plot_rate_scaled_mean_acgs(locking_acg, power_acg, bin_w, centers, group_order, group_title, locking_labels, power_labels, delta_lag, skip_abs_lag, smooth_n, curve_colors, cluster_alpha, cluster_perm)
locking_hz = rate_scaled_acg_stats(locking_acg, bin_w, group_order, group_title, ...
    locking_labels, 'locking', centers, delta_lag, skip_abs_lag, smooth_n, ...
    cluster_alpha, cluster_perm);
power_hz = rate_scaled_acg_stats(power_acg, bin_w, group_order, group_title, ...
    power_labels, 'delta power', centers, delta_lag, skip_abs_lag, smooth_n, ...
    cluster_alpha, cluster_perm);
acg_ylim = shared_mean_acg_ylim(locking_hz, power_hz, group_order, ...
    centers, skip_abs_lag, smooth_n, 'rate ratio (Poisson = 1)');
y_label = 'Ratio to Poisson rate';
name_suffix = ', rate ratio';
plot_summary_acgs(locking_hz, group_order, group_title, locking_labels, ...
    'locking', centers, delta_lag, skip_abs_lag, smooth_n, curve_colors, ...
    acg_ylim, y_label, name_suffix);
plot_summary_acgs(power_hz, group_order, group_title, power_labels, ...
    'delta power', centers, delta_lag, skip_abs_lag, smooth_n, curve_colors, ...
    acg_ylim, y_label, name_suffix);
end


function scaled = rate_scaled_acg_stats(summary, bin_w, group_order, group_title, curve_labels, split_name, centers, delta_lag, skip_abs_lag, smooth_n, cluster_alpha, cluster_perm)
scaled = rate_scale_summary(summary, bin_w, group_order);
names = [group_order, {'all_neurons'}];
titles = [group_title, {'All neurons'}];
curve_names = {'low', 'high', 'all'};
disp_block(['Mean autocorrelogram, rate ratio (Poisson = 1), ', split_name])
for g = 1:numel(names)
    if ~isfield(scaled, names{g}) || ~isfield(scaled.(names{g}), 'low')
        continue
    end
    disp_rule()
    scaled.(names{g}) = record_delta_mass(scaled.(names{g}), curve_names, curve_labels, ...
        [titles{g}, ', ', split_name, ', rate ratio'], centers, delta_lag, split_name, ...
        skip_abs_lag, smooth_n, cluster_alpha, cluster_perm);
end
end


function plot_rate_scaled_peaks_for(summary, bin_w, centers, smooth_n, group_order, group_title, score_label, split_name)
scaled = rate_scale_summary(summary, bin_w, group_order);
scaled = add_peak_scores(scaled, centers, smooth_n, group_order);
plot_paired_peak_scores(scaled, group_order, group_title, score_label, split_name, 'ratio');
end


function summary = rate_scale_summary(summary, bin_w, group_order)
names = [group_order, {'all_neurons'}];
curve_names = {'low', 'high', 'all'};
for k = 1:numel(names)
    if ~isfield(summary, names{k}) || ~isfield(summary.(names{k}), 'low')
        continue
    end
    group = summary.(names{k});
    for c = 1:numel(curve_names)
        curve = curve_names{c};
        group.(curve) = scale_acg_by_rate(group.(curve), ...
            group.(['mass_', curve]), group.(['rate_', curve]), bin_w);
    end
    summary.(names{k}) = group;
end
end


function Y = scale_acg_by_rate(Y, mass, rate, bin_w)
% (pair counts / overlap) / (tail firing rate ^ 2 * bin width).
% Y is the sum-to-1 curve and mass restores pair counts / overlap.
% The rate is the firing rate of the runs in that curve. A Poisson
% autocorrelogram then sits at 1.
factor = mass(:) ./ ((rate(:) .^ 2) * bin_w);
factor(~isfinite(factor) | rate(:) <= 0 | mass(:) <= 0) = NaN;
Y = Y .* factor;
end


function rates = tail_rate_struct(rt, pooled)
rates.low = tail_firing_rate(rt, pooled.low_index);
rates.high = tail_firing_rate(rt, pooled.high_index);
rates.all = neuron_run_rate(rt);
end


function rate = tail_firing_rate(rt, idx)
rate = NaN;
if ~ismember('RunNumSpikes', rt.Properties.VariableNames) || isempty(idx)
    return
end
idx = idx(:);
if numel(idx) ~= height(rt) || ~any(idx)
    return
end
n_spikes = double(rt.RunNumSpikes(:));
duration = acg_run_lengths(rt);
ok = idx & isfinite(n_spikes) & n_spikes >= 0 & isfinite(duration) & duration > 0;
if ~any(ok)
    return
end
rate = sum(n_spikes(ok)) / sum(duration(ok));
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


function name = canonical_kept_area(label, areas_to_keep)
% Map a depth-table label onto areas_to_keep. SupCol and SC are the
% superior colliculus. Tokens are split on commas so a combined label can
% match, while VLPAG does not match LPAG.
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
disp(sprintf('  %s two-sided signrank, high vs low 25%% %s: p = %.3g', group_label, split_name, p_low))
disp(sprintf('  %s two-sided signrank, high vs all runs: p = %.3g', group_label, p_all))
end


function report_acg_clusters(result, figure_label, split_name)
disp(sprintf(['  %s cluster test, high vs low 25%% %s, positive lags: ' ...
    'uncorrected signrank p < %.3g, %d sign flips'], ...
    figure_label, split_name, result.alpha, result.n_perm))
if isempty(result.clusters)
    disp('  No contiguous uncorrected difference.')
    return
end
for k = 1:numel(result.clusters)
    cluster = result.clusters(k);
    if cluster.survives
        fate = 'survives';
    else
        fate = 'does not survive';
    end
    disp(sprintf('  Cluster %.3g-%.3g s, cluster p = %.3g, %s', ...
        cluster.x0, cluster.x1, cluster.p, fate))
end
disp('  Only clusters that survive are shaded, on the positive lags and on their negative mirrors.')
end


function shade_acg_clusters(result, y_lim)
if isempty(result.clusters)
    return
end
span = y_lim(2) - y_lim(1);
y0 = y_lim(2) - 0.1 * span;
y1 = y_lim(2);
for k = 1:numel(result.clusters)
    if ~result.clusters(k).survives
        continue
    end
    x0 = result.clusters(k).x0;
    x1 = result.clusters(k).x1;
    patch([x0 x1 x1 x0], [y0 y0 y1 y1], [0.25 0.25 0.25], ...
        'FaceAlpha', 0.45, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    patch([-x1 -x0 -x0 -x1], [y0 y0 y1 y1], [0.25 0.25 0.25], ...
        'FaceAlpha', 0.45, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
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


function y_lim = shared_mean_acg_ylim(locking_acg, power_acg, group_order, centers, skip_abs_lag, smooth_n, scale_label)
% One vertical range for every mean autocorrelogram. The range is the
% mean curve plus and minus its standard error, not the delta-band shade.
if nargin < 7 || isempty(scale_label)
    scale_label = 'sum-to-1 normalized autocorrelogram';
end
curve_names = {'low', 'high', 'all'};
names = [group_order, {'all_neurons'}];
summaries = {locking_acg, power_acg};
lo = inf;
hi = -inf;
for s = 1:numel(summaries)
    summary = summaries{s};
    for k = 1:numel(names)
        if ~isfield(summary, names{k}) || ~isfield(summary.(names{k}), 'low')
            continue
        end
        group = summary.(names{k});
        for c = 1:numel(curve_names)
            shown = smooth_acg_away_from_zero( ...
                group.(curve_names{c}), centers, skip_abs_lag, smooth_n);
            mu = mean(shown, 1, 'omitnan');
            n = sum(isfinite(shown), 1);
            sem = std(shown, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
            band_lo = mu - sem;
            band_hi = mu + sem;
            band_lo = band_lo(isfinite(band_lo));
            band_hi = band_hi(isfinite(band_hi));
            if ~isempty(band_lo)
                lo = min(lo, min(band_lo));
            end
            if ~isempty(band_hi)
                hi = max(hi, max(band_hi));
            end
        end
    end
end
y_lim = readable_acg_scale(lo, hi);
disp_block(['Shared autocorrelogram scale, ', scale_label])
disp(sprintf(['Mean autocorrelograms span %.4g to %.4g. ' ...
    'Shared y limit is [%.4g %.4g].'], lo, hi, y_lim(1), y_lim(2)))
end


function y_lim = readable_acg_scale(lo, hi)
if ~isfinite(hi) || hi <= 0
    y_lim = [0 0.01];
    return
end
top = nice_decimal(hi, 'ceil');
if top <= hi
    top = nice_decimal(hi * 1.05, 'ceil');
end
if ~isfinite(lo) || lo >= 0
    bottom = 0;
else
    bottom = -nice_decimal(-lo, 'ceil');
    if bottom >= lo
        bottom = -nice_decimal(-lo * 1.05, 'ceil');
    end
end
y_lim = [bottom, top];
end


function ticks = readable_acg_ticks(y_lim)
span = y_lim(2) - y_lim(1);
step = nice_decimal(span / 4, 'ceil');
if step <= 0
    ticks = y_lim;
    return
end
ticks = y_lim(1):step:y_lim(2);
if numel(ticks) < 3
    step = nice_decimal(span / 4, 'floor');
    if step <= 0
        step = span / 4;
    end
    ticks = y_lim(1):step:y_lim(2);
end
if abs(ticks(end) - y_lim(2)) > step * 1e-6
    ticks(end + 1) = y_lim(2); %#ok<AGROW>
end
end


function y = nice_decimal(x, direction)
% Next 1, 2, or 5 times a power of ten, so tick labels stay short.
if ~isfinite(x) || x <= 0
    y = 0;
    return
end
exp10 = floor(log10(x));
base = 10 ^ exp10;
mantissa = x / base;
grid = [1 2 5 10];
if strcmp(direction, 'floor')
    pick = find(grid <= mantissa + 1e-9, 1, 'last');
    if isempty(pick)
        pick = 1;
    end
else
    pick = find(grid + 1e-9 >= mantissa, 1, 'first');
    if isempty(pick)
        pick = numel(grid);
    end
end
y = grid(pick) * base;
end


function shade = shade_delta_lags(lag_win)
shade = gobjects(1, 2);
for k = 1:2
    sign = 2 * (k - 1) - 1;
    x1 = sign * lag_win(1);
    x2 = sign * lag_win(2);
    shade(k) = patch([x1 x2 x2 x1], [0 0 0 0], [0.95 0.85 0.40], ...
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


function means = all_run_means(rt, acg, run_length, ppc, power)
% Means over the runs that enter the all-runs autocorrelogram.
all_idx = isfinite(run_length(:)) & run_length(:) > 0 & any(isfinite(acg), 2);
means.fr = tail_firing_rate(rt, all_idx);
means.ppc = mean_on(ppc, all_idx);
means.power = mean_on(power, all_idx);
end


function m = mean_on(values, idx)
values = double(values(:));
idx = idx(:) & isfinite(values);
if ~any(idx)
    m = NaN;
    return
end
m = mean(values(idx));
end


function model = fit_delta_power_contributions(summary, group_order, group_title, do_plot)
if nargin < 4
    do_plot = true;
end
% Mean delta power of the runs in the all-runs autocorrelogram, predicted
% by their mean firing rate, their mean PPC, the 0-50 ms burst, and the
% two delta-peak prominences of that autocorrelogram. One row per neuron.
model = struct();
if ~isfield(summary, 'all_neurons') || ~isfield(summary.all_neurons, 'all_peaks') ...
        || ~isfield(summary.all_neurons, 'model_power')
    model.all = empty_power_fit('All neurons');
    return
end
neurons = summary.all_neurons;
[X, y, names, kept] = power_model_design(neurons);
model.all = fit_one_power_model(X, y, names, 'All neurons');
model.groups = struct();
for g = 1:numel(group_order)
    in_group = in_lock_group(neurons.lock_group(kept), group_order{g});
    model.groups.(group_order{g}) = fit_one_power_model( ...
        X(in_group, :), y(in_group), names, group_title{g});
end
if do_plot
    plot_power_contributions(model, group_order, group_title, names);
end
end


function [X, y, names, kept] = power_model_design(neurons)
names = {'Firing rate', 'PPC', 'Burst, 0-50 ms', 'First-peak prominence', 'Second-cycle prominence'};
peaks = neurons.all_peaks;
y_all = neurons.model_power(:);
X_all = [neurons.model_fr(:), neurons.model_ppc(:), peaks.burst(:), ...
    peaks.first_prominence(:), peaks.second_prominence(:)];
kept = isfinite(y_all) & all(isfinite(X_all), 2) ...
    & peaks.first_measurable(:) & peaks.second_measurable(:);
X = X_all(kept, :);
y = y_all(kept);
end


function fit = fit_one_power_model(X, y, names, label)
fit = empty_power_fit(label);
fit.n = size(X, 1);
n_pred = numel(names);
if fit.n < n_pred + 5
    disp_block(['Delta power model, ', label])
    disp(['Delta power model, ', label, ': n = ', num2str(fit.n), ...
        ' is too small for ', num2str(n_pred), ' predictors.'])
    return
end
sd = std(X, 0, 1);
if any(~isfinite(sd) | sd <= 0) || ~isfinite(std(y)) || std(y) <= 0
    disp_block(['Delta power model, ', label])
    disp(['Delta power model, ', label, ': a variable has no spread.'])
    return
end
fit.ok = true;
fit.names = names;
Xs = (X - mean(X, 1)) ./ sd;
ys = (y - mean(y)) / std(y, 0);
[fit.beta, fit.rsq, fit.F, fit.p_model, fit.t_beta, fit.p_beta] = standardized_ols(Xs, ys);
fit.unique_rsq = zeros(n_pred, 1);
for k = 1:n_pred
    keep = true(1, n_pred);
    keep(k) = false;
    [~, rsq_without] = standardized_ols(Xs(:, keep), ys);
    fit.unique_rsq(k) = fit.rsq - rsq_without;
end
[fit.weight, fit.weight_rsq] = johnson_relative_weights(Xs, ys);
if fit.rsq > 0
    fit.percent_of_rsq = 100 * fit.weight / fit.rsq;
else
    fit.percent_of_rsq = nan(n_pred, 1);
end
fit.predictor_corr = corr(X);
report_power_model(fit);
end


function fit = empty_power_fit(label)
fit.ok = false;
fit.label = label;
fit.n = 0;
fit.names = {};
fit.beta = [];
fit.rsq = NaN;
fit.F = NaN;
fit.p_model = NaN;
fit.t_beta = [];
fit.p_beta = [];
fit.unique_rsq = [];
fit.weight = [];
fit.weight_rsq = NaN;
fit.percent_of_rsq = [];
fit.predictor_corr = [];
end


function [beta, rsq, Fstat, p_model, t_beta, p_beta] = standardized_ols(X, y)
n = size(X, 1);
k = size(X, 2);
A = [ones(n, 1), X];
coef = pinv(A) * y;
beta = coef(2:end);
t_beta = nan(k, 1);
p_beta = nan(k, 1);
resid = y - A * coef;
ss_res = sum(resid .^ 2);
ss_tot = sum((y - mean(y)) .^ 2);
if ss_tot <= 0
    rsq = NaN;
    Fstat = NaN;
    p_model = NaN;
    return
end
rsq = 1 - ss_res / ss_tot;
dof1 = k;
dof2 = n - k - 1;
if dof2 < 1 || ss_res < 0
    Fstat = NaN;
    p_model = NaN;
    return
end
Fstat = ((ss_tot - ss_res) / dof1) / (ss_res / dof2);
p_model = betainc(dof2 / (dof2 + dof1 * Fstat), dof2 / 2, dof1 / 2);
mse = ss_res / dof2;
covb = mse * pinv(A' * A);
se = sqrt(max(real(diag(covb)), 0));
tstat = coef ./ se;
tstat(~isfinite(se) | se <= 0) = NaN;
pstat = t_two_tailed(tstat, dof2);
t_beta = tstat(2:end);
p_beta = pstat(2:end);
end


function p = t_two_tailed(tstat, dof)
p = nan(size(tstat));
ok = isfinite(tstat) & isfinite(dof) & dof > 0;
ratio = dof ./ (dof + tstat(ok) .^ 2);
p(ok) = betainc(ratio, dof / 2, 0.5);
end


function [weight, rsq] = johnson_relative_weights(X, y)
% Relative weights for correlated predictors. The symmetric factor of the
% predictor correlation matrix splits R^2, including variance the
% predictors share, so the weights sum to the model R^2.
Rxx = corr(X);
Rxy = corr(X, y);
[Q, D] = eig((Rxx + Rxx') / 2);
eigvals = real(diag(D));
[eigvals, order] = sort(eigvals, 'descend');
Q = Q(:, order);
eigvals(eigvals < 1e-10) = 1e-10;
Lam = Q * diag(sqrt(eigvals)) * Q';
beta = Lam \ Rxy;
weight = (Lam .^ 2) * (beta .^ 2);
weight = real(weight(:));
rsq = sum(weight);
end


function report_power_model(fit)
disp_block(['Delta power model, ', fit.label])
disp(['Delta power model, ', fit.label, ', n = ', num2str(fit.n)])
disp(sprintf('  Full R^2 = %.3g, F = %.3g, p = %.3g', fit.rsq, fit.F, fit.p_model))
disp(['  Predictor correlations (', strjoin(fit.names, ', '), ')'])
disp(fit.predictor_corr)
for k = 1:numel(fit.names)
    disp(sprintf(['  %s: standardized beta = %.3g, t = %.3g, p = %.3g, ' ...
        'unique dR^2 = %.3g, relative weight = %.3g (%.1f%% of R^2)'], ...
        fit.names{k}, fit.beta(k), fit.t_beta(k), fit.p_beta(k), ...
        fit.unique_rsq(k), fit.weight(k), fit.percent_of_rsq(k)))
end
end


function plot_power_contributions(model, group_order, group_title, names)
if ~model.all.ok
    return
end
labels = [{'All neurons'}, group_title];
keys = [{'all'}, group_order];
pct = nan(numel(keys), numel(names));
for i = 1:numel(keys)
    if i == 1
        fit = model.all;
    else
        fit = model.groups.(keys{i});
    end
    if ~fit.ok
        continue
    end
    pct(i, :) = fit.percent_of_rsq(:)';
end
figure('Name', 'Delta power contributions', 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.15 0.2 0.7 0.5])
bar(pct)
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 20)
ylabel('Share of explained delta-power variance (%)')
legend(names, 'Location', 'eastoutside')
title(sprintf('All-runs autocorrelogram, R^2 = %.2f, n = %d', model.all.rsq, model.all.n))
box off
end
