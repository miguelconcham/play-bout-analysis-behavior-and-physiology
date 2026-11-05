%% Analyze_locking_period_acg
% Mean autocorrelogram of spikes in low-locking periods, high-locking
% periods, and every firing run. Reads the delta-lock table and the firing
% run tables written by Estimate_spike_train_parsing. Does not recompute
% phase locking.
%
% Each neuron is labeled from EntireSession in delta_all_neurons_v2.mat:
% trough-locked, peak-locked (PPC p < alpha, split at +/- pi/2), or
% non-entrained (PPC p >= 0.1). For each neuron the saved runs are split
% at that neuron's own PPC: bottom 25% and top 25%. A neuron is
% plotted only when at least four of its runs have a PPC, so the two
% tails are separated. The third curve pools
% every firing run of that neuron. The estimate stores autocorrelograms
% inside firing runs, so spikes that fall in the gaps between runs are not
% in these curves.
%
% Run autocorrelograms are overlap-weighted (pair counts / (T-|lag|)), the
% lag-0 bins are blanked, and each neuron is scaled to sum to 1 so a cell
% with more spikes does not dominate the mean. One figure per lock group.
% Shading marks lags of 0.2-0.5 s (2-5 Hz). A 1 Hz peak sits outside the
% +/-0.5 s window saved with the runs.

%% Paths and parameters
data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
run('\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
phase_folder   = fullfile(data_root, 'Analysis results', 'phase locking data');

alpha            = 0.01;
non_entrained_p  = 0.1;
low_percent      = 25;
high_percent     = 75;
delta_lag        = [0.2 0.5];   % s, 5 Hz to 2 Hz
smooth_half      = 0.04;        % s, display smoothing away from zero
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
n_in_parsing = zeros(1, numel(group_order));
n_kept = zeros(1, numel(group_order));
n_unsplit = zeros(1, numel(group_order));
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
        if isempty(edges)
            edges = ps.config.autocorr_bin_edges(:)';
            bin_w = mean(diff(edges));
            centers = edges(1:end - 1) + bin_w / 2;
            n_bins = numel(centers);
            acg_raw = empty_group_store(group_order, n_bins);
        elseif numel(ps.config.autocorr_bin_edges) ~= numel(edges) ...
                || max(abs(ps.config.autocorr_bin_edges(:)' - edges)) > 1e-9
            warning('Analyze_locking_period_acg:edges', ...
                'ACG bins differ in %s; session skipped.', animal_code);
            clear ps
            continue
        end

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
                continue
            end
            if size(rt.phase_entrainment, 2) < ppc_col
                n_unsplit(g) = n_unsplit(g) + 1;
                continue
            end
            acg = rt.autoccorrelograms;
            if size(acg, 2) ~= n_bins || size(acg, 1) ~= height(rt)
                warning('Analyze_locking_period_acg:width', ...
                    'ACG size mismatch in %s id %g; cell skipped.', animal_code, this_id);
                n_unsplit(g) = n_unsplit(g) + 1;
                continue
            end
            pooled = pool_locking_period_acgs(acg, edges, acg_run_lengths(rt), ...
                rt.phase_entrainment(:, ppc_col), already_corrected, low_percent, high_percent);
            if pooled.n_low < 1 || pooled.n_high < 1 ...
                    || ~any(isfinite(pooled.low)) || ~any(isfinite(pooled.high)) ...
                    || ~any(isfinite(pooled.all))
                n_unsplit(g) = n_unsplit(g) + 1;
                continue
            end

            row = acg_raw.(group_order{g}).n + 1;
            acg_raw.(group_order{g}).n = row;
            acg_raw.(group_order{g}).low(row, :) = pooled.low;
            acg_raw.(group_order{g}).high(row, :) = pooled.high;
            acg_raw.(group_order{g}).all(row, :) = pooled.all;
            acg_raw.(group_order{g}).n_low(row, 1) = pooled.n_low;
            acg_raw.(group_order{g}).n_high(row, 1) = pooled.n_high;
            acg_raw.(group_order{g}).n_all(row, 1) = pooled.n_all;
            acg_raw.(group_order{g}).session{row, 1} = animal_code;
            acg_raw.(group_order{g}).cluster_id(row, 1) = this_id;
            acg_raw.(group_order{g}).area{row, 1} = cluster_area_char(ci, nn);
            n_kept(g) = n_kept(g) + 1;
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
curve_names = {'low', 'high', 'all'};
curve_labels = { ...
    sprintf('Low %g%% locking', low_percent), ...
    sprintf('High %g%% locking', 100 - high_percent), ...
    'All firing runs'};
curve_colors = [0.20 0.45 0.75; 0.85 0.25 0.20; 0.15 0.15 0.15];
locking_acg = struct();

for g = 1:numel(group_order)
    name = group_order{g};
    n = acg_raw.(name).n;
    disp([group_title{g}, ' in parsing: ', num2str(n_in_parsing(g)), ...
        '   plotted: ', num2str(n), ...
        '   skipped (no separated 25/75 split or no spikes): ', num2str(n_unsplit(g))])
    if n < 1
        warning('Analyze_locking_period_acg:group', ...
            'No %s neurons with both locking tails.', group_title{g});
        continue
    end
    keep = 1:n;
    for c = 1:numel(curve_names)
        raw = acg_raw.(name).(curve_names{c})(keep, :);
        raw = blank_acg_zero_lag(raw, centers, acg_refractory_s);
        locking_acg.(name).(curve_names{c}) = normalize_acg_rows(raw);
    end
    locking_acg.(name).session = acg_raw.(name).session(keep);
    locking_acg.(name).cluster_id = acg_raw.(name).cluster_id(keep);
    locking_acg.(name).area = acg_raw.(name).area(keep);
    locking_acg.(name).n_low = acg_raw.(name).n_low(keep);
    locking_acg.(name).n_high = acg_raw.(name).n_high(keep);
    locking_acg.(name).n_all = acg_raw.(name).n_all(keep);

    figure('Name', group_title{g}, 'Color', 'w', ...
        'Units', 'normalized', 'OuterPosition', [0.15 0.15 0.55 0.55])
    hold on
    shade = shade_delta_lags(delta_lag);
    plotted = gobjects(1, numel(curve_names));
    for c = 1:numel(curve_names)
        shown = smooth_acg_away_from_zero( ...
            locking_acg.(name).(curve_names{c}), centers, skip_abs_lag, smooth_n);
        plotted(c) = plot_mean_sem(centers, shown, curve_colors(c, :));
        mass = delta_side_mass(locking_acg.(name).(curve_names{c}), centers, delta_lag);
        disp(sprintf('  %s  %s   mean ACG in [%g %g] s = %.4g   n = %d', ...
            group_title{g}, curve_labels{c}, delta_lag(1), delta_lag(2), ...
            mean(mass, 'omitnan'), sum(isfinite(mass))))
        locking_acg.(name).([curve_names{c}, '_delta_mass']) = mass;
    end
    report_locking_contrast(locking_acg.(name), group_title{g});
    xline(0, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off')
    xlim([centers(1) centers(end)])
    xlabel('Lag (s)')
    ylabel('Mean normalized autocorrelogram')
    title({group_title{g}, ...
        sprintf('n = %d    shaded: 2-5 Hz', n)})
    legend(plotted, curve_labels, 'Location', 'best')
    place_delta_shade(shade)
    box off
end

locking_acg.centers = centers;
locking_acg.bin_edges = edges;
locking_acg.parameters = struct( ...
    'alpha', alpha, ...
    'non_entrained_p', non_entrained_p, ...
    'low_percent', low_percent, ...
    'high_percent', high_percent, ...
    'delta_lag', delta_lag);

summary_path = fullfile(parsing_folder, 'locking_period_acg_summary.mat');
save(summary_path, 'locking_acg');
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
    block.session = {};
    block.cluster_id = zeros(0, 1);
    block.area = {};
    store.(group_order{g}) = block;
end
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
        y = smooth(X(j, pos), smooth_n);
        Y(j, pos) = y(:)';
    end
    if any(neg)
        y = smooth(X(j, neg), smooth_n);
        Y(j, neg) = y(:)';
    end
end
end


function mass = delta_side_mass(X, centers, lag_win)
lag_idx = abs(centers) >= lag_win(1) & abs(centers) <= lag_win(2);
mass = mean(X(:, lag_idx), 2, 'omitnan');
end


function report_locking_contrast(group_acg, group_label)
p_low = paired_signrank(group_acg.high_delta_mass, group_acg.low_delta_mass);
p_all = paired_signrank(group_acg.high_delta_mass, group_acg.all_delta_mass);
disp(sprintf('  %s signrank, high vs low 25%%: p = %.3g', group_label, p_low))
disp(sprintf('  %s signrank, high vs all runs: p = %.3g', group_label, p_all))
end


function p = paired_signrank(a, b)
ok = isfinite(a) & isfinite(b);
a = a(ok);
b = b(ok);
if numel(a) < 2 || all(a == b)
    p = NaN;
    return
end
p = signrank(a, b);
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
