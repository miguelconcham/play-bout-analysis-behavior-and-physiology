%% Analyze_spike_train_parsing
% Load firing-run parsing structs one file at a time (pooled .mat first,
% then per-animal files). Keep delta-locked neurons from the all_neurons
% table, split each cell's runs into low vs high PPC, edge-correct run
% autocorrelograms (counts / (T-|lag|); NaN if lag >= run length), sum,
% normalize, and compare. Existing GENERATE files without
% config.acg_edge_corrected are corrected here.

%% Paths and parameters
data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
run('\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
phase_folder   = fullfile(data_root, 'Analysis results', 'phase locking data');

alpha         = 0.01;    % delta lock: EntireSession PPC p-value
ppc_low_max   = 0.05;    % runs with PPC below this (try 0.1)
ppc_high_min  = 0.20;    % runs with PPC above this
min_runs      = 1;       % minimum runs in each PPC group to keep the cell
lag_win       = [0.15 0.40];  % lags (s) used as a delta-side-peak summary
smooth_half   = 0.04;    % s, ACG smoothing away from zero (as in parsing plots)
skip_abs_lag  = 0.02;    % s, do not smooth near zero lag

%% Load delta-locked labels (same table as Fig 3)
load(fullfile(phase_folder, 'delta_all_neurons_v2.mat'), 'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) = {'isRt'};
delta_lock = all_neurons.EntireSession.PPCPval < alpha & ~isnan(all_neurons.EntireSession.PPC);
id_all = cluster_id_column(all_neurons);
disp(['Delta-locked neurons in all_neurons: ', num2str(sum(delta_lock)), ' / ', num2str(height(all_neurons))])

%% Load + merge parsing files one at a time
% Pooled file is ~12 GB uncompressed; each later session ~3.6–3.9 GB.
% Do not keep every parsing_struct in the workspace.
parsing_files = list_parsing_mat_files(parsing_folder);
if isempty(parsing_files)
    error('Analyze_spike_train_parsing:empty', 'No parsing structs found.');
end

n_est = sum(delta_lock);
loaded_codes = {};
n_sessions = 0;
edges = [];
n_bins = [];
row = 0;
row_u = 0;
n_delta_in_parsing = 0;
n_unlock_in_parsing = 0;
n_missing_in_all_neurons = 0;
acg_low_raw = [];
acg_high_raw = [];
acg_unlock_raw = [];
acg_unlock_low_raw = [];
acg_unlock_high_raw = [];
n_runs_low = [];
n_runs_high = [];
n_runs_unlock = [];
n_runs_unlock_low = [];
n_runs_unlock_high = [];
neuron_id = [];
neuron_sess = {};
neuron_area = {};

for f = 1:numel(parsing_files)
    disp(['Loading ', parsing_files(f).name, ' ...'])
    S = load(parsing_files(f).path, 'parsing_struct');
    if ~isfield(S, 'parsing_struct')
        warning('No parsing_struct in %s', parsing_files(f).name);
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
            disp(['  already loaded, skip ', animal_code, ' (', parsing_files(f).name, ')'])
            clear ps
            continue
        end
        if isempty(edges)
            edges = ps.config.autocorr_bin_edges;
            bin_w = mean(diff(edges));
            centers = edges(2:end) - bin_w / 2;
            n_bins = numel(centers);
            smooth_n = max(1, round(smooth_half / bin_w));
            n_unl_est = max(sum(~delta_lock), 1);
            acg_low_raw  = nan(n_est, n_bins);
            acg_high_raw = nan(n_est, n_bins);
            acg_unlock_raw = nan(n_unl_est, n_bins);
            acg_unlock_low_raw = nan(n_unl_est, n_bins);
            acg_unlock_high_raw = nan(n_unl_est, n_bins);
            n_runs_low   = nan(n_est, 1);
            n_runs_high  = nan(n_est, 1);
            n_runs_unlock = nan(n_unl_est, 1);
            n_runs_unlock_low = nan(n_unl_est, 1);
            n_runs_unlock_high = nan(n_unl_est, 1);
            neuron_id    = nan(n_est, 1);
            neuron_sess  = cell(n_est, 1);
            neuron_area  = cell(n_est, 1);
        end

        ci = ps.cluster_info;
        id_col = cluster_id_column(ci);
        run_cell = ps.run_events_table_per_cell;
        ppc_col = ppc_column_index(ps);
        if isempty(ppc_col)
            warning('No PPC column in %s; skip session.', animal_code);
            clear ps
            continue
        end

        for nn = 1:height(ci)
            this_id = double(id_col(nn));
            hit = session_matches(all_neurons.session, animal_code) & id_all == this_id;
            if ~any(hit)
                n_missing_in_all_neurons = n_missing_in_all_neurons + 1;
                continue
            end
            is_locked = any(delta_lock(hit));
            if is_locked
                n_delta_in_parsing = n_delta_in_parsing + 1;
            else
                n_unlock_in_parsing = n_unlock_in_parsing + 1;
            end

            rt = run_cell{nn};
            if isempty(rt) || ~ismember('autoccorrelograms', rt.Properties.VariableNames)
                continue
            end
            pe = rt.phase_entrainment;
            if size(pe, 2) < ppc_col
                continue
            end
            ppc = pe(:, ppc_col);
            acg = rt.autoccorrelograms;
            if size(acg, 2) ~= n_bins
                warning('ACG width mismatch in %s id %g; skip cell.', animal_code, this_id);
                continue
            end
            already_corrected = isfield(ps, 'config') && isfield(ps.config, 'acg_edge_corrected') ...
                && ps.config.acg_edge_corrected;
            T = acg_run_lengths(rt);
            if numel(T) ~= size(acg, 1)
                warning('Run length mismatch in %s id %g; ACG not overlap-pooled.', animal_code, this_id);
                T = nan(size(acg, 1), 1);
            end

            low_idx  = ppc < ppc_low_max & ~isnan(ppc);
            high_idx = ppc > ppc_high_min & ~isnan(ppc);
            if is_locked
                if sum(low_idx) < min_runs && sum(high_idx) < min_runs
                    continue
                end
                row = row + 1;
                if row > size(acg_low_raw, 1)
                    extra = 200;
                    acg_low_raw  = [acg_low_raw; nan(extra, n_bins)]; %#ok<AGROW>
                    acg_high_raw = [acg_high_raw; nan(extra, n_bins)]; %#ok<AGROW>
                    n_runs_low   = [n_runs_low; nan(extra, 1)]; %#ok<AGROW>
                    n_runs_high  = [n_runs_high; nan(extra, 1)]; %#ok<AGROW>
                    neuron_id    = [neuron_id; nan(extra, 1)]; %#ok<AGROW>
                    neuron_sess  = [neuron_sess; cell(extra, 1)]; %#ok<AGROW>
                    neuron_area  = [neuron_area; cell(extra, 1)]; %#ok<AGROW>
                end
                neuron_id(row)   = this_id;
                neuron_sess{row} = animal_code;
                neuron_area{row} = cluster_area_char(ci, nn);
                n_runs_low(row)  = sum(low_idx);
                n_runs_high(row) = sum(high_idx);
                if sum(low_idx) >= min_runs
                    acg_low_raw(row, :) = pool_acg_runs(acg, edges, T, low_idx, already_corrected);
                end
                if sum(high_idx) >= min_runs
                    acg_high_raw(row, :) = pool_acg_runs(acg, edges, T, high_idx, already_corrected);
                end
            else
                use_idx = any(isfinite(acg), 2);
                if sum(use_idx) < min_runs
                    continue
                end
                row_u = row_u + 1;
                if row_u > size(acg_unlock_raw, 1)
                    extra = 200;
                    acg_unlock_raw = [acg_unlock_raw; nan(extra, n_bins)]; %#ok<AGROW>
                    acg_unlock_low_raw = [acg_unlock_low_raw; nan(extra, n_bins)]; %#ok<AGROW>
                    acg_unlock_high_raw = [acg_unlock_high_raw; nan(extra, n_bins)]; %#ok<AGROW>
                    n_runs_unlock  = [n_runs_unlock; nan(extra, 1)]; %#ok<AGROW>
                    n_runs_unlock_low = [n_runs_unlock_low; nan(extra, 1)]; %#ok<AGROW>
                    n_runs_unlock_high = [n_runs_unlock_high; nan(extra, 1)]; %#ok<AGROW>
                end
                n_runs_unlock(row_u) = sum(use_idx);
                n_runs_unlock_low(row_u) = sum(low_idx);
                n_runs_unlock_high(row_u) = sum(high_idx);
                acg_unlock_raw(row_u, :) = pool_acg_runs(acg, edges, T, use_idx, already_corrected);
                if sum(low_idx) >= min_runs
                    acg_unlock_low_raw(row_u, :) = pool_acg_runs(acg, edges, T, low_idx, already_corrected);
                end
                if sum(high_idx) >= min_runs
                    acg_unlock_high_raw(row_u, :) = pool_acg_runs(acg, edges, T, high_idx, already_corrected);
                end
            end
        end

        loaded_codes{end + 1} = animal_code; %#ok<AGROW>
        n_sessions = n_sessions + 1;
        disp(['  processed ', animal_code, '  (', parsing_files(f).kind, ')'])
        clear ps ci run_cell
    end
    clear S
end

disp(['Parsing sessions processed: ', num2str(n_sessions)])
if isempty(edges)
    error('Analyze_spike_train_parsing:empty', 'No parsing structs found.');
end

keep = 1:row;
acg_low_raw  = acg_low_raw(keep, :);
acg_high_raw = acg_high_raw(keep, :);
n_runs_low   = n_runs_low(keep);
n_runs_high  = n_runs_high(keep);
neuron_id    = neuron_id(keep);
neuron_sess  = neuron_sess(keep);
neuron_area  = neuron_area(keep);
if isempty(acg_unlock_raw)
    acg_unlock_raw = nan(0, n_bins);
    acg_unlock_low_raw = nan(0, n_bins);
    acg_unlock_high_raw = nan(0, n_bins);
    n_runs_unlock = nan(0, 1);
    n_runs_unlock_low = nan(0, 1);
    n_runs_unlock_high = nan(0, 1);
else
    acg_unlock_raw = acg_unlock_raw(1:row_u, :);
    acg_unlock_low_raw = acg_unlock_low_raw(1:row_u, :);
    acg_unlock_high_raw = acg_unlock_high_raw(1:row_u, :);
    n_runs_unlock  = n_runs_unlock(1:row_u);
    n_runs_unlock_low = n_runs_unlock_low(1:row_u);
    n_runs_unlock_high = n_runs_unlock_high(1:row_u);
end

if row < 1
    warning('Analyze_spike_train_parsing:none', ...
        'No delta-locked cells with low/high PPC runs. Check ppc_low_max / ppc_high_min and that parsing files exist.');
    return
end

acg_low_norm  = normalize_rows(acg_low_raw);
acg_high_norm = normalize_rows(acg_high_raw);
acg_unlock_norm = normalize_rows(acg_unlock_raw);
acg_unlock_low_norm = normalize_rows(acg_unlock_low_raw);
acg_unlock_high_norm = normalize_rows(acg_unlock_high_raw);

disp(['Delta-locked cells found in parsing: ', num2str(n_delta_in_parsing)])
disp(['Cells with usable low and/or high PPC runs: ', num2str(row)])
disp(['  both groups: ', num2str(sum(n_runs_low >= min_runs & n_runs_high >= min_runs))])
disp(['Unlocked cells with usable runs: ', num2str(row_u), ' / ', num2str(n_unlock_in_parsing)])
disp(['  unlocked low-PPC runs:  ', num2str(sum(n_runs_unlock_low >= min_runs))])
disp(['  unlocked high-PPC runs: ', num2str(sum(n_runs_unlock_high >= min_runs))])
disp(['all_neurons misses (parsing cell not in table): ', num2str(n_missing_in_all_neurons)])

%% Plots: mean ACG low vs high PPC
paired = n_runs_low >= min_runs & n_runs_high >= min_runs ...
    & any(isfinite(acg_low_norm), 2) & any(isfinite(acg_high_norm), 2);

low_sm  = smooth_acg_away_from_zero(acg_low_norm, centers, skip_abs_lag, smooth_n);
high_sm = smooth_acg_away_from_zero(acg_high_norm, centers, skip_abs_lag, smooth_n);
unlock_low_sm = smooth_acg_away_from_zero(acg_unlock_low_norm, centers, skip_abs_lag, smooth_n);
unlock_high_sm = smooth_acg_away_from_zero(acg_unlock_high_norm, centers, skip_abs_lag, smooth_n);
unlock_low_ok = n_runs_unlock_low >= min_runs;
unlock_high_ok = n_runs_unlock_high >= min_runs;
col_unlock = [0.55 0.55 0.55];
col_unlock_high = [0.25 0.25 0.25];

if ~any(paired)
    warning('Analyze_spike_train_parsing:paired', ...
        'No cells had both PPC < %g and PPC > %g runs. Relax ppc_low_max / ppc_high_min or min_runs.', ...
        ppc_low_max, ppc_high_min);
end

figure('Name', 'ACG low vs high PPC', 'units', 'normalized', 'outerposition', [0.1 0.1 0.7 0.8])
subplot(2, 2, 1)
hold on
leg = {};
if any(unlock_low_ok)
    plot_mean_sem(centers, unlock_low_sm(unlock_low_ok, :), col_unlock)
    leg{end + 1} = sprintf('Unlocked, PPC < %g (n=%d)', ppc_low_max, sum(unlock_low_ok));
end
if any(unlock_high_ok)
    plot_mean_sem(centers, unlock_high_sm(unlock_high_ok, :), col_unlock_high, '--')
    leg{end + 1} = sprintf('Unlocked, PPC > %g (n=%d)', ppc_high_min, sum(unlock_high_ok));
end
plot_mean_sem(centers, low_sm(paired, :), [0.2 0.2 0.8])
plot_mean_sem(centers, high_sm(paired, :), [0.85 0.2 0.2])
leg{end + 1} = sprintf('Locked, PPC < %g (n=%d)', ppc_low_max, sum(paired));
leg{end + 1} = sprintf('Locked, PPC > %g', ppc_high_min);
xline(0, ':k')
xlabel('Lag (s)')
ylabel('Normalized ACG')
legend(leg, 'Location', 'best')
title('Delta-locked vs unlocked (split by run PPC), mean ± SEM')
xlim([centers(1) centers(end)])

subplot(2, 2, 2)
hold on
if any(unlock_low_ok)
    plot_mean_sem(centers, zscore_rows(unlock_low_sm(unlock_low_ok, :)), col_unlock)
end
if any(unlock_high_ok)
    plot_mean_sem(centers, zscore_rows(unlock_high_sm(unlock_high_ok, :)), col_unlock_high, '--')
end
if any(paired)
    low_z  = zscore_rows(low_sm(paired, :));
    high_z = zscore_rows(high_sm(paired, :));
    plot_mean_sem(centers, low_z, [0.2 0.2 0.8])
    plot_mean_sem(centers, high_z, [0.85 0.2 0.2])
end
xline(0, ':k')
xlabel('Lag (s)')
ylabel('z-scored ACG')
title('Shape (z-scored per cell)')
xlim([centers(1) centers(end)])

subplot(2, 2, 3)
if any(paired)
    imagesc(centers, 1:sum(paired), low_sm(paired, :))
    axis xy
    colorbar
end
clim([0 0.02])
xlabel('Lag (s)')
ylabel('Neuron')
title('Low PPC')

subplot(2, 2, 4)
if any(paired)
    imagesc(centers, 1:sum(paired), high_sm(paired, :))
    axis xy
    colorbar
end
clim([0 0.02])
xlabel('Lag (s)')
ylabel('Neuron')
title('High PPC')

%% Side-peak summary (paired)
lag_idx = abs(centers) >= lag_win(1) & abs(centers) <= lag_win(2);
mod_low  = mean(acg_low_norm(paired, lag_idx), 2, 'omitnan');
mod_high = mean(acg_high_norm(paired, lag_idx), 2, 'omitnan');
if sum(paired) >= 2
    p_sr = signrank(mod_high, mod_low);
else
    p_sr = NaN;
end
disp(['Paired signrank, ACG mass in ', num2str(lag_win), ' s: p = ', num2str(p_sr)])
disp(['  mean low = ', num2str(mean(mod_low, 'omitnan')), '   mean high = ', num2str(mean(mod_high, 'omitnan'))])
mod_unlock_low = mean(acg_unlock_low_norm(unlock_low_ok, lag_idx), 2, 'omitnan');
mod_unlock_high = mean(acg_unlock_high_norm(unlock_high_ok, lag_idx), 2, 'omitnan');
disp(['  mean unlocked low-PPC  = ', num2str(mean(mod_unlock_low, 'omitnan')), '   n = ', num2str(sum(unlock_low_ok))])
disp(['  mean unlocked high-PPC = ', num2str(mean(mod_unlock_high, 'omitnan')), '   n = ', num2str(sum(unlock_high_ok))])

figure('Name', 'ACG lag-window paired')
hold on
if any(unlock_low_ok)
    yline(mean(mod_unlock_low, 'omitnan'), 'Color', col_unlock, 'LineStyle', '-', 'LineWidth', 1.5)
end
if any(unlock_high_ok)
    yline(mean(mod_unlock_high, 'omitnan'), 'Color', col_unlock_high, 'LineStyle', '--', 'LineWidth', 1.5)
end
if ~any(paired)
    warning('Analyze_spike_train_parsing:paired', ...
        'No cells had both PPC < %g and PPC > %g runs.', ppc_low_max, ppc_high_min);
else
    plot([1 2], [mod_low mod_high]', '-', 'Color', [0.7 0.7 0.7])
    plot(1, mean(mod_low, 'omitnan'), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 8)
    plot(2, mean(mod_high, 'omitnan'), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 8)
end
xlim([0.5 2.5])
set(gca, 'XTick', [1 2], 'XTickLabel', {sprintf('PPC < %g', ppc_low_max), sprintf('PPC > %g', ppc_high_min)})
ylabel(sprintf('Mean normalized ACG in [%g %g] s', lag_win(1), lag_win(2)))
title({['Paired signrank p = ', num2str(p_sr, '%.3g'), '   n = ', num2str(sum(paired))], ...
       'Gray: unlocked low-PPC (solid) / high-PPC (dashed)'})

%% Workspace summary table
acg_summary = table(neuron_sess, neuron_id, neuron_area, n_runs_low, n_runs_high, ...
    'VariableNames', {'session', 'cluster_id', 'area', 'n_runs_low_ppc', 'n_runs_high_ppc'});


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


function mu = pool_acg_runs(acg, bin_edges, run_length, idx, already_corrected)
% Sum raw pair counts / sum overlap (T-|lag|). Short runs do not inflate
% long lags the way averaging per-run (counts/overlap) does.
    n_bins = size(acg, 2);
    mu = nan(1, n_bins);
    if ~any(idx)
        return
    end
    bin_edges = bin_edges(:)';
    bin_w = mean(diff(bin_edges));
    centers = bin_edges(1:end-1) + 0.5 * bin_w;
    overlap = run_length(:) - abs(centers);
    overlap(overlap < bin_w | ~isfinite(overlap)) = NaN;
    C = acg(idx, :);
    W = overlap(idx, :);
    if already_corrected
        C = C .* W;
    end
    den = sum(W, 1, 'omitnan');
    num = sum(C, 1, 'omitnan');
    ok = den > 0;
    mu(ok) = num(ok) ./ den(ok);
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
        error('cluster_info has neither cluster_id nor id');
    end
end


function col = ppc_column_index(ps)
    col = [];
    if isfield(ps, 'config') && isfield(ps.config, 'phase_entrainment_fields')
        f = ps.config.phase_entrainment_fields;
        hit = find(strcmp(f, 'PPC'), 1);
        if ~isempty(hit)
            col = hit;
            return
        end
    end
    w = [];
    if isfield(ps, 'ALL_TABLES') && ~isempty(ps.ALL_TABLES) ...
            && ismember('phase_entrainment', ps.ALL_TABLES.Properties.VariableNames)
        w = size(ps.ALL_TABLES.phase_entrainment, 2);
    elseif isfield(ps, 'run_events_table_per_cell')
        for j = 1:numel(ps.run_events_table_per_cell)
            rt = ps.run_events_table_per_cell{j};
            if ~isempty(rt) && istable(rt) && ismember('phase_entrainment', rt.Properties.VariableNames)
                w = size(rt.phase_entrainment, 2);
                break
            end
        end
    end
    if ~isempty(w) && w >= 5
        col = 5;
    end
end


function tf = session_matches(session_col, animal_code)
    if iscell(session_col)
        tf = strcmp(session_col, animal_code);
        empty = cellfun(@isempty, session_col);
        if any(empty)
            tf(empty) = false;
        end
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


function Y = normalize_rows(X)
    s = sum(X, 2, 'omitnan');
    Y = X ./ s;
    Y(~isfinite(Y)) = NaN;
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


function Z = zscore_rows(X)
    if isempty(X)
        Z = X;
        return
    end
    mu = mean(X, 2, 'omitnan');
    sd = std(X, 0, 2, 'omitnan');
    sd(sd == 0 | isnan(sd)) = 1;
    Z = (X - mu) ./ sd;
end


function plot_mean_sem(t, M, col, ls)
    if nargin < 4 || isempty(ls)
        ls = '-';
    end
    if isempty(M)
        return
    end
    mu = mean(M, 1, 'omitnan');
    n = sum(isfinite(M), 1);
    sem = std(M, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
    fill([t fliplr(t)], [mu - sem, fliplr(mu + sem)], col, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off')
    hold on
    plot(t, mu, 'Color', col, 'LineStyle', ls, 'LineWidth', 2)
end
