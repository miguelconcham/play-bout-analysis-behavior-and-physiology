%% Analyze_spike_train_parsing
% Load firing-run parsing structs (first sessions from the pooled .mat,
% later sessions from per-animal files), keep delta-locked neurons from
% the all_neurons table, split each cell's runs into low vs high PPC, sum
% run autocorrelograms (counts), normalize, and compare.

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

%% Load parsing structs: pooled file first, then per-session files
[parsing_all, loaded_from] = load_parsing_structs(parsing_folder);
disp(['Parsing sessions loaded: ', num2str(numel(parsing_all))])
for k = 1:numel(parsing_all)
    disp(['  ', parsing_all{k}.animal_code, '  (', loaded_from{k}, ')'])
end
if isempty(parsing_all)
    error('Analyze_spike_train_parsing:empty', 'No parsing structs found.');
end

edges = parsing_all{1}.config.autocorr_bin_edges;
bin_w = mean(diff(edges));
centers = edges(2:end) - bin_w / 2;
n_bins = numel(centers);
smooth_n = max(1, round(smooth_half / bin_w));

%% Per delta-locked neuron: sum ACG of low-PPC vs high-PPC runs
n_est = sum(delta_lock);
acg_low_raw  = nan(n_est, n_bins);
acg_high_raw = nan(n_est, n_bins);
n_runs_low   = nan(n_est, 1);
n_runs_high  = nan(n_est, 1);
neuron_id    = nan(n_est, 1);
neuron_sess  = cell(n_est, 1);
neuron_area  = cell(n_est, 1);
row = 0;
n_delta_in_parsing = 0;
n_missing_in_all_neurons = 0;

for ps = 1:numel(parsing_all)
    animal_code = char(string(parsing_all{ps}.animal_code));
    ci = parsing_all{ps}.cluster_info;
    id_col = cluster_id_column(ci);
    run_cell = parsing_all{ps}.run_events_table_per_cell;
    ppc_col = ppc_column_index(parsing_all{ps});
    if isempty(ppc_col)
        warning('No PPC column in %s; skip session.', animal_code);
        continue
    end

    for nn = 1:height(ci)
        this_id = double(id_col(nn));
        hit = session_matches(all_neurons.session, animal_code) & id_all == this_id;
        if ~any(hit)
            n_missing_in_all_neurons = n_missing_in_all_neurons + 1;
            continue
        end
        if ~any(delta_lock(hit))
            continue
        end
        n_delta_in_parsing = n_delta_in_parsing + 1;

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

        low_idx  = ppc < ppc_low_max & ~isnan(ppc);
        high_idx = ppc > ppc_high_min & ~isnan(ppc);
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
        end;
        neuron_id(row)   = this_id;
        neuron_sess{row} = animal_code;
        neuron_area{row} = cluster_area_char(ci, nn);
        n_runs_low(row)  = sum(low_idx);
        n_runs_high(row) = sum(high_idx);
        if sum(low_idx) >= min_runs
            acg_low_raw(row, :) = sum(acg(low_idx, :), 1, 'omitnan');
        end
        if sum(high_idx) >= min_runs
            acg_high_raw(row, :) = sum(acg(high_idx, :), 1, 'omitnan');
        end
    end
end

keep = 1:row;
acg_low_raw  = acg_low_raw(keep, :);
acg_high_raw = acg_high_raw(keep, :);
n_runs_low   = n_runs_low(keep);
n_runs_high  = n_runs_high(keep);
neuron_id    = neuron_id(keep);
neuron_sess  = neuron_sess(keep);
neuron_area  = neuron_area(keep);

if row < 1
    warning('Analyze_spike_train_parsing:none', ...
        'No delta-locked cells with low/high PPC runs. Check ppc_low_max / ppc_high_min and that parsing files exist.');
    return
end

acg_low_norm  = normalize_rows(acg_low_raw);
acg_high_norm = normalize_rows(acg_high_raw);

disp(['Delta-locked cells found in parsing: ', num2str(n_delta_in_parsing)])
disp(['Cells with usable low and/or high PPC runs: ', num2str(row)])
disp(['  both groups: ', num2str(sum(n_runs_low >= min_runs & n_runs_high >= min_runs))])
disp(['all_neurons misses (parsing cell not in table): ', num2str(n_missing_in_all_neurons)])

%% Plots: mean ACG low vs high PPC
paired = n_runs_low >= min_runs & n_runs_high >= min_runs ...
    & all(isfinite(acg_low_norm), 2) & all(isfinite(acg_high_norm), 2);

low_sm  = smooth_acg_away_from_zero(acg_low_norm, centers, skip_abs_lag, smooth_n);
high_sm = smooth_acg_away_from_zero(acg_high_norm, centers, skip_abs_lag, smooth_n);

if ~any(paired)
    warning('Analyze_spike_train_parsing:paired', ...
        'No cells had both PPC < %g and PPC > %g runs. Relax ppc_low_max / ppc_high_min or min_runs.', ...
        ppc_low_max, ppc_high_min);
end

figure('Name', 'ACG low vs high PPC', 'units', 'normalized', 'outerposition', [0.1 0.1 0.7 0.8])
subplot(2, 2, 1)
hold on
plot_mean_sem(centers, low_sm(paired, :), [0.2 0.2 0.8])
plot_mean_sem(centers, high_sm(paired, :), [0.85 0.2 0.2])
xline(0, ':k')
xlabel('Lag (s)')
ylabel('Normalized ACG')
legend({sprintf('PPC < %g (n=%d)', ppc_low_max, sum(paired)), ...
        sprintf('PPC > %g', ppc_high_min)}, 'Location', 'best')
title('Delta-locked cells, mean ± SEM')
xlim([centers(1) centers(end)])

subplot(2, 2, 2)
hold on
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
xlabel('Lag (s)')
ylabel('Neuron')
title('Low PPC')

subplot(2, 2, 4)
if any(paired)
    imagesc(centers, 1:sum(paired), high_sm(paired, :))
    axis xy
    colorbar
end
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

figure('Name', 'ACG lag-window paired')
hold on
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
title(['Paired signrank p = ', num2str(p_sr, '%.3g'), '   n = ', num2str(sum(paired))])

%% Workspace summary table
acg_summary = table(neuron_sess, neuron_id, neuron_area, n_runs_low, n_runs_high, ...
    'VariableNames', {'session', 'cluster_id', 'area', 'n_runs_low_ppc', 'n_runs_high_ppc'});


function [parsing_all, loaded_from] = load_parsing_structs(parsing_folder)
    parsing_all = {};
    loaded_from = {};
    loaded_codes = {};

    combined = fullfile(parsing_folder, 'spike_train_parsing_structure.mat');
    if exist(combined, 'file')
        S = load(combined, 'parsing_struct');
        for k = 1:numel(S.parsing_struct)
            parsing_all{end + 1} = S.parsing_struct(k); %#ok<AGROW>
            loaded_codes{end + 1} = char(string(S.parsing_struct(k).animal_code)); %#ok<AGROW>
            loaded_from{end + 1} = 'pooled'; %#ok<AGROW>
        end
        disp(['Loaded ', num2str(numel(S.parsing_struct)), ' session(s) from pooled file'])
        clear S
    end

    flist = dir(fullfile(parsing_folder, '*spike_train_parsing*.mat'));
    skip = {'spike_train_parsing_structure.mat', 'spike_train_parsing_animal_names.mat'};
    for f = 1:numel(flist)
        if ismember(flist(f).name, skip)
            continue
        end
        S = load(fullfile(flist(f).folder, flist(f).name), 'parsing_struct');
        if ~isfield(S, 'parsing_struct')
            warning('No parsing_struct in %s', flist(f).name);
            continue
        end
        for k = 1:numel(S.parsing_struct)
            code = char(string(S.parsing_struct(k).animal_code));
            if ismember(code, loaded_codes)
                disp(['  already in pooled file, skip ', code, ' (', flist(f).name, ')'])
                continue
            end
            parsing_all{end + 1} = S.parsing_struct(k); %#ok<AGROW>
            loaded_codes{end + 1} = code; %#ok<AGROW>
            loaded_from{end + 1} = flist(f).name; %#ok<AGROW>
            disp(['  loaded ', code, ' from ', flist(f).name])
        end
        clear S
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
    mu = mean(X, 2, 'omitnan');
    sd = std(X, 0, 2, 'omitnan');
    sd(sd == 0 | isnan(sd)) = 1;
    Z = (X - mu) ./ sd;
end


function plot_mean_sem(t, M, col)
    if isempty(M)
        return
    end
    mu = mean(M, 1, 'omitnan');
    sem = std(M, 0, 1, 'omitnan') / sqrt(max(sum(all(isfinite(M), 2)), 1));
    fill([t fliplr(t)], [mu - sem, fliplr(mu + sem)], col, 'FaceAlpha', 0.25, 'EdgeColor', 'none')
    hold on
    plot(t, mu, 'Color', col, 'LineWidth', 2)
end
