function result = locking_acg_within_neuron(source, options)
%LOCKING_ACG_WITHIN_NEURON  Common within-neuron LFP slope, by lock group.
%
% source is the locking summary, its power_strata struct, or the path of
% a summary file. options.min_strata is the number of finite strata a
% neuron needs. options.min_clean is how many must remain after a stratum
% with modified z-score above 3.5 is removed. The common slope and the
% per-neuron slopes use those cleaned strata. Spearman correlations keep
% every finite stratum.
%
% The common slope is one fit. Neuron intercepts are removed by
% subtracting each neuron's mean, which is the same as giving every neuron
% its own intercept. Cluster-robust standard errors are grouped by neuron.

options = fill_options(options);
result = empty_result(options);
[summary, strata] = unpack_source(source);
assert_requested_method(summary, options);
if ~strata_have_rows(strata)
    disp('Power strata are not in this summary. Run Analyze_locking_period_acg once.')
    return
end
check_strata_fields(strata);
if ~(options.min_strata >= 2)
    error('locking_acg_within_neuron:minStrata', ...
        'min_strata must be at least 2.');
end
if ~(options.min_clean >= 4)
    error('locking_acg_within_neuron:minClean', ...
        ['min_clean must be at least 4 so a within-neuron slope that ' ...
        'includes firing rate has a residual degree of freedom.']);
end

diary_cleanup = start_report(options.figure_folder);
print_model_header(options);
for o = 1:numel(result.outcome_fields)
    block = fit_one_outcome(strata, result.outcome_fields{o}, options, result);
    if o == 1
        result.per_neuron = block.per_neuron;
        result.model = block.model;
    else
        result.per_neuron(o) = block.per_neuron;
        result.model(o) = block.model;
    end
    report_outcome(result, o);
end
result = summarize_paired(result);
report_paired_change(result);
result.aligned = aligned_properties(strata, result);
result = summarize_examples(result);
report_examples(result);
clear diary_cleanup
if options.do_plot
    plot_violins(result, options, 'slope', 'slope_p', ...
        'Neuron slopes, violins', 'Within-neuron slope');
    plot_violins(result, options, 'rho', 'rho_p', ...
        'Neuron slopes, rank violins', 'Spearman rho');
    plot_violins(result, options, 'partial_slope', 'partial_slope_p', ...
        'Neuron slopes after firing rate, violins', 'Slope after firing rate');
    plot_violins(result, options, 'partial_rho', 'partial_rho_p', ...
        'Neuron rank correlations after firing rate, violins', ...
        'Partial Spearman rho');
    plot_median_r_bars(result, options);
    plot_slope_change(result, options);
    plot_common_slopes(result, options);
    plot_trough_examples(result, options, 'median_index', ...
        'Trough correlations near the median');
    plot_trough_examples(result, options, 'significant_index', ...
        'Trough significant correlations');
    plot_rho_matrices(result, options);
    plot_trough_property_link(result, options);
end
end


function options = fill_options(options)
defaults = struct( ...
    'min_strata', 10, ...
    'min_clean', 6, ...
    'alpha', 0.05, ...
    'do_plot', true, ...
    'figure_folder', '', ...
    'figure_visible', 'on', ...
    'acg_method', '');
if nargin < 1 || isempty(options)
    options = defaults;
    return
end
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(options, names{i}) || isempty(options.(names{i}))
        options.(names{i}) = defaults.(names{i});
    end
end
end


function result = empty_result(options)
result.min_strata = options.min_strata;
result.min_clean = options.min_clean;
result.alpha = options.alpha;
result.outcome_fields = {'acg_delta', 'first_prominence', 'second_prominence'};
result.outcome_titles = {'ACG delta power', 'First-peak prominence', ...
    'Second-cycle prominence'};
result.type_order = {'trough', 'peak', 'nonentrained'};
result.type_title = {'Trough-locked', 'Peak-locked', 'Non-entrained'};
result.panel_names = {'all', 'trough', 'peak', 'delta', 'nonentrained'};
result.panel_titles = {'All neurons', 'Trough-locked', 'Peak-locked', ...
    'Delta-locked', 'Non-entrained'};
result.cell_order = {'trough', 'peak', 'delta', 'nonentrained'};
result.cell_title = {'Trough-locked', 'Peak-locked', 'Delta-locked', ...
    'Non-entrained'};
result.per_neuron = struct([]);
result.model = struct([]);
end


function [summary, strata] = unpack_source(source)
summary = struct();
if ischar(source)
    loaded = load(source, 'locking_acg');
    if ~isfield(loaded, 'locking_acg')
        error('locking_acg_within_neuron:summary', ...
            'The file %s has no locking summary.', source);
    end
    source = loaded.locking_acg;
end
if isstruct(source) && isfield(source, 'power_strata')
    summary = source;
    strata = source.power_strata;
else
    strata = source;
end
end


function assert_requested_method(summary, options)
if ~isfield(options, 'acg_method') || isempty(options.acg_method)
    return
end
if ~isstruct(summary) || ~isfield(summary, 'parameters')
    return
end
saved = 'overlap';
if isfield(summary.parameters, 'acg_method') && ~isempty(summary.parameters.acg_method)
    saved = canonical_acg_method(summary.parameters.acg_method);
end
requested = canonical_acg_method(options.acg_method);
if ~strcmp(saved, requested)
    error('locking_acg_within_neuron:method', ...
        ['This summary was saved with acg_method ''%s'', and the ' ...
        'script is set to ''%s''.'], saved, requested);
end
end


function tf = strata_have_rows(strata)
tf = isstruct(strata) && isfield(strata, 'neuron') && ~isempty(strata.neuron) ...
    && isfield(strata, 'power');
end


function check_strata_fields(strata)
required = {'neuron', 'lock_group', 'power', 'rate', 'acg_delta', ...
    'first_prominence', 'second_prominence'};
for i = 1:numel(required)
    if ~isfield(strata, required{i})
        error('locking_acg_within_neuron:strata', ...
            'Power strata have no %s column.', required{i});
    end
end
end


function block = fit_one_outcome(strata, outcome_field, options, result)
neuron = as_cellstr(strata.neuron);
lock_group = as_cellstr(strata.lock_group);
power = double(strata.power(:));
rate = double(strata.rate(:));
outcome = double(strata.(outcome_field)(:));
[rec, n_other] = neuron_records(neuron, lock_group, power, rate, outcome, ...
    options, result.type_order);
per_neuron = blank_panels(result);
for i = 1:numel(rec)
    for p = 1:numel(result.panel_names)
        if ~panel_has(rec(i).group, result.panel_names{p})
            continue
        end
        per_neuron = append_record(per_neuron, p, rec(i));
    end
end
per_neuron.n_other_group = n_other;
model = fit_common_models(rec, outcome_field, result);
block.per_neuron = per_neuron;
block.model = model;
end


function [rec, n_other] = neuron_records(neuron, lock_group, power, rate, outcome, options, type_order)
ids = unique(neuron, 'stable');
rec = struct('id', {}, 'group', {}, 'slope', {}, 'slope_p', {}, ...
    'partial_slope', {}, 'partial_slope_p', {}, 'rho', {}, 'rho_p', {}, ...
    'partial_rho', {}, 'partial_rho_p', {}, 'n_dropped', {}, ...
    'P', {}, 'y', {}, 'rate', {}, 'P_corr', {}, 'y_corr', {});
n_other = 0;
n_rec = 0;
for i = 1:numel(ids)
    row_idx = find(strcmp(neuron, ids{i}));
    group = lock_group{row_idx(1)};
    if ~any(strcmp(group, type_order))
        n_other = n_other + 1;
        continue
    end
    P = power(row_idx);
    r = rate(row_idx);
    y = outcome(row_idx);
    ok = isfinite(P) & isfinite(r) & isfinite(y);
    P = P(ok);
    r = r(ok);
    y = y(ok);
    if numel(P) < options.min_strata
        continue
    end
    [rho, rho_p] = spearman_with_p(P, y);
    [partial_rho, partial_rho_p] = partial_spearman(P, y, r);
    [slope, slope_p, partial_slope, partial_slope_p, n_drop, Pc, yc, rc] = ...
        cleaned_slopes(P, y, r, options.min_clean);
    if ~isfinite(rho) && ~isfinite(slope)
        continue
    end
    n_rec = n_rec + 1;
    rec(n_rec).id = ids{i}; %#ok<AGROW>
    rec(n_rec).group = group;
    rec(n_rec).slope = slope;
    rec(n_rec).slope_p = slope_p;
    rec(n_rec).partial_slope = partial_slope;
    rec(n_rec).partial_slope_p = partial_slope_p;
    rec(n_rec).rho = rho;
    rec(n_rec).rho_p = rho_p;
    rec(n_rec).partial_rho = partial_rho;
    rec(n_rec).partial_rho_p = partial_rho_p;
    rec(n_rec).n_dropped = n_drop;
    rec(n_rec).P = Pc;
    rec(n_rec).y = yc;
    rec(n_rec).rate = rc;
    rec(n_rec).P_corr = P;
    rec(n_rec).y_corr = y;
end
end


function [slope, slope_p, partial, partial_p, n_drop, x, y, rate] = ...
        cleaned_slopes(P, y_all, rate_all, min_clean)
slope = NaN;
slope_p = NaN;
partial = NaN;
partial_p = NaN;
x = zeros(0, 1);
y = zeros(0, 1);
rate = zeros(0, 1);
bad = modified_z_rows([P(:), rate_all(:), y_all(:)]);
n_drop = nnz(bad);
x = P(~bad);
y = y_all(~bad);
rate = rate_all(~bad);
if numel(x) < min_clean || column_std(x) <= 0 || column_std(y) <= 0
    x = zeros(0, 1);
    y = zeros(0, 1);
    rate = zeros(0, 1);
    return
end
[slope, slope_p] = simple_slope(x, y);
if column_std(rate) <= 0
    partial = slope;
    partial_p = slope_p;
    return
end
[partial, partial_p] = rate_adjusted_slope(x, y, rate);
end


function [slope, p] = simple_slope(x, y)
slope = NaN;
p = NaN;
n = numel(x);
A = [ones(n, 1), x(:)];
coef = A \ y(:);
resid = y(:) - A * coef;
df = n - 2;
if df < 1
    return
end
ss = sum(resid .^ 2) / df;
C = pinv(A' * A);
se = sqrt(max(ss * C(2, 2), 0));
slope = coef(2);
if ~isfinite(se)
    slope = NaN;
    return
end
if se <= 0
    p = double(slope == 0);
    return
end
p = t_two_tailed(slope / se, df);
end


function [slope, p] = rate_adjusted_slope(x, y, rate)
slope = NaN;
p = NaN;
n = numel(x);
A = [ones(n, 1), x(:), rate(:)];
coef = A \ y(:);
resid = y(:) - A * coef;
df = n - 3;
if df < 1
    return
end
ss = sum(resid .^ 2) / df;
covb = ss * pinv(A' * A);
se = sqrt(max(covb(2, 2), 0));
slope = coef(2);
if ~isfinite(se)
    slope = NaN;
    return
end
if se <= 0
    p = double(slope == 0);
    return
end
p = t_two_tailed(slope / se, df);
end


function [rho, p] = spearman_with_p(x, y)
rho = rank_corr(x, y);
df = numel(x) - 2;
p = corr_p(rho, df);
end


function [rho, p] = partial_spearman(x, y, z)
rho = NaN;
p = NaN;
n = numel(x);
if n < 4
    return
end
rx = average_rank(x);
ry = average_rank(y);
rz = average_rank(z);
rxy = pearson_corr(rx, ry);
rxz = pearson_corr(rx, rz);
ryz = pearson_corr(ry, rz);
den = sqrt((1 - rxz ^ 2) * (1 - ryz ^ 2));
if ~isfinite(rxy) || ~isfinite(den) || den <= 0
    return
end
rho = (rxy - rxz * ryz) / den;
rho = max(min(rho, 1), -1);
p = corr_p(rho, n - 3);
end


function per_neuron = blank_panels(result)
n_panel = numel(result.panel_names);
per_neuron.slope = empty_cells(n_panel);
per_neuron.slope_p = empty_cells(n_panel);
per_neuron.partial_slope = empty_cells(n_panel);
per_neuron.partial_slope_p = empty_cells(n_panel);
per_neuron.rho = empty_cells(n_panel);
per_neuron.rho_p = empty_cells(n_panel);
per_neuron.partial_rho = empty_cells(n_panel);
per_neuron.partial_rho_p = empty_cells(n_panel);
per_neuron.n_dropped = empty_cells(n_panel);
per_neuron.slope_drop = empty_cells(n_panel);
per_neuron.rho_drop = empty_cells(n_panel);
per_neuron.paired_slope_id = empty_id_cells(n_panel);
per_neuron.rho_id = empty_id_cells(n_panel);
per_neuron.rho_P = empty_id_cells(n_panel);
per_neuron.rho_y = empty_id_cells(n_panel);
per_neuron.n_other_group = 0;
end


function cells = empty_cells(n_panel)
cells = cell(1, n_panel);
for p = 1:n_panel
    cells{p} = zeros(0, 1);
end
end


function cells = empty_id_cells(n_panel)
cells = cell(1, n_panel);
for p = 1:n_panel
    cells{p} = {};
end
end


function per_neuron = append_record(per_neuron, p, rec)
if isfinite(rec.rho)
    per_neuron.rho{p}(end + 1, 1) = rec.rho;
    per_neuron.rho_p{p}(end + 1, 1) = rec.rho_p;
    per_neuron.rho_id{p}{end + 1, 1} = rec.id;
    per_neuron.rho_P{p}{end + 1, 1} = rec.P_corr;
    per_neuron.rho_y{p}{end + 1, 1} = rec.y_corr;
end
if isfinite(rec.partial_rho)
    per_neuron.partial_rho{p}(end + 1, 1) = rec.partial_rho;
    per_neuron.partial_rho_p{p}(end + 1, 1) = rec.partial_rho_p;
end
if isfinite(rec.slope)
    per_neuron.slope{p}(end + 1, 1) = rec.slope;
    per_neuron.slope_p{p}(end + 1, 1) = rec.slope_p;
    per_neuron.n_dropped{p}(end + 1, 1) = rec.n_dropped;
end
if isfinite(rec.partial_slope)
    per_neuron.partial_slope{p}(end + 1, 1) = rec.partial_slope;
    per_neuron.partial_slope_p{p}(end + 1, 1) = rec.partial_slope_p;
end
if isfinite(rec.slope) && isfinite(rec.partial_slope)
    per_neuron.slope_drop{p}(end + 1, 1) = rec.slope - rec.partial_slope;
    per_neuron.paired_slope_id{p}{end + 1, 1} = rec.id;
end
if isfinite(rec.rho) && isfinite(rec.partial_rho)
    per_neuron.rho_drop{p}(end + 1, 1) = rec.rho - rec.partial_rho;
end
end


function model = fit_common_models(rec, outcome_field, result)
[P, y, rate, neuron, group] = cleaned_rows(rec);
model.outcome = outcome_field;
model.simple = cluster_group_slopes(P, y, rate, neuron, group, ...
    result.type_order, result.type_title, false);
model.with_rate = cluster_group_slopes(P, y, rate, neuron, group, ...
    result.type_order, result.type_title, true);
model.share = slope_share(P, y, rate, neuron, group, result, model);
end


function [P, y, rate, neuron, group] = cleaned_rows(rec)
P = zeros(0, 1);
y = zeros(0, 1);
rate = zeros(0, 1);
neuron = {};
group = {};
for i = 1:numel(rec)
    if ~isfinite(rec(i).slope) || isempty(rec(i).P)
        continue
    end
    n = numel(rec(i).P);
    P = [P; rec(i).P(:)]; %#ok<AGROW>
    y = [y; rec(i).y(:)]; %#ok<AGROW>
    rate = [rate; rec(i).rate(:)]; %#ok<AGROW>
    neuron = [neuron; repmat({rec(i).id}, n, 1)]; %#ok<AGROW>
    group = [group; repmat({rec(i).group}, n, 1)]; %#ok<AGROW>
end
end


function fit = cluster_group_slopes(P, y, rate, neuron, group, type_order, type_title, include_rate)
% First coefficients are the LFP slopes. When firing rate is included,
% those slopes come next, in the same group order.
n_type = numel(type_order);
names = cell(1, n_type);
for g = 1:n_type
    names{g} = ['LFP delta power, ', type_title{g}];
end
if include_rate
    for g = 1:n_type
        names{end + 1} = ['Firing rate, ', type_title{g}]; %#ok<AGROW>
    end
end
fit = blank_fit(names);
fit.include_rate = include_rate;
if numel(y) < 4 || isempty(neuron)
    fit.reason = 'Too few strata remain.';
    return
end
[demeaned, neuron] = demean_columns([P(:), rate(:), y(:)], neuron);
Pc = demeaned(:, 1);
Rc = demeaned(:, 2);
yc = demeaned(:, 3);
n_col = n_type * (1 + double(include_rate));
Xc = zeros(numel(yc), n_col);
group = group(:);
for g = 1:n_type
    in_group = strcmp(group, type_order{g});
    Xc(:, g) = Pc .* double(in_group);
    if include_rate
        Xc(:, n_type + g) = Rc .* double(in_group);
    end
end
fit = fit_cluster(Xc, yc, neuron, names, fit);
fit = add_power_contrasts(fit, n_type);
fit.n_by_group = zeros(1, n_type);
fit.n_strata_by_group = zeros(1, n_type);
ids = unique(neuron(:), 'stable');
for g = 1:n_type
    in_group = strcmp(group, type_order{g});
    fit.n_strata_by_group(g) = nnz(in_group);
    n_neurons = 0;
    for i = 1:numel(ids)
        if any(in_group & strcmp(neuron(:), ids{i}))
            n_neurons = n_neurons + 1;
        end
    end
    fit.n_by_group(g) = n_neurons;
end
end


function fit = blank_fit(names)
n = numel(names);
fit.ok = false;
fit.names = names;
fit.beta = nan(n, 1);
fit.se = nan(n, 1);
fit.t = nan(n, 1);
fit.p = nan(n, 1);
fit.V = nan(n, n);
fit.df = NaN;
fit.n = 0;
fit.n_neurons = 0;
fit.rsq = NaN;
fit.reason = 'The common slope was not fit.';
fit.include_rate = false;
fit.n_by_group = zeros(1, 0);
fit.n_strata_by_group = zeros(1, 0);
fit.contrast_name = {};
fit.contrast_beta = zeros(0, 1);
fit.contrast_se = zeros(0, 1);
fit.contrast_p = zeros(0, 1);
fit.contrast_p_above = zeros(0, 1);
end


function fit = fit_cluster(X, y, neuron, names, fit)
spread = std(X, 0, 1);
keep = isfinite(spread) & spread > 1e-12;
if ~any(keep) || column_std(y) <= 0
    fit.reason = 'No predictor varies within a neuron.';
    return
end
Xk = X(:, keep);
beta_k = Xk \ y(:);
resid = y(:) - Xk * beta_k;
ids = unique(neuron(:), 'stable');
G = numel(ids);
p = size(Xk, 2);
bread = pinv(Xk' * Xk);
meat = zeros(p, p);
for i = 1:G
    in = strcmp(neuron(:), ids{i});
    score = Xk(in, :)' * resid(in);
    meat = meat + score * score';
end
if G > 1
    meat = meat * (G / (G - 1));
end
V_k = bread * meat * bread;
V_k = (V_k + V_k') / 2;
se_k = sqrt(max(real(diag(V_k)), 0));
df = G - p;
t_k = beta_k ./ se_k;
t_k(~isfinite(se_k) | se_k <= 0) = NaN;
p_k = t_two_tailed(t_k, df);
kept = find(keep);
fit.ok = true;
fit.reason = '';
fit.beta(kept) = beta_k;
fit.se(kept) = se_k;
fit.t(kept) = t_k;
fit.p(kept) = p_k;
fit.V(kept, kept) = V_k;
fit.df = df;
fit.n = numel(y);
fit.n_neurons = G;
fit.names = names;
ss = sum(y(:) .^ 2);
if ss > 0
    fit.rsq = 1 - sum(resid .^ 2) / ss;
end
end


function fit = add_power_contrasts(fit, n_type)
pairs = [1 2; 1 3];
labels = {'Trough minus peak', 'Trough minus non-entrained'};
if n_type < 3
    pairs = pairs(pairs(:, 2) <= n_type, :);
    labels = labels(1:size(pairs, 1));
end
n_pairs = size(pairs, 1);
fit.contrast_name = labels;
fit.contrast_beta = nan(n_pairs, 1);
fit.contrast_se = nan(n_pairs, 1);
fit.contrast_p = nan(n_pairs, 1);
fit.contrast_p_above = nan(n_pairs, 1);
if ~fit.ok
    return
end
for k = 1:n_pairs
    i = pairs(k, 1);
    j = pairs(k, 2);
    if ~isfinite(fit.beta(i)) || ~isfinite(fit.beta(j))
        continue
    end
    difference = fit.beta(i) - fit.beta(j);
    se = sqrt(max(fit.V(i, i) + fit.V(j, j) - 2 * fit.V(i, j), 0));
    fit.contrast_beta(k) = difference;
    fit.contrast_se(k) = se;
    if isfinite(se) && se > 0
        fit.contrast_p(k) = t_two_tailed(difference / se, fit.df);
        fit.contrast_p_above(k) = one_sided_above(difference, fit.contrast_p(k));
    end
end
end


function share = slope_share(P, y, rate, neuron, group, result, model)
n_type = numel(result.type_order);
share = blank_share(n_type);
if isempty(y)
    return
end
[demeaned, ~] = demean_columns([P(:), rate(:), y(:)], neuron);
Pc = demeaned(:, 1);
Rc = demeaned(:, 2);
yc = demeaned(:, 3);
for g = 1:n_type
    in_group = strcmp(group(:), result.type_order{g});
    r2_power = column_rsq(Pc(in_group), yc(in_group));
    r2_rate = column_rsq(Rc(in_group), yc(in_group));
    r2_both = matrix_rsq([Pc(in_group), Rc(in_group)], yc(in_group));
    unique_power = r2_both - r2_rate;
    fraction_remaining_r2 = NaN;
    fraction_due = NaN;
    if isfinite(r2_power) && r2_power > 1e-8 && isfinite(unique_power)
        fraction_remaining_r2 = unique_power / r2_power;
        fraction_due = 1 - fraction_remaining_r2;
    end
    beta_simple = model.simple.beta(g);
    beta_rate = model.with_rate.beta(g);
    fraction_slope = NaN;
    if isfinite(beta_simple) && abs(beta_simple) > 1e-8 && isfinite(beta_rate)
        fraction_slope = beta_rate / beta_simple;
    end
    share(g).r2_power = r2_power;
    share(g).r2_rate = r2_rate;
    share(g).r2_both = r2_both;
    share(g).fraction_remaining_r2 = fraction_remaining_r2;
    share(g).fraction_due_to_rate = fraction_due;
    share(g).beta_simple = beta_simple;
    share(g).beta_with_rate = beta_rate;
    share(g).fraction_of_slope_remaining = fraction_slope;
end
end


function share = blank_share(n_type)
share = struct('r2_power', {}, 'r2_rate', {}, 'r2_both', {}, ...
    'fraction_remaining_r2', {}, 'fraction_due_to_rate', {}, ...
    'beta_simple', {}, 'beta_with_rate', {}, 'fraction_of_slope_remaining', {});
for g = 1:n_type
    share(g).r2_power = NaN; %#ok<AGROW>
    share(g).r2_rate = NaN;
    share(g).r2_both = NaN;
    share(g).fraction_remaining_r2 = NaN;
    share(g).fraction_due_to_rate = NaN;
    share(g).beta_simple = NaN;
    share(g).beta_with_rate = NaN;
    share(g).fraction_of_slope_remaining = NaN;
end
end


function cleanup = start_report(folder)
cleanup = [];
if isempty(folder)
    return
end
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
report_path = fullfile(folder, 'within_neuron_acg_slopes.txt');
if exist(report_path, 'file') == 2
    delete(report_path);
end
diary(report_path);
cleanup = onCleanup(@() diary('off'));
end


function print_model_header(options)
disp('WITHIN-NEURON LFP SLOPE')
disp(['One fit uses every neuron. Each neuron has its own intercept. ' ...
    'The slope of LFP delta power is separate for trough-locked, ' ...
    'peak-locked, and non-entrained neurons.'])
disp(sprintf(['A neuron needs %d finite strata of LFP power, firing rate, ' ...
    'and the autocorrelogram measure. Slopes then drop strata with modified ' ...
    'z-score above 3.5 and keep a neuron when %d strata remain. Spearman ' ...
    'correlations use every finite stratum.'], options.min_strata, options.min_clean))
disp(['The second fit adds firing rate, with its own slope in each lock ' ...
    'group. The rate-adjusted LFP slope is the test of whether the ' ...
    'association remains. Delta-locked neurons appear in the per-neuron ' ...
    'distributions as trough and peak together. They are not a fourth slope.'])
end


function report_outcome(result, o)
disp(repmat('-', 1, 72))
disp(result.outcome_titles{o})
if o == 1 && result.per_neuron(o).n_other_group > 0
    disp(sprintf(['  %d neurons had a lock label outside trough, peak, and ' ...
        'non-entrained, and were left out.'], result.per_neuron(o).n_other_group))
end
disp('  Common slope, one fit for every neuron')
report_fit(result, o, result.model(o).simple, 'slope', 'Without firing rate');
report_fit(result, o, result.model(o).with_rate, 'partial_slope', 'With firing rate');
report_share(result, o);
disp('  One slope and one correlation per neuron')
for p = 1:numel(result.panel_titles)
    report_panel(result.panel_titles{p}, result.per_neuron(o), p, result.alpha);
end
end


function report_fit(result, o, fit, slope_field, label)
disp(['    ', label])
if ~fit.ok
    if isempty(fit.reason)
        disp('      The common slope was not fit.')
    else
        disp(['      ', fit.reason])
    end
    return
end
n_type = numel(result.type_order);
for g = 1:n_type
    median_slope = panel_median(result, o, slope_field, result.type_order{g});
    disp(sprintf(['      %s: beta = %.3g, SE = %.3g, two-sided p = %.3g, ' ...
        'one-sided above 0 p = %.3g, neurons = %d, strata = %d, ' ...
        'median neuron slope = %.3g'], ...
        result.type_title{g}, fit.beta(g), fit.se(g), fit.p(g), ...
        one_sided_above(fit.beta(g), fit.p(g)), ...
        fit.n_by_group(g), fit.n_strata_by_group(g), median_slope))
end
if fit.include_rate
    disp('      Firing-rate slopes in the same fit')
    for g = 1:n_type
        k = n_type + g;
        disp(sprintf('      %s: beta = %.3g, SE = %.3g, two-sided p = %.3g', ...
            result.type_title{g}, fit.beta(k), fit.se(k), fit.p(k)))
    end
end
for k = 1:numel(fit.contrast_name)
    disp(sprintf(['      %s: difference = %.3g, SE = %.3g, two-sided p = %.3g, ' ...
        'one-sided trough larger p = %.3g'], ...
        fit.contrast_name{k}, fit.contrast_beta(k), fit.contrast_se(k), ...
        fit.contrast_p(k), fit.contrast_p_above(k)))
end
disp(sprintf(['      Within R^2 = %.3g, neurons = %d, strata = %d, ' ...
    'cluster df = %.3g'], fit.rsq, fit.n_neurons, fit.n, fit.df))
end


function value = panel_median(result, o, field_name, panel_name)
k = find(strcmp(result.panel_names, panel_name), 1);
value = NaN;
if isempty(k)
    return
end
value = finite_median(result.per_neuron(o).(field_name){k});
end


function report_share(result, o)
disp('  How much of the association firing rate accounts for')
share = result.model(o).share;
simple = result.model(o).simple;
for g = 1:numel(result.type_title)
    use_fraction = isfinite(share(g).fraction_of_slope_remaining) ...
        && (~isfinite(simple.p(g)) || simple.p(g) < 0.05);
    if use_fraction
        disp(sprintf(['    %s: slope without rate = %.3g, slope with rate = %.3g, ' ...
            'fraction of the slope that remains = %.3g'], ...
            result.type_title{g}, share(g).beta_simple, share(g).beta_with_rate, ...
            share(g).fraction_of_slope_remaining))
    else
        disp(sprintf(['    %s: slope without rate = %.3g, slope with rate = %.3g. ' ...
            'The simple slope is not distinguishable from zero, so the ' ...
            'fraction removed by firing rate is not used.'], ...
            result.type_title{g}, share(g).beta_simple, share(g).beta_with_rate))
    end
    stable_r2 = isfinite(share(g).r2_power) && share(g).r2_power >= 0.02 ...
        && isfinite(share(g).fraction_remaining_r2) ...
        && share(g).fraction_remaining_r2 >= -0.05 ...
        && share(g).fraction_remaining_r2 <= 1.05;
    if use_fraction && stable_r2
        disp(sprintf(['      Within R^2 of LFP power = %.3g, of firing rate = %.3g, ' ...
            'of both = %.3g'], share(g).r2_power, share(g).r2_rate, share(g).r2_both))
        disp(sprintf(['      Fraction of the LFP-power R^2 removed by firing rate ' ...
            '= %.3g. Fraction of that R^2 that remains = %.3g.'], ...
            share(g).fraction_due_to_rate, share(g).fraction_remaining_r2))
    elseif use_fraction && ~isfinite(share(g).r2_rate)
        disp(['      Firing rate does not vary within these neurons, so the ' ...
            'LFP slope is unchanged.'])
    elseif use_fraction
        disp(['      The within-neuron R^2 of LFP power is too small for a ' ...
            'stable split with firing rate.'])
    end
    if use_fraction && isfinite(share(g).fraction_of_slope_remaining) ...
            && share(g).fraction_of_slope_remaining < 0
        disp('      The LFP slope changes sign once firing rate is included.')
    end
end
disp(['    The fraction of the slope is the rate-adjusted beta divided by the ' ...
    'simple beta. The R^2 fraction is the part of the within-neuron fit that ' ...
    'LFP power still adds after firing rate. It is smaller when firing rate ' ...
    'and LFP power rise together, even if the slope itself remains.'])
end


function report_panel(label, per, p, alpha)
disp(['    ', label])
report_value_line('Slopes', per.slope{p}, per.slope_p{p}, per.n_dropped{p}, alpha, true);
report_value_line('Rate-adjusted slopes', per.partial_slope{p}, ...
    per.partial_slope_p{p}, [], alpha, false);
report_shift('Simple slope minus rate-adjusted slope', per.slope_drop{p});
report_value_line('Spearman', per.rho{p}, per.rho_p{p}, [], alpha, false);
report_value_line('Partial Spearman given firing rate', per.partial_rho{p}, ...
    per.partial_rho_p{p}, [], alpha, false);
report_shift('Simple Spearman minus partial Spearman', per.rho_drop{p});
end


function report_value_line(label, values, pvals, n_dropped, alpha, with_drop)
[n_pos, n_neg] = signed_counts(values, pvals, alpha);
extra = '';
if with_drop && ~isempty(n_dropped)
    extra = sprintf(', strata removed as outliers = %d', sum(n_dropped));
end
disp(sprintf(['      %s: n = %d, median = %.3g%s, significant positive = %d, ' ...
    'significant negative = %d'], label, nnz(isfinite(values)), ...
    finite_median(values), extra, n_pos, n_neg))
disp(sprintf(['        Versus 0: two-sided signrank p = %.3g, ' ...
    'one-sided above 0 p = %.3g'], ...
    signed_rank_p(values, 'both'), signed_rank_p(values, 'right')))
end


function report_shift(label, values)
[mean_change, paired_p] = mean_vs_zero(values);
disp(sprintf(['      %s: n = %d, median change = %.3g, mean change = %.3g, ' ...
    'paired t p = %.3g, signrank p = %.3g'], ...
    label, nnz(isfinite(values)), finite_median(values), mean_change, ...
    paired_p, signed_rank_p(values, 'both')))
end


function result = summarize_paired(result)
% Median Spearman rho by cell type, and the paired slope change when
% firing rate is added. A bar is filled when the mean differs from zero.
n_out = numel(result.outcome_fields);
n_cell = numel(result.cell_order);
panel = zeros(1, n_cell);
for g = 1:n_cell
    panel(g) = find(strcmp(result.panel_names, result.cell_order{g}), 1);
end
result.cell_panel = panel;
for o = 1:n_out
    bars.median_r = nan(1, n_cell);
    bars.mean_r = nan(1, n_cell);
    bars.p_mean = nan(1, n_cell);
    bars.filled = false(1, n_cell);
    change.median = nan(1, n_cell);
    change.mean = nan(1, n_cell);
    change.p_mean = nan(1, n_cell);
    change.p_signrank = nan(1, n_cell);
    change.filled = false(1, n_cell);
    change.n = zeros(1, n_cell);
    for g = 1:n_cell
        rho = result.per_neuron(o).rho{panel(g)};
        [bars.mean_r(g), bars.p_mean(g)] = mean_vs_zero(rho);
        bars.median_r(g) = finite_median(rho);
        bars.filled(g) = isfinite(bars.p_mean(g)) && bars.p_mean(g) < result.alpha;
        delta = result.per_neuron(o).slope_drop{panel(g)};
        [change.mean(g), change.p_mean(g)] = mean_vs_zero(delta);
        change.median(g) = finite_median(delta);
        change.p_signrank(g) = signed_rank_p(delta, 'both');
        change.n(g) = nnz(isfinite(delta));
        change.filled(g) = isfinite(change.p_mean(g)) && change.p_mean(g) < result.alpha;
    end
    if o == 1
        result.correlation_bars = bars;
        result.slope_change = change;
    else
        result.correlation_bars(o) = bars;
        result.slope_change(o) = change;
    end
end
for g = 1:n_cell
    p_var = nan(1, n_out);
    for o = 1:n_out
        p_var(o) = result.slope_change(o).p_mean(g);
    end
    finite_p = p_var(isfinite(p_var));
    joint.p_each = p_var;
    joint.p_any = NaN;
    joint.changes = false;
    if ~isempty(finite_p)
        joint.p_any = min(1, n_out * min(finite_p));
        joint.changes = isfinite(joint.p_any) && joint.p_any < result.alpha;
    end
    [joint.F, joint.p_hotelling, joint.df1, joint.df2, joint.n] = ...
        hotelling_slope_change(result, panel(g));
    if g == 1
        result.joint_slope_change = joint;
    else
        result.joint_slope_change(g) = joint;
    end
end
end


function report_paired_change(result)
disp(repmat('=', 1, 72))
disp('PAIRED CHANGE, WITHOUT FIRING RATE MINUS WITH FIRING RATE')
disp(['Each neuron contributes the difference of its two values. The paired ' ...
    't test asks whether the mean difference differs from zero.'])
for g = 1:numel(result.cell_title)
    disp(['  ', result.cell_title{g}])
    for o = 1:numel(result.outcome_titles)
        change = result.slope_change(o);
        disp(sprintf(['    %s slopes: n = %d, median change = %.3g, ' ...
            'mean change = %.3g, paired t p = %.3g, signrank p = %.3g'], ...
            result.outcome_titles{o}, change.n(g), change.median(g), ...
            change.mean(g), change.p_mean(g), change.p_signrank(g)))
    end
    joint = result.joint_slope_change(g);
    if joint.changes
        changed = 'yes';
    else
        changed = 'no';
    end
    disp(sprintf(['    Any slope changes when firing rate is added: %s ' ...
        '(Bonferroni p = %.3g).'], changed, joint.p_any))
    disp(sprintf(['    Joint Hotelling test of the three slope changes: ' ...
        'F = %.3g, df = %d %d, p = %.3g, neurons = %d'], ...
        joint.F, joint.df1, joint.df2, joint.p_hotelling, joint.n))
end
end


function [Fstat, p, df1, df2, n] = hotelling_slope_change(result, panel)
Fstat = NaN;
p = NaN;
df1 = numel(result.outcome_fields);
df2 = NaN;
n = 0;
ids = result.per_neuron(1).paired_slope_id{panel};
ids = ids(:);
for o = 2:numel(result.outcome_fields)
    ids = intersect_stable(ids, result.per_neuron(o).paired_slope_id{panel});
end
n = numel(ids);
pdim = numel(result.outcome_fields);
df1 = pdim;
df2 = n - pdim;
if n <= pdim
    return
end
D = nan(n, pdim);
for o = 1:pdim
    these = result.per_neuron(o).paired_slope_id{panel};
    drops = result.per_neuron(o).slope_drop{panel};
    for i = 1:n
        hit = find(strcmp(these(:), ids{i}), 1);
        if ~isempty(hit)
            D(i, o) = drops(hit);
        end
    end
end
if any(~isfinite(D(:)))
    return
end
mu = mean(D, 1)';
C = cov(D);
if rank(C) < pdim
    return
end
T2 = n * (mu' * (C \ mu));
Fstat = ((n - pdim) / (pdim * (n - 1))) * T2;
p = f_upper_tail(Fstat, pdim, n - pdim);
end


function kept = intersect_stable(first, second)
kept = {};
second = second(:);
for i = 1:numel(first)
    if any(strcmp(second, first{i}))
        kept{end + 1, 1} = first{i}; %#ok<AGROW>
    end
end
end


function [mean_value, p] = mean_vs_zero(values)
values = values(isfinite(values));
mean_value = NaN;
p = NaN;
n = numel(values);
if n < 1
    return
end
mean_value = mean(values);
if n < 2
    return
end
spread = std(values, 0);
if ~isfinite(spread) || spread <= 0
    p = double(mean_value == 0);
    return
end
tstat = mean_value / (spread / sqrt(n));
p = t_two_tailed(tstat, n - 1);
end


function p = f_upper_tail(Fstat, df1, df2)
p = NaN;
if ~isfinite(Fstat) || df1 < 1 || df2 < 1
    return
end
if Fstat < 0
    return
end
p = betainc(df2 / (df2 + df1 * Fstat), df2 / 2, df1 / 2);
end


function aligned = aligned_properties(strata, result)
neuron = as_cellstr(strata.neuron);
lock_group = as_cellstr(strata.lock_group);
power = double(strata.power(:));
rate = double(strata.rate(:));
Y = [double(strata.acg_delta(:)), double(strata.first_prominence(:)), ...
    double(strata.second_prominence(:))];
ids = unique(neuron, 'stable');
aligned.id = {};
aligned.group = {};
aligned.P = {};
aligned.Y = {};
n = 0;
for i = 1:numel(ids)
    rows = find(strcmp(neuron, ids{i}));
    group = lock_group{rows(1)};
    if ~any(strcmp(group, result.type_order))
        continue
    end
    P = power(rows);
    r = rate(rows);
    y = Y(rows, :);
    ok = isfinite(P) & isfinite(r);
    if nnz(ok) < 2
        continue
    end
    n = n + 1;
    aligned.id{n, 1} = ids{i};
    aligned.group{n, 1} = group;
    aligned.P{n, 1} = P(ok);
    aligned.Y{n, 1} = y(ok, :);
end
end


function result = summarize_examples(result)
n_pick = 5;
panel = find(strcmp(result.panel_names, 'trough'), 1);
examples.n_pick = n_pick;
examples.panel = panel;
n_out = numel(result.outcome_fields);
examples.median_index = cell(1, n_out);
examples.significant_index = cell(1, n_out);
for o = 1:n_out
    rho = result.per_neuron(o).rho{panel};
    pvals = result.per_neuron(o).rho_p{panel};
    examples.median_index{o} = nearest_to_value(rho, finite_median(rho), n_pick);
    examples.significant_index{o} = strongest_significant(rho, pvals, result.alpha, n_pick);
end
result.examples = examples;
result.rho_link = link_matrices(result);
end


function idx = nearest_to_value(values, target, n_pick)
ok = find(isfinite(values(:)));
idx = zeros(0, 1);
if isempty(ok) || ~isfinite(target)
    return
end
distance = abs(values(ok) - target);
[~, order] = sortrows([distance, ok]);
n_pick = min(n_pick, numel(order));
idx = ok(order(1:n_pick));
[~, show] = sort(values(idx));
idx = idx(show);
end


function idx = strongest_significant(values, pvals, alpha, n_pick)
values = values(:);
pvals = pvals(:);
ok = find(isfinite(values) & isfinite(pvals) & pvals < alpha);
idx = zeros(0, 1);
if isempty(ok)
    return
end
[~, order] = sortrows([pvals(ok), -abs(values(ok)), ok]);
n_pick = min(n_pick, numel(order));
idx = ok(order(1:n_pick));
[~, show] = sort(values(idx));
idx = idx(show);
end


function link = link_matrices(result)
n_cell = numel(result.cell_order);
n_out = numel(result.outcome_fields);
for g = 1:n_cell
    panel = result.cell_panel(g);
    [ids, R] = rho_table(result, panel);
    C = nan(n_out, n_out);
    n_pair = zeros(n_out, n_out);
    for a = 1:n_out
        for b = 1:n_out
            ok = isfinite(R(:, a)) & isfinite(R(:, b));
            n_pair(a, b) = nnz(ok);
            if a == b && n_pair(a, b) >= 2
                C(a, b) = 1;
            elseif n_pair(a, b) >= 4
                C(a, b) = pearson_corr(R(ok, a), R(ok, b));
            end
        end
    end
    block.ids = ids;
    block.R = R;
    block.C = C;
    block.n_pair = n_pair;
    if g == 1
        link = block;
    else
        link(g) = block; %#ok<AGROW>
    end
end
end


function [ids, R] = rho_table(result, panel)
n_out = numel(result.outcome_fields);
ids = {};
for o = 1:n_out
    ids = union_stable(ids, result.per_neuron(o).rho_id{panel});
end
R = nan(numel(ids), n_out);
for o = 1:n_out
    these = result.per_neuron(o).rho_id{panel};
    vals = result.per_neuron(o).rho{panel};
    for i = 1:numel(these)
        hit = find(strcmp(ids, these{i}), 1);
        if ~isempty(hit)
            R(hit, o) = vals(i);
        end
    end
end
end


function kept = union_stable(first, second)
kept = first(:);
second = second(:);
for i = 1:numel(second)
    if ~any(strcmp(kept, second{i}))
        kept{end + 1, 1} = second{i}; %#ok<AGROW>
    end
end
end


function report_examples(result)
panel = result.examples.panel;
disp(repmat('=', 1, 72))
disp('TROUGH EXAMPLE CORRELATIONS WITH LFP DELTA POWER')
disp(['Five neurons closest to the median Spearman correlation, and up to ' ...
    'five with a significant correlation. r is the within-neuron Spearman ' ...
    'correlation.'])
for o = 1:numel(result.outcome_titles)
    rho = result.per_neuron(o).rho{panel};
    pvals = result.per_neuron(o).rho_p{panel};
    ids = result.per_neuron(o).rho_id{panel};
    disp(['  ', result.outcome_titles{o}])
    disp(sprintf('    Median r = %s', plain_number(finite_median(rho))))
    report_example_list('Near the median', ids, rho, pvals, result.examples.median_index{o});
    report_example_list('Significant', ids, rho, pvals, result.examples.significant_index{o});
end
disp('CORRELATION OF THE PER-NEURON CORRELATIONS')
disp(['Each entry is the Pearson correlation, across neurons, of the ' ...
    'Spearman values for two autocorrelogram measures. The same neurons ' ...
    'pull both values up when the entry is positive.'])
short_outcome = {'ACG delta power', 'First-peak prominence', 'Second-cycle prominence'};
for g = 1:numel(result.cell_title)
    disp(['  ', result.cell_title{g}])
    C = result.rho_link(g).C;
    n_pair = result.rho_link(g).n_pair;
    for a = 1:numel(short_outcome)
        for b = a + 1:numel(short_outcome)
            disp(sprintf('    %s with %s: r = %s, neurons = %d', ...
                short_outcome{a}, short_outcome{b}, plain_number(C(a, b)), n_pair(a, b)))
        end
    end
end
end


function report_example_list(label, ids, rho, pvals, idx)
disp(['    ', label])
if isempty(idx)
    disp('      none')
    return
end
for i = 1:numel(idx)
    disp(sprintf('      %s  r = %s  p = %s', ids{idx(i)}, ...
        plain_number(rho(idx(i))), plain_number(pvals(idx(i)))))
end
end


function plot_trough_examples(result, options, index_field, figure_name)
n_out = numel(result.outcome_titles);
n_col = result.examples.n_pick;
panel = result.examples.panel;
short_outcome = {'ACG delta power', 'First-peak prominence', 'Second-cycle prominence'};
figure('Name', figure_name, 'Color', 'w', 'Visible', options.figure_visible)
set(gcf, 'Units', 'pixels', 'Position', [20 20 1680 980])
for o = 1:n_out
    idx = result.examples.(index_field){o};
    for c = 1:n_col
        subplot(n_out, n_col, (o - 1) * n_col + c)
        if c > numel(idx)
            axis off
            if c == 1
                title(sprintf('%s  none', short_outcome{o}))
            end
            continue
        end
        k = idx(c);
        draw_example_scatter(result.per_neuron(o).rho_P{panel}{k}, ...
            result.per_neuron(o).rho_y{panel}{k})
        title(sprintf('%s  r=%s', short_label(result.per_neuron(o).rho_id{panel}{k}), ...
            plain_number(result.per_neuron(o).rho{panel}(k))))
        if c == 1
            ylabel(short_outcome{o})
        end
        if o == n_out
            xlabel('LFP delta power')
        end
        box off
    end
end
save_figure(options.figure_folder, figure_name);
end


function draw_example_scatter(x, y)
x = x(:);
y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);
plot(x, y, '.', 'Color', [0.20 0.40 0.70], 'MarkerSize', 10)
hold on
if numel(x) >= 2 && std(x, 0) > 0 && std(y, 0) > 0
    coef = [ones(numel(x), 1), x] \ y;
    xx = [min(x); max(x)];
    plot(xx, coef(1) + coef(2) * xx, '-', 'Color', [0.75 0.12 0.14], 'LineWidth', 1.2)
end
end


function label = short_label(id)
label = char(id);
label = strrep(label, ',', ' ');
label = strrep(label, '-', ' ');
if numel(label) > 18
    label = label(end - 17:end);
end
end


function plot_rho_matrices(result, options)
n_cell = numel(result.cell_title);
short_outcome = {'ACG', 'Peak', 'Cycle'};
figure('Name', 'Correlation values by cell type', 'Color', 'w', ...
    'Visible', options.figure_visible)
set(gcf, 'Units', 'pixels', 'Position', [40 40 1400 420])
for g = 1:n_cell
    subplot(1, n_cell, g)
    draw_corr_matrix(result.rho_link(g).C, short_outcome)
    title(result.cell_title{g})
end
save_figure(options.figure_folder, 'Correlation values by cell type');
end


function draw_corr_matrix(C, labels)
hold on
n = size(C, 1);
for i = 1:n
    for j = 1:n
        patch([j - 0.5, j + 0.5, j + 0.5, j - 0.5], ...
            [i - 0.5, i - 0.5, i + 0.5, i + 0.5], diverging_color(C(i, j)), ...
            'EdgeColor', [0.85 0.85 0.85])
        if isfinite(C(i, j))
            text(j, i, sprintf('%.2f', C(i, j)), 'HorizontalAlignment', 'center', ...
                'FontSize', 10)
        else
            text(j, i, 'n/a', 'HorizontalAlignment', 'center', 'FontSize', 10, ...
                'Color', [0.35 0.35 0.35])
        end
    end
end
set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'YTick', 1:n, 'YTickLabel', labels, ...
    'YDir', 'reverse')
xlim([0.5, n + 0.5])
ylim([0.5, n + 0.5])
axis square
box on
end


function color = diverging_color(value)
if ~isfinite(value)
    color = [0.93 0.93 0.93];
    return
end
value = max(min(value, 1), -1);
if value >= 0
    color = [1, 1 - 0.75 * value, 1 - 0.75 * value];
else
    color = [1 + 0.75 * value, 1 + 0.55 * value, 1];
end
end


function plot_trough_property_link(result, options)
pairs = [1 2; 1 3; 2 3];
pair_names = {'ACG delta', 'First peak', 'Second cycle'};
panel = result.examples.panel;
figure('Name', 'Trough correlation value scatters', 'Color', 'w', ...
    'Visible', options.figure_visible)
set(gcf, 'Units', 'pixels', 'Position', [40 40 980 980])
for row = 1:size(pairs, 1)
    a = pairs(row, 1);
    b = pairs(row, 2);
    [xa, xb] = paired_rhos(result, panel, a, b);
    subplot(size(pairs, 1), 2, (row - 1) * 2 + 1)
    draw_pair_scatter(xa, xb, [0.20 0.45 0.80])
    title(sprintf('%s vs %s correlations  r=%s', pair_names{a}, pair_names{b}, ...
        plain_number(pairwise_pearson(xa, xb))))
    xlabel([pair_names{a}, ' correlation'])
    ylabel([pair_names{b}, ' correlation'])
    box off
    [ya, yb] = pooled_properties(result, a, b);
    subplot(size(pairs, 1), 2, (row - 1) * 2 + 2)
    draw_pair_scatter(ya, yb, [0.35 0.35 0.35])
    title(sprintf('%s vs %s values  r=%s', pair_names{a}, pair_names{b}, ...
        plain_number(pairwise_pearson(ya, yb))))
    xlabel([pair_names{a}, ' value'])
    ylabel([pair_names{b}, ' value'])
    box off
end
save_figure(options.figure_folder, 'Trough correlation value scatters');
end


function [xa, xb] = paired_rhos(result, panel, a, b)
ids = intersect_stable(result.per_neuron(a).rho_id{panel}, ...
    result.per_neuron(b).rho_id{panel});
xa = nan(numel(ids), 1);
xb = nan(numel(ids), 1);
ida = result.per_neuron(a).rho_id{panel};
idb = result.per_neuron(b).rho_id{panel};
ra = result.per_neuron(a).rho{panel};
rb = result.per_neuron(b).rho{panel};
for i = 1:numel(ids)
    ha = find(strcmp(ida, ids{i}), 1);
    hb = find(strcmp(idb, ids{i}), 1);
    if ~isempty(ha)
        xa(i) = ra(ha);
    end
    if ~isempty(hb)
        xb(i) = rb(hb);
    end
end
end


function [x, y] = pooled_properties(result, col_a, col_b)
panel = find(strcmp(result.panel_names, 'trough'), 1);
ids = result.per_neuron(1).rho_id{panel};
for o = 2:numel(result.outcome_fields)
    ids = union_stable(ids, result.per_neuron(o).rho_id{panel});
end
x = zeros(0, 1);
y = zeros(0, 1);
for i = 1:numel(result.aligned.id)
    if ~any(strcmp(ids, result.aligned.id{i}))
        continue
    end
    Y = result.aligned.Y{i};
    ok = isfinite(Y(:, col_a)) & isfinite(Y(:, col_b));
    x = [x; Y(ok, col_a)]; %#ok<AGROW>
    y = [y; Y(ok, col_b)]; %#ok<AGROW>
end
end


function draw_pair_scatter(x, y, color)
ok = isfinite(x(:)) & isfinite(y(:));
x = x(ok);
y = y(ok);
plot(x, y, '.', 'Color', color, 'MarkerSize', 8)
hold on
if numel(x) >= 2 && std(x, 0) > 0 && std(y, 0) > 0
    coef = [ones(numel(x), 1), x] \ y;
    xx = [min(x); max(x)];
    plot(xx, coef(1) + coef(2) * xx, '-', 'Color', [0.75 0.12 0.14], 'LineWidth', 1.2)
end
end


function r = pairwise_pearson(a, b)
ok = isfinite(a(:)) & isfinite(b(:));
if nnz(ok) < 4
    r = NaN;
    return
end
r = pearson_corr(a(ok), b(ok));
end


function plot_violins(result, options, value_field, p_field, figure_name, value_name)
n_out = numel(result.outcome_titles);
n_panel = numel(result.panel_titles);
figure('Name', figure_name, 'Color', 'w', 'Visible', options.figure_visible)
set(gcf, 'Units', 'pixels', 'Position', [40 40 1680 980])
short_panels = {'All', 'Trough', 'Peak', 'Delta', 'Non'};
for o = 1:n_out
    span = shared_value_span(result.per_neuron(o).(value_field));
    for p = 1:n_panel
        subplot(n_out, n_panel, (o - 1) * n_panel + p)
        values = result.per_neuron(o).(value_field){p};
        pvals = result.per_neuron(o).(p_field){p};
        draw_violin_scatter(values, pvals, result.alpha, span)
        title(sprintf('%s  median %s', short_panels{p}, ...
            plain_number(finite_median(values))))
        if p == 1
            ylabel(result.outcome_titles{o})
        end
        if o == n_out && p == 1
            xlabel(value_name)
        end
    end
end
if exist('sgtitle', 'file') == 2
    sgtitle([value_name, '. Red is significantly positive. Blue is significantly negative. Gray is not.'])
end
save_figure(options.figure_folder, figure_name);
end


function span = shared_value_span(panels)
values = zeros(0, 1);
for p = 1:numel(panels)
    values = [values; panels{p}(:)]; %#ok<AGROW>
end
values = values(isfinite(values));
if isempty(values)
    span = [-1 1];
    return
end
lo = min([values; 0]);
hi = max([values; 0]);
pad = 0.08 * max(hi - lo, 1);
span = [lo - pad, hi + pad];
end


function draw_violin_scatter(values, pvals, alpha, span)
hold on
plot([0.55 1.45], [0 0], '--', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
values = values(:);
pvals = pvals(:);
n = min(numel(values), numel(pvals));
values = values(1:n);
pvals = pvals(1:n);
ok = isfinite(values);
values = values(ok);
pvals = pvals(ok);
if numel(values) >= 2 && std(values, 0) > 0
    grid = linspace(span(1), span(2), 80)';
    scale = 0.12 * max(std(values, 0), eps);
    width = exp(-0.5 * (bsxfun(@minus, grid, values(:)') ./ scale) .^ 2);
    density = sum(width, 2);
    density = 0.32 * density / max(density);
    fill([1 - density; flipud(1 + density)], [grid; flipud(grid)], ...
        [0.75 0.80 0.86], 'EdgeColor', 'none', 'HandleVisibility', 'off')
end
for i = 1:numel(values)
    x = 1 + 0.12 * sin(i * 12.9898);
    if isfinite(pvals(i)) && pvals(i) < alpha && values(i) > 0
        col = [0.75 0.12 0.14];
    elseif isfinite(pvals(i)) && pvals(i) < alpha && values(i) < 0
        col = [0.20 0.40 0.72];
    else
        col = [0.45 0.45 0.45];
    end
    plot(x, values(i), '.', 'Color', col, 'MarkerSize', 10, 'HandleVisibility', 'off')
end
xlim([0.45 1.55])
ylim(span)
set(gca, 'XTick', [])
box off
end


function plot_median_r_bars(result, options)
n_out = numel(result.outcome_titles);
n_cell = numel(result.cell_title);
colors = cell_colors(result.cell_order);
short_cell = {'Trough', 'Peak', 'Delta', 'Non'};
figure('Name', 'Median Spearman by lock group', 'Color', 'w', ...
    'Visible', options.figure_visible)
set(gcf, 'Units', 'pixels', 'Position', [40 40 1400 480])
for o = 1:n_out
    subplot(1, n_out, o)
    hold on
    plot([0.4, n_cell + 0.6], [0 0], '--', 'Color', [0.35 0.35 0.35], ...
        'HandleVisibility', 'off')
    bars = result.correlation_bars(o);
    for g = 1:n_cell
        draw_sig_bar(g, bars.median_r(g), 0.7, colors(g, :), bars.filled(g))
    end
    set(gca, 'XTick', 1:n_cell, 'XTickLabel', short_cell)
    xlim([0.4, n_cell + 0.6])
    ylabel('Median Spearman rho')
    title(result.outcome_titles{o})
    box off
end
save_figure(options.figure_folder, 'Median Spearman by lock group');
end


function plot_slope_change(result, options)
n_out = numel(result.outcome_titles);
n_cell = numel(result.cell_title);
colors = [0.20 0.45 0.80; 0.85 0.45 0.15; 0.30 0.55 0.35];
short_outcome = {'ACG delta', 'First peak', 'Second cycle'};
short_cell = {'Trough', 'Peak', 'Delta', 'Non'};
figure('Name', 'Slope change after firing rate', 'Color', 'w', ...
    'Visible', options.figure_visible)
set(gcf, 'Units', 'pixels', 'Position', [40 40 920 560])
hold on
plot([0.4, n_cell + 0.6], [0 0], '--', 'Color', [0.35 0.35 0.35], ...
    'HandleVisibility', 'off')
offsets = [-0.24 0 0.24];
for o = 1:n_out
    for g = 1:n_cell
        draw_sig_bar(g + offsets(o), result.slope_change(o).median(g), 0.22, ...
            colors(o, :), result.slope_change(o).filled(g))
    end
    plot(nan, nan, 's', 'Color', colors(o, :), 'MarkerFaceColor', colors(o, :), ...
        'MarkerSize', 8)
end
y_top = ylim;
for g = 1:n_cell
    if result.joint_slope_change(g).changes
        plot(g, y_top(2), 'k*', 'MarkerSize', 8, 'HandleVisibility', 'off')
    end
end
set(gca, 'XTick', 1:n_cell, 'XTickLabel', short_cell)
xlim([0.4, n_cell + 0.6])
ylabel('Median slope change')
title('Filled when the mean slope change differs from zero')
legend(short_outcome, 'Location', 'northeast')
box off
save_figure(options.figure_folder, 'Slope change after firing rate');
end


function draw_sig_bar(x, height, width, color, filled)
if ~isfinite(height)
    return
end
h = bar(x, height, width);
if filled
    set(h, 'FaceColor', color, 'EdgeColor', color, 'LineWidth', 1)
else
    set(h, 'FaceColor', 'none', 'EdgeColor', color, 'LineWidth', 1.6)
end
set(h, 'HandleVisibility', 'off')
end


function label = plain_number(value)
if ~isfinite(value)
    label = 'NaN';
    return
end
label = sprintf('%.3g', abs(value));
if value < 0
    label = ['neg ', label];
end
end


function colors = cell_colors(names)
palette.trough = [0.20 0.45 0.80];
palette.peak = [0.85 0.25 0.20];
palette.delta = [0.10 0.50 0.35];
palette.nonentrained = [0.45 0.45 0.45];
colors = zeros(numel(names), 3);
for g = 1:numel(names)
    colors(g, :) = palette.(names{g});
end
end


function plot_common_slopes(result, options)
n_out = numel(result.outcome_titles);
n_type = numel(result.type_order);
simple_color = [0.20 0.45 0.80];
rate_color = [0.85 0.33 0.20];
short_title = {'Trough', 'Peak', 'Non-entrained'};
short_outcome = {'ACG delta power', 'First peak', 'Second cycle'};
figure('Name', 'Common within-neuron slopes', 'Color', 'w', ...
    'Visible', options.figure_visible)
set(gcf, 'Units', 'pixels', 'Position', [40 40 1500 520])
for o = 1:n_out
    subplot(1, n_out, o)
    hold on
    plot([0.5, n_type + 0.5], [0 0], '--', 'Color', [0.35 0.35 0.35], ...
        'HandleVisibility', 'off')
    draw_slope_points((1:n_type) - 0.12, result.model(o).simple, simple_color)
    draw_slope_points((1:n_type) + 0.12, result.model(o).with_rate, rate_color)
    set(gca, 'XTick', 1:n_type, 'XTickLabel', short_title)
    xlim([0.5, n_type + 0.5])
    title(short_outcome{o})
    ylabel('Common within-neuron slope')
    box off
    if o == 1
        plot(nan, nan, '--', 'Color', [0.35 0.35 0.35])
        plot(nan, nan, 'o', 'Color', simple_color, 'MarkerFaceColor', simple_color)
        plot(nan, nan, 'o', 'Color', rate_color, 'MarkerFaceColor', rate_color)
        legend({'Zero', 'Without firing rate', 'With firing rate'}, ...
            'Location', 'southeast')
    end
end
save_figure(options.figure_folder, 'Common within-neuron slopes');
end


function draw_slope_points(x, fit, color)
for g = 1:numel(x)
    if g > numel(fit.beta) || ~isfinite(fit.beta(g)) || ~isfinite(fit.se(g))
        continue
    end
    h = errorbar(x(g), fit.beta(g), 1.96 * fit.se(g), 'o');
    set(h, 'Color', color, 'MarkerFaceColor', color, 'LineWidth', 1.2)
    try
        set(h, 'HandleVisibility', 'off')
    catch
    end
end
end


function save_figure(folder, name)
if isempty(folder)
    return
end
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [16 9], 'PaperPosition', [0 0 16 9])
print(gcf, '-dpng', '-r120', fullfile(folder, [name, '.png']));
end


function [Xc, neuron] = demean_columns(X, neuron)
neuron = neuron(:);
Xc = X;
if isempty(X)
    return
end
ids = unique(neuron, 'stable');
for i = 1:numel(ids)
    in = strcmp(neuron, ids{i});
    Xc(in, :) = X(in, :) - mean(X(in, :), 1);
end
end


function s = column_std(x)
x = x(isfinite(x));
if numel(x) < 2
    s = 0;
    return
end
s = std(x, 0);
end


function r = column_rsq(x, y)
r = NaN;
x = x(:);
y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);
if numel(x) < 3 || column_std(x) <= 0 || column_std(y) <= 0
    return
end
beta = x \ y;
resid = y - x * beta;
ss = sum(y .^ 2);
if ss <= 0
    return
end
r = 1 - sum(resid .^ 2) / ss;
end


function r = matrix_rsq(X, y)
r = NaN;
y = y(:);
ok = isfinite(y) & all(isfinite(X), 2);
X = X(ok, :);
y = y(ok);
if size(X, 1) < size(X, 2) + 2 || column_std(y) <= 0
    return
end
spread = std(X, 0, 1);
if any(~isfinite(spread) | spread <= 0)
    return
end
beta = X \ y;
resid = y - X * beta;
ss = sum(y .^ 2);
if ss <= 0
    return
end
r = 1 - sum(resid .^ 2) / ss;
end


function bad = modified_z_rows(M)
bad = false(size(M, 1), 1);
for k = 1:size(M, 2)
    col = M(:, k);
    ok = isfinite(col);
    if nnz(ok) < 5
        continue
    end
    center = median(col(ok));
    mad = median(abs(col(ok) - center));
    if ~isfinite(mad) || mad <= 0
        continue
    end
    modified_z = 0.6745 * (col - center) / mad;
    bad = bad | (ok & abs(modified_z) > 3.5);
end
end


function rho = rank_corr(x, y)
rho = pearson_corr(average_rank(x), average_rank(y));
end


function r = pearson_corr(a, b)
a = a(:) - mean(a(:));
b = b(:) - mean(b(:));
den = sqrt(sum(a .^ 2) * sum(b .^ 2));
if ~isfinite(den) || den <= 0
    r = NaN;
    return
end
r = sum(a .* b) / den;
end


function ranks = average_rank(x)
[sorted, order] = sort(x(:));
ranks = zeros(numel(x), 1);
i = 1;
n = numel(sorted);
while i <= n
    j = i;
    while j < n && sorted(j + 1) == sorted(i)
        j = j + 1;
    end
    ranks(order(i:j)) = (i + j) / 2;
    i = j + 1;
end
end


function p = corr_p(r, df)
p = NaN;
if ~isfinite(r) || ~isfinite(df) || df < 1
    return
end
if abs(r) >= 1
    p = 0;
    return
end
tstat = r * sqrt(df / (1 - r ^ 2));
p = t_two_tailed(tstat, df);
end


function p = t_two_tailed(tstat, dof)
p = nan(size(tstat));
ok = isfinite(tstat) & isfinite(dof) & dof > 0;
ratio = dof ./ (dof + tstat(ok) .^ 2);
p(ok) = betainc(ratio, dof / 2, 0.5);
end


function p = one_sided_above(beta, p_two)
p = NaN;
if ~isfinite(beta) || ~isfinite(p_two)
    return
end
if beta >= 0
    p = p_two / 2;
else
    p = 1 - p_two / 2;
end
end


function p = signed_rank_p(x, tail)
p = NaN;
x = x(isfinite(x));
x = x(x ~= 0);
n = numel(x);
if n < 1
    return
end
ranks = average_rank(abs(x));
w_plus = sum(ranks(x > 0));
mu = n * (n + 1) / 4;
sigma2 = n * (n + 1) * (2 * n + 1) / 24 - tie_cubic(abs(x)) / 48;
if ~isfinite(sigma2) || sigma2 <= 0
    return
end
sigma = sqrt(sigma2);
z_hi = (w_plus - 0.5 - mu) / sigma;
z_lo = (w_plus + 0.5 - mu) / sigma;
p_right = min(max(0.5 * erfc(z_hi / sqrt(2)), 0), 1);
p_left = min(max(0.5 * erfc(-z_lo / sqrt(2)), 0), 1);
if strcmp(tail, 'right')
    p = p_right;
else
    p = min(1, 2 * min(p_right, p_left));
end
end


function t3 = tie_cubic(abs_x)
abs_x = sort(abs_x(:));
t3 = 0;
i = 1;
n = numel(abs_x);
while i <= n
    j = i;
    while j < n && abs_x(j + 1) == abs_x(i)
        j = j + 1;
    end
    t = j - i + 1;
    t3 = t3 + (t ^ 3 - t);
    i = j + 1;
end
end


function [n_pos, n_neg] = signed_counts(values, pvals, alpha)
sig = isfinite(values) & isfinite(pvals) & pvals < alpha;
n_pos = sum(sig & values > 0);
n_neg = sum(sig & values < 0);
end


function value = finite_median(values)
values = values(isfinite(values));
if isempty(values)
    value = NaN;
    return
end
value = median(values);
end


function labels = as_cellstr(values)
if iscell(values)
    labels = values(:);
    return
end
if ischar(values)
    labels = cellstr(values);
    return
end
error('locking_acg_within_neuron:labels', ...
    'Neuron and lock-group labels must be text.');
end


function tf = panel_has(group, panel)
if strcmp(panel, 'all')
    tf = true;
elseif strcmp(panel, 'delta')
    tf = any(strcmp(group, {'trough', 'peak'}));
else
    tf = strcmp(group, panel);
end
end
