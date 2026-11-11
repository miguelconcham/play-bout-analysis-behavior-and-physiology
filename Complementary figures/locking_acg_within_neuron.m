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
clear diary_cleanup
if options.do_plot
    plot_slope_distributions(result, options, 'slope', 'Within-neuron slopes');
    plot_slope_distributions(result, options, 'partial_slope', ...
        'Within-neuron slopes after firing rate');
    plot_common_slopes(result, options);
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
    'P', {}, 'y', {}, 'rate', {});
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
per_neuron.n_other_group = 0;
end


function cells = empty_cells(n_panel)
cells = cell(1, n_panel);
for p = 1:n_panel
    cells{p} = zeros(0, 1);
end
end


function per_neuron = append_record(per_neuron, p, rec)
if isfinite(rec.rho)
    per_neuron.rho{p}(end + 1, 1) = rec.rho;
    per_neuron.rho_p{p}(end + 1, 1) = rec.rho_p;
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
disp(sprintf(['      %s: n = %d, median = %.3g, two-sided signrank p = %.3g, ' ...
    'one-sided positive drop p = %.3g'], ...
    label, nnz(isfinite(values)), finite_median(values), ...
    signed_rank_p(values, 'both'), signed_rank_p(values, 'right')))
end


function plot_slope_distributions(result, options, field_name, figure_name)
n_out = numel(result.outcome_titles);
n_panel = numel(result.panel_titles);
figure('Name', figure_name, 'Color', 'w', 'Visible', options.figure_visible)
set(gcf, 'Units', 'pixels', 'Position', [40 40 1680 980])
short_panels = {'All', 'Trough', 'Peak', 'Delta', 'Non-entrained'};
for o = 1:n_out
    edges = shared_edges(result.per_neuron(o).(field_name));
    for p = 1:n_panel
        subplot(n_out, n_panel, (o - 1) * n_panel + p)
        values = result.per_neuron(o).(field_name){p};
        draw_hist(values, edges, [0.25 0.38 0.62])
        hold on
        mark_zero()
        mark_median(finite_median(values))
        title(sprintf('%s n=%d', short_panels{p}, nnz(isfinite(values))))
        if o == n_out
            xlabel('Slope')
        end
        if p == 1
            ylabel(result.outcome_titles{o})
        end
        box off
    end
end
save_figure(options.figure_folder, figure_name);
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


function draw_hist(values, edges, color)
values = values(isfinite(values));
n_bin = numel(edges) - 1;
if n_bin < 1
    return
end
counts = zeros(1, n_bin);
for k = 1:n_bin
    if k < n_bin
        counts(k) = sum(values >= edges(k) & values < edges(k + 1));
    else
        counts(k) = sum(values >= edges(k) & values <= edges(k + 1));
    end
end
centers = (edges(1:end - 1) + edges(2:end)) / 2;
h = bar(centers, counts, 1);
set(h, 'FaceColor', color, 'EdgeColor', 'none')
xlim(edges([1 end]))
end


function mark_zero()
yl = ylim;
plot([0 0], yl, '--', 'Color', [0.35 0.35 0.35], 'HandleVisibility', 'off')
ylim(yl)
end


function mark_median(value)
if ~isfinite(value)
    return
end
yl = ylim;
plot([value value], yl, '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 1.2, ...
    'HandleVisibility', 'off')
ylim(yl)
end


function edges = shared_edges(panels)
values = zeros(0, 1);
for p = 1:numel(panels)
    values = [values; panels{p}(:)]; %#ok<AGROW>
end
values = values(isfinite(values));
if numel(values) < 2
    edges = linspace(-1, 1, 21);
    return
end
lo = min([values; 0]);
hi = max([values; 0]);
pad = 0.05 * (hi - lo);
if pad <= 0
    pad = 1;
end
edges = linspace(lo - pad, hi + pad, 21);
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
