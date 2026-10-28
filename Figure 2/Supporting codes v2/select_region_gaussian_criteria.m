function selection = select_region_gaussian_criteria(model_selection, options)
%SELECT_REGION_GAUSSIAN_CRITERIA  Compare stored mixture-selection rules.
%
% model_selection is one region's per-count table from
% fit_gaussian_center_mixture. Rules that need only AIC, BIC, and the
% negative log-likelihood run on an older save. Entropy, overlap, and
% goodness-of-fit rules run when those columns are present. Lower is
% better for the information criteria. A tie keeps the smaller count.

if nargin < 2 || isempty(options)
    options = struct();
end
options = fill_selection_option(options, 'elbow_fraction', 0.05);
options = fill_selection_option(options, 'ashman_d_minimum', 2);
options = fill_selection_option(options, 'bhattacharyya_maximum', 0.5);
options = fill_selection_option(options, 'gof_relative_tolerance', 0.1);
options = fill_selection_option(options, 'nec_one_component_threshold', 1);
options = fill_selection_option(options, 'n_observations', NaN);

criterion = cell(0, 1);
selected_n = zeros(0, 1);
criterion_value = zeros(0, 1);
available = false(0, 1);
note = cell(0, 1);

if isempty(model_selection) || ~istable(model_selection) || ...
        height(model_selection) == 0 || ...
        ~ismember('NComponents', model_selection.Properties.VariableNames)
    selection = selection_table(criterion, selected_n, ...
        criterion_value, available, note);
    return
end

k = double(model_selection.NComponents(:));
[has_aic, aic] = numeric_column(model_selection, 'AIC');
[has_bic, bic] = numeric_column(model_selection, 'BIC');
[has_nll, nll] = numeric_column(model_selection, 'NegativeLogLikelihood');
[has_entropy, entropy] = numeric_column(model_selection, 'Entropy');
[has_parameters, n_parameters] = numeric_column( ...
    model_selection, 'NParameters');
[has_observations, n_observations] = numeric_column( ...
    model_selection, 'NObservations');
[has_completed, completed] = numeric_column( ...
    model_selection, 'CompletedLogLikelihood');
[has_hard, hard_completed] = numeric_column( ...
    model_selection, 'HardCompletedLogLikelihood');
[has_ashman, ashman] = numeric_column( ...
    model_selection, 'MinimumNeighborAshmanD');
[has_bhattacharyya, bhattacharyya] = numeric_column( ...
    model_selection, 'MaximumPairBhattacharyya');
[has_ks, ks_distance] = numeric_column(model_selection, 'KSDistance');
[has_ad, anderson] = numeric_column(model_selection, 'AndersonDarling');
[has_chi, chi_square] = numeric_column(model_selection, 'ChiSquare');
[has_l1, histogram_l1] = numeric_column(model_selection, 'HistogramL1');

if ismember('Eligible', model_selection.Properties.VariableNames)
    eligible = logical(model_selection.Eligible(:));
else
    eligible = isfinite(nll);
end
if ~has_observations && isfinite(options.n_observations)
    n_observations = repmat(double(options.n_observations), numel(k), 1);
    has_observations = true;
end
if ~has_parameters && has_aic && has_nll
    n_parameters = (aic - 2 * nll) / 2;
    close_to_integer = abs(n_parameters - round(n_parameters)) < 1e-4;
    n_parameters(close_to_integer) = round(n_parameters(close_to_integer));
    has_parameters = true;
end
if ~has_completed && has_nll && has_entropy
    completed = -nll - entropy;
    has_completed = true;
end

remember('AIC', aic, has_aic, ...
    'minimum AIC; a tie keeps the smaller count');
remember('AICc', corrected_aic(aic, n_parameters, n_observations), ...
    has_aic && has_parameters && has_observations, ...
    'AIC with the small-sample penalty');
remember('BIC', bic, has_bic, ...
    'minimum BIC; a tie keeps the smaller count');
remember('HannanQuinn', hannan_quinn(nll, n_parameters, n_observations), ...
    has_nll && has_parameters && has_observations, ...
    'Hannan-Quinn penalty');
remember('SampleSizeAdjustedBIC', sample_size_bic( ...
    nll, n_parameters, n_observations), ...
    has_nll && has_parameters && has_observations, ...
    'BIC with log((n+2)/24)');
remember('ICL', integrated_completed_likelihood(bic, entropy), ...
    has_bic && has_entropy, ...
    'BIC plus twice the soft-assignment entropy');
remember('CLC', classification_criterion(completed), has_completed, ...
    'minus twice the soft completed log-likelihood');
remember('HardCLC', classification_criterion(hard_completed), has_hard, ...
    'minus twice the hard-assignment completed log-likelihood');
remember('AWE', approximate_weight_of_evidence( ...
    completed, n_parameters, n_observations), ...
    has_completed && has_parameters && has_observations, ...
    'completed likelihood with the AWE penalty');
remember_nec(k, nll, entropy, eligible, has_nll && has_entropy, ...
    options.nec_one_component_threshold);
remember_slope('SlopeHeuristicJump', k, nll, n_parameters, eligible, ...
    has_nll && has_parameters, 'jump');
remember_slope('SlopeHeuristicRegression', k, nll, n_parameters, ...
    eligible, has_nll && has_parameters, 'regression');
remember_elbow('ElbowBICSuccessive', k, bic, eligible, has_bic, ...
    options.elbow_fraction, 'successive');
remember_elbow('ElbowBICCumulative', k, bic, eligible, has_bic, ...
    options.elbow_fraction, 'cumulative');
remember_elbow('ElbowLogLikelihoodSuccessive', k, nll, eligible, ...
    has_nll, options.elbow_fraction, 'successive');
remember_elbow('ElbowLogLikelihoodCumulative', k, nll, eligible, ...
    has_nll, options.elbow_fraction, 'cumulative');
remember_elbow('ElbowChiSquareSuccessive', k, chi_square, eligible, ...
    has_chi, options.elbow_fraction, 'successive');
remember_within('ChiSquareWithinTolerance', k, chi_square, eligible, ...
    has_chi, options.gof_relative_tolerance, ...
    'smallest count within tolerance of the best chi-square');
remember_elbow('ElbowKSSuccessive', k, ks_distance, eligible, ...
    has_ks, options.elbow_fraction, 'successive');
remember_within('KSDistanceWithinTolerance', k, ks_distance, eligible, ...
    has_ks, options.gof_relative_tolerance, ...
    'smallest count within tolerance of the best KS distance');
remember_elbow('ElbowAndersonDarlingSuccessive', k, anderson, eligible, ...
    has_ad, options.elbow_fraction, 'successive');
remember_within('AndersonDarlingWithinTolerance', k, anderson, eligible, ...
    has_ad, options.gof_relative_tolerance, ...
    'smallest count within tolerance of the best Anderson-Darling');
remember_elbow('ElbowHistogramL1Successive', k, histogram_l1, eligible, ...
    has_l1, options.elbow_fraction, 'successive');
remember_within('HistogramL1WithinTolerance', k, histogram_l1, eligible, ...
    has_l1, options.gof_relative_tolerance, ...
    'smallest count within tolerance of the best histogram L1');
remember_ashman(k, ashman, eligible, has_ashman, options.ashman_d_minimum);
remember_overlap(k, bhattacharyya, eligible, has_bhattacharyya, ...
    options.bhattacharyya_maximum);

selection = selection_table(criterion, selected_n, ...
    criterion_value, available, note);

    function remember(name, values, is_present, note_text)
        if ~is_present
            push(name, NaN, NaN, false, ['unavailable: ' note_text]);
            return
        end
        [best_k, best_value] = minimize_eligible(k, values, eligible);
        if ~isfinite(best_k)
            push(name, NaN, NaN, true, 'no eligible finite value');
            return
        end
        push(name, best_k, best_value, true, note_text);
    end

    function remember_elbow(name, counts, values, keep, is_present, ...
            fraction, kind)
        if ~is_present
            push(name, NaN, NaN, false, ...
                'unavailable: the source column was not stored');
            return
        end
        if strcmp(kind, 'cumulative')
            best_k = elbow_cumulative(counts, values, keep, fraction);
            note_text = sprintf(['first count that reaches %.0f%% of ' ...
                'the total drop'], 100 * (1 - fraction));
        else
            best_k = elbow_successive(counts, values, keep, fraction);
            note_text = sprintf(['stop when the next gain is below ' ...
                '%.0f%% of the first gain'], 100 * fraction);
        end
        push(name, best_k, value_at(counts, values, best_k), ...
            true, note_text);
    end

    function remember_within(name, counts, values, keep, is_present, ...
            fraction, note_text)
        if ~is_present
            push(name, NaN, NaN, false, ...
                'unavailable: the source column was not stored');
            return
        end
        best_k = smallest_within_tolerance(counts, values, keep, fraction);
        push(name, best_k, value_at(counts, values, best_k), ...
            true, note_text);
    end

    function remember_slope(name, counts, values, dimensions, keep, ...
            is_present, kind)
        if ~is_present
            push(name, NaN, NaN, false, ...
                'unavailable: likelihood or parameter count was not stored');
            return
        end
        use = keep & isfinite(values) & isfinite(dimensions) & ...
            isfinite(counts);
        if nnz(use) < 2
            push(name, NaN, NaN, true, ...
                'fewer than two eligible models');
            return
        end
        if strcmp(kind, 'regression')
            [best_k, kappa_min] = slope_regression( ...
                counts(use), values(use), dimensions(use));
            note_text = sprintf(['twice the tail slope; kappa_min ' ...
                '= %.3g'], kappa_min);
        else
            [best_k, kappa_min] = slope_dimension_jump( ...
                counts(use), values(use), dimensions(use));
            note_text = sprintf(['twice the penalty at the largest ' ...
                'dimension jump; kappa_min = %.3g'], kappa_min);
        end
        if ~isfinite(best_k)
            push(name, NaN, NaN, true, note_text);
            return
        end
        score = values + (2 * kappa_min) * dimensions;
        push(name, best_k, value_at(counts, score, best_k), true, note_text);
    end

    function remember_nec(counts, values, entropies, keep, is_present, ...
            threshold)
        if ~is_present
            push('NEC', NaN, NaN, false, ...
                'unavailable: entropy was not stored');
            return
        end
        [best_k, best_value, used_one] = normalized_entropy( ...
            counts, values, entropies, keep, threshold);
        if used_one
            note_text = sprintf(['minimum NEC is above %.3g, so the ' ...
                'one-component model is kept'], threshold);
        else
            note_text = 'minimum normalized entropy for two or more components';
        end
        push('NEC', best_k, best_value, true, note_text);
    end

    function remember_ashman(counts, separation, keep, is_present, minimum_d)
        if ~is_present
            push('AshmanSeparation', NaN, NaN, false, ...
                'unavailable: neighbor separation was not stored');
            return
        end
        pass = keep & (counts == 1 | ...
            (isfinite(separation) & separation >= minimum_d));
        best_k = largest_count(counts, pass);
        push('AshmanSeparation', best_k, ...
            value_at(counts, separation, best_k), true, ...
            sprintf('largest count whose neighbor Ashman D is at least %.3g', ...
            minimum_d));
    end

    function remember_overlap(counts, overlap, keep, is_present, maximum_overlap)
        if ~is_present
            push('BhattacharyyaOverlap', NaN, NaN, false, ...
                'unavailable: component overlap was not stored');
            return
        end
        pass = keep & (counts == 1 | ...
            (isfinite(overlap) & overlap <= maximum_overlap));
        best_k = largest_count(counts, pass);
        push('BhattacharyyaOverlap', best_k, ...
            value_at(counts, overlap, best_k), true, ...
            sprintf(['largest count whose maximum Bhattacharyya ' ...
            'coefficient is at most %.3g'], maximum_overlap));
    end

    function push(name, best_k, best_value, is_available, note_text)
        criterion{end + 1, 1} = name; %#ok<AGROW>
        selected_n(end + 1, 1) = best_k; %#ok<AGROW>
        criterion_value(end + 1, 1) = best_value; %#ok<AGROW>
        available(end + 1, 1) = logical(is_available); %#ok<AGROW>
        note{end + 1, 1} = note_text; %#ok<AGROW>
    end
end

function selection = selection_table(criterion, selected_n, ...
        criterion_value, available, note)
selection = table(criterion, selected_n, criterion_value, available, note, ...
    'VariableNames', {'Criterion', 'SelectedN', 'CriterionValue', ...
    'Available', 'Note'});
end

function [present, values] = numeric_column(model_selection, name)
present = ismember(name, model_selection.Properties.VariableNames);
if present
    values = double(model_selection.(name)(:));
else
    values = nan(height(model_selection), 1);
end
end

function options = fill_selection_option(options, name, default_value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = default_value;
end
end

function [best_k, best_value] = minimize_eligible(k, value, eligible)
keep = eligible(:) & isfinite(k(:)) & isfinite(value(:));
if ~any(keep)
    best_k = NaN;
    best_value = NaN;
    return
end
[best_value, index] = min(value(keep));
kept_k = k(keep);
best_k = kept_k(index);
end

function value = value_at(k, values, selected)
value = NaN;
if ~isfinite(selected)
    return
end
match = find(k(:)== selected, 1, 'first');
if isempty(match) || ~isfinite(values(match))
    return
end
value = values(match);
end

function best_k = largest_count(k, pass)
pass = pass(:) & isfinite(k(:));
if ~any(pass)
    best_k = NaN;
    return
end
best_k = max(k(pass));
end

function values = corrected_aic(aic, n_parameters, n_observations)
values = aic + 2 .* n_parameters .* (n_parameters + 1) ./ ...
    (n_observations - n_parameters - 1);
invalid = ~isfinite(aic) | ~isfinite(n_parameters) | ...
    ~isfinite(n_observations) | n_observations <= n_parameters + 1;
values(invalid) = NaN;
end

function values = hannan_quinn(nll, n_parameters, n_observations)
values = 2 .* nll + 2 .* n_parameters .* log(log(n_observations));
invalid = ~isfinite(nll) | ~isfinite(n_parameters) | ...
    ~isfinite(n_observations) | n_observations <= exp(1);
values(invalid) = NaN;
end

function values = sample_size_bic(nll, n_parameters, n_observations)
values = 2 .* nll + n_parameters .* log((n_observations + 2) ./ 24);
invalid = ~isfinite(nll) | ~isfinite(n_parameters) | ...
    ~isfinite(n_observations) | n_observations <= 0;
values(invalid) = NaN;
end

function values = integrated_completed_likelihood(bic, entropy)
values = bic + 2 .* entropy;
values(~isfinite(bic) | ~isfinite(entropy)) = NaN;
end

function values = classification_criterion(completed_log_likelihood)
values = -2 .* completed_log_likelihood;
values(~isfinite(completed_log_likelihood)) = NaN;
end

function values = approximate_weight_of_evidence( ...
        completed_log_likelihood, n_parameters, n_observations)
values = -2 .* completed_log_likelihood + 2 .* n_parameters .* ...
    (1.5 + log(n_observations));
invalid = ~isfinite(completed_log_likelihood) | ...
    ~isfinite(n_parameters) | ~isfinite(n_observations) | ...
    n_observations <= 0;
values(invalid) = NaN;
end

function [best_k, best_value, used_one] = normalized_entropy( ...
        k, nll, entropy, eligible, threshold)
used_one = false;
best_k = NaN;
best_value = NaN;
one = eligible(:) & k(:) == 1 & isfinite(nll(:));
if ~any(one)
    return
end
nll_one = nll(find(one, 1, 'first'));
keep = eligible(:) & k(:) >= 2 & isfinite(entropy(:)) & isfinite(nll(:));
if ~any(keep)
    return
end
gain = -nll(keep) + nll_one;
nec = nan(nnz(keep), 1);
positive = gain > 0;
kept_entropy = entropy(keep);
nec(positive) = kept_entropy(positive) ./ gain(positive);
finite_nec = isfinite(nec);
if ~any(finite_nec)
    return
end
[best_value, finite_index] = min(nec(finite_nec));
finite_positions = find(finite_nec);
index = finite_positions(finite_index);
kept_k = k(keep);
if best_value > threshold
    best_k = 1;
    used_one = true;
    return
end
best_k = kept_k(index);
end

function best_k = elbow_successive(k, value, eligible, fraction)
[kk, vv] = eligible_series(k, value, eligible);
if isempty(kk)
    best_k = NaN;
    return
end
if numel(kk) < 2
    best_k = kk(1);
    return
end
first_gain = vv(1) - vv(2);
if first_gain <= 0
    best_k = kk(1);
    return
end
best_k = kk(1);
for index = 2:numel(kk)
    gain = vv(index - 1) - vv(index);
    if gain < fraction * first_gain
        return
    end
    best_k = kk(index);
end
end

function best_k = elbow_cumulative(k, value, eligible, fraction)
[kk, vv] = eligible_series(k, value, eligible);
if isempty(kk)
    best_k = NaN;
    return
end
total = vv(1) - min(vv);
if total <= 0
    best_k = kk(1);
    return
end
achieved = vv(1) - vv;
index = find(achieved >= (1 - fraction) * total, 1, 'first');
best_k = kk(index);
end

function best_k = smallest_within_tolerance(k, value, eligible, fraction)
keep = eligible(:) & isfinite(k(:)) & isfinite(value(:));
if ~any(keep)
    best_k = NaN;
    return
end
best = min(value(keep));
if best > 0
    limit = best * (1 + fraction);
else
    limit = best + abs(fraction);
end
pass = keep & value(:) <= limit + 1e-12;
best_k = min(k(pass));
end

function [kk, vv] = eligible_series(k, value, eligible)
keep = eligible(:) & isfinite(k(:)) & isfinite(value(:));
kk = k(keep);
vv = value(keep);
[kk, order] = sort(kk(:));
vv = vv(order);
end

function [best_k, kappa_min] = slope_dimension_jump(k, nll, dim)
[k, order] = sort(k(:));
nll = nll(order);
dim = dim(order);
best_k = NaN;
kappa_min = NaN;
pair_kappa = [];
for first_n = 1:numel(k) - 1
    for second_n = first_n + 1:numel(k)
        delta_dim = dim(second_n) - dim(first_n);
        delta_nll = nll(first_n) - nll(second_n);
        if delta_dim > 0 && delta_nll > 0
            pair_kappa(end + 1) = delta_nll / delta_dim; %#ok<AGROW>
        end
    end
end
if isempty(pair_kappa)
    [~, index] = min(nll);
    best_k = k(index);
    kappa_min = 0;
    return
end
kappa_max = max(pair_kappa) * 1.01;
grid = linspace(0, kappa_max, 2000);
selected_dim = zeros(numel(grid), 1);
selected_k = zeros(numel(grid), 1);
for grid_n = 1:numel(grid)
    score = nll + grid(grid_n) * dim;
    [~, index] = min(score);
    selected_dim(grid_n) = dim(index);
    selected_k(grid_n) = k(index);
end
jump = selected_dim(1:end - 1) - selected_dim(2:end);
[jump_size, jump_at] = max(jump);
if jump_size <= 0
    best_k = selected_k(end);
    kappa_min = grid(end);
    return
end
kappa_min = grid(jump_at + 1);
score = nll + (2 * kappa_min) * dim;
[~, index] = min(score);
best_k = k(index);
end

function [best_k, kappa_min] = slope_regression(k, nll, dim)
[k, order] = sort(k(:));
nll = nll(order);
dim = dim(order);
best_k = NaN;
kappa_min = NaN;
cutoff = median(dim);
use = dim >= cutoff;
if nnz(use) < 2
    use(:) = true;
end
if nnz(use) < 2
    return
end
design = [ones(nnz(use), 1), dim(use)];
beta = design \ nll(use);
slope = beta(2);
if ~isfinite(slope) || slope >= 0
    return
end
kappa_min = -slope;
score = nll + (2 * kappa_min) * dim;
[~, index] = min(score);
best_k = k(index);
end
