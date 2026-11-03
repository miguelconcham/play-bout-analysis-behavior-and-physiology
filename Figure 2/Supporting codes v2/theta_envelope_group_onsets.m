function grouped = theta_envelope_group_onsets(recordings, min_events)
%THETA_ENVELOPE_GROUP_ONSETS  Animal-level onset locking from recording summaries.
%
% Within an animal, recordings are combined with weights equal to the number
% of valid pre-onset windows. The group null is the mean across animals of
% one onset-shift draw from each animal. Animals with fewer than min_events
% valid windows stay in the output and are left out of the group p-values.
%
% The primary window is the first pre window, -2 to 0 s. Phase fraction and
% leave-one-out correlation were computed inside each recording.

grouped = empty_grouped();
if isempty(recordings)
    return
end
if nargin < 2 || isempty(min_events)
    min_events = 15;
end
names = {recordings.animal_area_table_name};
animal_names = unique(names, 'stable');
n_animals = numel(animal_names);
frequencies = recordings(1).frequencies_test;
mask = recordings(1).confirmatory_mask;
n_freq = numel(frequencies);
grouped.animals = animal_names(:);
grouped.frequencies = frequencies;
grouped.confirmatory_mask = mask;
grouped.epoch_time = recordings(1).onset.epoch_time;
grouped.n_trials = zeros(n_animals, 1);
grouped.included = false(n_animals, 1);
grouped.locked_amplitude = nan(n_animals, n_freq);
grouped.ppc = nan(n_animals, n_freq);
grouped.loo_mean_r = nan(n_animals, 1);
grouped.phase_fraction = nan(n_animals, 1);
grouped.mean_epoch = nan(n_animals, numel(grouped.epoch_time));
grouped.clean_mean_epoch = nan(n_animals, numel(grouped.epoch_time));
grouped.confirmatory_locked = nan(n_animals, 1);
grouped.confirmatory_ppc = nan(n_animals, 1);
grouped.animal_locked_p = nan(n_animals, 1);
grouped.animal_ppc_p = nan(n_animals, 1);
grouped.animal_loo_p = nan(n_animals, 1);
grouped.animal_phase_p = nan(n_animals, 1);

null_locked = cell(n_animals, 1);
null_ppc = cell(n_animals, 1);
null_loo = cell(n_animals, 1);
null_phase = cell(n_animals, 1);
for animal_n = 1:n_animals
    rows = recordings(strcmp(names, animal_names{animal_n}));
    [summary, nulls] = one_animal(rows, frequencies, mask);
    grouped.n_trials(animal_n) = summary.n_trials;
    grouped.included(animal_n) = summary.n_trials >= min_events;
    grouped.locked_amplitude(animal_n, :) = summary.locked_amplitude;
    grouped.ppc(animal_n, :) = summary.ppc;
    grouped.loo_mean_r(animal_n) = summary.loo_mean_r;
    grouped.phase_fraction(animal_n) = summary.phase_fraction;
    grouped.mean_epoch(animal_n, :) = summary.mean_epoch;
    grouped.clean_mean_epoch(animal_n, :) = summary.clean_mean_epoch;
    grouped.confirmatory_locked(animal_n) = summary.confirmatory_locked;
    grouped.confirmatory_ppc(animal_n) = summary.confirmatory_ppc;
    grouped.animal_locked_p(animal_n) = tail_p( ...
        summary.confirmatory_locked, nulls.confirmatory_locked);
    grouped.animal_ppc_p(animal_n) = tail_p( ...
        summary.confirmatory_ppc, nulls.confirmatory_ppc);
    grouped.animal_loo_p(animal_n) = tail_p(summary.loo_mean_r, nulls.loo);
    grouped.animal_phase_p(animal_n) = tail_p(summary.phase_fraction, nulls.phase);
    null_locked{animal_n} = nulls.confirmatory_locked;
    null_ppc{animal_n} = nulls.confirmatory_ppc;
    null_loo{animal_n} = nulls.loo;
    null_phase{animal_n} = nulls.phase;
end
included = grouped.included;
grouped.locked_p = group_p(grouped.confirmatory_locked, null_locked, included);
grouped.ppc_p = group_p(grouped.confirmatory_ppc, null_ppc, included);
grouped.loo_p = group_p(grouped.loo_mean_r, null_loo, included);
grouped.phase_fraction_p = group_p(grouped.phase_fraction, null_phase, included);
if any(included)
    grouped.grand_mean_epoch = finite_mean(grouped.mean_epoch(included, :));
    grouped.grand_clean_epoch = finite_mean(grouped.clean_mean_epoch(included, :));
else
    grouped.grand_mean_epoch = nan(1, numel(grouped.epoch_time));
    grouped.grand_clean_epoch = nan(1, numel(grouped.epoch_time));
end
end

function [summary, nulls] = one_animal(rows, frequencies, mask)
n_freq = numel(frequencies);
n_recordings = numel(rows);
weights = zeros(n_recordings, 1);
pair_weights = zeros(n_recordings, 1);
observed = nan(n_recordings, n_freq);
ppc = nan(n_recordings, n_freq);
loo = nan(n_recordings, 1);
phase = nan(n_recordings, 1);
n_perm = inf;
for row_n = 1:n_recordings
    pre = rows(row_n).onset.pre(1);
    weights(row_n) = pre.n_valid;
    pair_weights(row_n) = pre.n_valid * max(pre.n_valid - 1, 0) / 2;
    observed(row_n, :) = pre.mean_coefficient(:).';
    ppc(row_n, :) = pre.ppc(:).';
    loo(row_n) = pre.loo_mean_r;
    phase(row_n) = pre.phase_fraction;
    n_perm = min(n_perm, size(pre.null_mean_coefficient, 1));
end
if isinf(n_perm)
    n_perm = 0;
end
summary.n_trials = sum(weights);
summary.locked_amplitude = abs(weighted_mean(observed, weights));
summary.ppc = weighted_mean(ppc, pair_weights);
summary.loo_mean_r = weighted_mean(loo, weights);
summary.phase_fraction = weighted_mean(phase, weights);
summary.mean_epoch = pooled_epoch(rows, false);
summary.clean_mean_epoch = pooled_epoch(rows, true);
summary.confirmatory_locked = mean_mask(summary.locked_amplitude, mask);
summary.confirmatory_ppc = mean_mask(summary.ppc, mask);

nulls.confirmatory_locked = nan(n_perm, 1);
nulls.confirmatory_ppc = nan(n_perm, 1);
nulls.loo = nan(n_perm, 1);
nulls.phase = nan(n_perm, 1);
for perm = 1:n_perm
    null_coefficient = nan(n_recordings, n_freq);
    null_ppc_rows = nan(n_recordings, n_freq);
    null_loo_rows = nan(n_recordings, 1);
    null_phase_rows = nan(n_recordings, 1);
    for row_n = 1:n_recordings
        pre = rows(row_n).onset.pre(1);
        null_coefficient(row_n, :) = pre.null_mean_coefficient(perm, :);
        null_ppc_rows(row_n, :) = pre.null_ppc(perm, :);
        if ~isempty(pre.null_loo_mean_r)
            null_loo_rows(row_n) = pre.null_loo_mean_r(perm);
        end
        if ~isempty(pre.null_phase_fraction)
            null_phase_rows(row_n) = pre.null_phase_fraction(perm);
        end
    end
    nulls.confirmatory_locked(perm) = mean_mask( ...
        abs(weighted_mean(null_coefficient, weights)), mask);
    nulls.confirmatory_ppc(perm) = mean_mask( ...
        weighted_mean(null_ppc_rows, pair_weights), mask);
    nulls.loo(perm) = weighted_mean(null_loo_rows, weights);
    nulls.phase(perm) = weighted_mean(null_phase_rows, weights);
end
end

function value = pooled_epoch(rows, clean_only)
value = nan(1, numel(rows(1).onset.epoch_time));
stack = [];
for row_n = 1:numel(rows)
    epochs = rows(row_n).onset.epochs;
    valid = all(isfinite(epochs), 2);
    if clean_only
        valid = valid & rows(row_n).onset.clean_pre(:);
    end
    if any(valid)
        stack = [stack; epochs(valid, :)]; %#ok<AGROW>
    end
end
if ~isempty(stack)
    value = mean(stack, 1);
end
end

function value = weighted_mean(values, weights)
values = double(values);
weights = double(weights(:));
if isvector(values) && numel(values) == numel(weights) && size(values, 2) ~= 1
    values = values(:);
end
if size(values, 1) ~= numel(weights)
    error('theta_envelope_group_onsets:weights', ...
        'Weights must have one entry per recording row.');
end
n_columns = size(values, 2);
value = nan(1, n_columns);
for column = 1:n_columns
    finite = isfinite(values(:, column)) & isfinite(weights) & weights > 0;
    if any(finite)
        value(column) = sum(weights(finite) .* values(finite, column)) / sum(weights(finite));
    end
end
if isvector(values) && n_columns == 1
    value = value(1);
end
end

function value = finite_mean(matrix)
value = nan(1, size(matrix, 2));
for column = 1:size(matrix, 2)
    finite = matrix(isfinite(matrix(:, column)), column);
    if ~isempty(finite)
        value(column) = mean(finite);
    end
end
end

function value = mean_mask(values, mask)
selected = values(logical(mask));
selected = selected(isfinite(selected));
if isempty(selected)
    value = NaN;
else
    value = mean(selected);
end
end

function p = tail_p(observed, null_draws)
draws = null_draws(isfinite(null_draws));
if isempty(draws) || ~isfinite(observed)
    p = NaN;
else
    p = mean(draws >= observed);
end
end

function p = group_p(observed, null_cells, included)
observed = observed(included);
null_cells = null_cells(included);
if isempty(observed) || all(~isfinite(observed))
    p = NaN;
    return
end
n_perm = numel(null_cells{1});
null_mean = nan(n_perm, 1);
for perm = 1:n_perm
    draw = nan(numel(observed), 1);
    for animal_n = 1:numel(observed)
        if numel(null_cells{animal_n}) >= perm
            draw(animal_n) = null_cells{animal_n}(perm);
        end
    end
    finite = isfinite(draw) & isfinite(observed);
    if any(finite)
        null_mean(perm) = mean(draw(finite));
    end
end
finite_observed = observed(isfinite(observed));
finite_null = null_mean(isfinite(null_mean));
if isempty(finite_null)
    p = NaN;
else
    p = mean(finite_null >= mean(finite_observed));
end
end

function grouped = empty_grouped()
grouped = struct( ...
    'animals', {{}}, ...
    'frequencies', [], ...
    'confirmatory_mask', [], ...
    'epoch_time', [], ...
    'n_trials', [], ...
    'included', [], ...
    'locked_amplitude', [], ...
    'ppc', [], ...
    'loo_mean_r', [], ...
    'phase_fraction', [], ...
    'mean_epoch', [], ...
    'clean_mean_epoch', [], ...
    'grand_mean_epoch', [], ...
    'grand_clean_epoch', [], ...
    'confirmatory_locked', [], ...
    'confirmatory_ppc', [], ...
    'animal_locked_p', [], ...
    'animal_ppc_p', [], ...
    'animal_loo_p', [], ...
    'animal_phase_p', [], ...
    'locked_p', NaN, ...
    'ppc_p', NaN, ...
    'loo_p', NaN, ...
    'phase_fraction_p', NaN);
end
