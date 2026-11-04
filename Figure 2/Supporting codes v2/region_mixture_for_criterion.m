function [region_result, outcome] = region_mixture_for_criterion( ...
        region_result, criterion_name, frequency_limits, options)
%REGION_MIXTURE_FOR_CRITERION  Draw the mixture a stored rule selects.
%
% criterion_name is one of the names printed by
% Select_play_session_region_gaussian_centers.m. The component count is
% read from the saved measurements. Curves are rebuilt from
% component_parameters when that count was stored. A count that matches
% the mixture already saved is kept even when the parameter table is
% absent. outcome.applied is false when the rule cannot be drawn.

if nargin < 4 || isempty(options)
    options = struct();
end
criterion_name = strtrim(char(string(criterion_name)));
outcome = struct( ...
    'applied', false, ...
    'criterion', criterion_name, ...
    'selected_n', NaN, ...
    'note', '');
fit = region_result.fit;
if ~isstruct(fit) || ~isfield(fit, 'model_selection') || ...
        ~istable(fit.model_selection) || height(fit.model_selection) == 0
    outcome.note = 'no component-count measurements were stored';
    return
end
if ~isfield(options, 'n_observations') || isempty(options.n_observations)
    if isfield(fit, 'center_frequencies')
        options.n_observations = numel(fit.center_frequencies);
    end
end
selection = select_region_gaussian_criteria(fit.model_selection, options);
match = find(strcmpi(selection.Criterion, criterion_name), 1, 'first');
if isempty(match)
    outcome.note = sprintf('unknown selection criterion "%s"', ...
        criterion_name);
    return
end
outcome.criterion = char(selection.Criterion{match});
outcome.note = char(selection.Note{match});
chosen = selection.SelectedN(match);
outcome.selected_n = chosen;
if ~selection.Available(match) || ~isfinite(chosen)
    return
end

parameters = stored_component_parameters(fit);
parameter_rows = false(0, 1);
if ~isempty(parameters)
    parameter_rows = parameters.NComponents == chosen;
end
saved_count_matches = isfield(fit, 'best_n_components') && ...
    isequal(fit.best_n_components, chosen) && ...
    isfield(fit, 'mixture_density') && ~isempty(fit.mixture_density);
if any(parameter_rows)
    means = parameters.MeanHz(parameter_rows);
    sds = parameters.SDHz(parameter_rows);
    weights = parameters.Weight(parameter_rows);
    effective_n = parameters.EffectiveN(parameter_rows);
    hard_count = parameters.HardCount(parameter_rows);
    [component_table, density_frequency, component_density, ...
        mixture_density] = gaussian_mixture_density_curves( ...
        means, sds, weights, frequency_limits, effective_n, hard_count);
    fit.best_component_table = component_table;
    fit.density_frequency = density_frequency;
    fit.component_density = component_density;
    fit.mixture_density = mixture_density;
elseif ~saved_count_matches
    outcome.note = sprintf(['parameters for %g Gaussians were not ' ...
        'stored, so %s cannot be drawn'], chosen, outcome.criterion);
    return
end
fit.best_n_components = chosen;
fit.selection_criterion = outcome.criterion;
fit.status = 'fitted';
region_result.fit = fit;
outcome.applied = true;
end

function parameters = stored_component_parameters(fit)
parameters = [];
if ~isfield(fit, 'component_parameters') || ...
        ~istable(fit.component_parameters) || ...
        height(fit.component_parameters) == 0
    return
end
required = {'NComponents', 'MeanHz', 'SDHz', 'Weight'};
if ~all(ismember(required, ...
        fit.component_parameters.Properties.VariableNames))
    return
end
parameters = fit.component_parameters;
if ~ismember('EffectiveN', parameters.Properties.VariableNames)
    parameters.EffectiveN = nan(height(parameters), 1);
end
if ~ismember('HardCount', parameters.Properties.VariableNames)
    parameters.HardCount = nan(height(parameters), 1);
end
end
