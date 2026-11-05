function plan = gaussian_mixture_search_plan(n_centers, n_unique, options)
%GAUSSIAN_MIXTURE_SEARCH_PLAN  Component cap, restarts, and workload.
%
% max_components_tested never exceeds options.max_components. The default
% cap is 20. Restarts stay at options.replicates unless
% large_sample_replicates is smaller and n_centers exceeds
% large_sample_centers. work_units counts one unit per center, per
% restart, per component count from 1 through the cap.
% compared_with_reference is that cost divided by one mid-LPAG-sized
% search: 300 centers, the full restart count, and the same cap.

if nargin < 3 || isempty(options)
    options = struct();
end
options = gaussian_mixture_options(options);
n_centers = double(n_centers);
if nargin < 2 || isempty(n_unique) || ~isfinite(n_unique)
    n_unique = n_centers;
end
n_unique = double(n_unique);

max_components = min([options.max_components, n_centers, n_unique, ...
    floor(n_centers / options.minimum_observations_per_component)]);
if ~isfinite(max_components) || max_components < 1
    max_components = 0;
end

replicates = options.replicates;
restarts_reduced = false;
if n_centers > options.large_sample_centers && ...
        options.large_sample_replicates < options.replicates
    replicates = options.large_sample_replicates;
    restarts_reduced = true;
end

if max_components < 1
    work_units = 0;
else
    work_units = n_centers * replicates * sum(1:max_components);
end
reference_work_units = reference_search_work(options);
if reference_work_units > 0
    compared_with_reference = work_units / reference_work_units;
else
    compared_with_reference = NaN;
end

if max_components < 1
    runtime_note = 'too few centers to fit';
elseif restarts_reduced && compared_with_reference > 10
    runtime_note = ['fewer restarts than requested; the first ' ...
        'component fit will project the clock time'];
elseif restarts_reduced
    runtime_note = 'fewer restarts than requested for this sample size';
elseif compared_with_reference > 10
    runtime_note = ['full restart search through the component cap; ' ...
        'the first component fit will project the clock time'];
else
    runtime_note = 'within the size of the original mid-LPAG search';
end

plan = struct( ...
    'n_centers', n_centers, ...
    'n_unique', n_unique, ...
    'max_components_cap', options.max_components, ...
    'max_components_tested', max_components, ...
    'replicates_requested', options.replicates, ...
    'replicates', replicates, ...
    'restarts_reduced', restarts_reduced, ...
    'work_units', work_units, ...
    'reference_work_units', reference_work_units, ...
    'compared_with_reference', compared_with_reference, ...
    'runtime_note', runtime_note);
end

function work_units = reference_search_work(options)
% One mid-LPAG-sized search: a few hundred centers and the full restarts.
reference_centers = 300;
component_cap = min(options.max_components, floor( ...
    reference_centers / options.minimum_observations_per_component));
if component_cap < 1
    work_units = 0;
    return
end
work_units = reference_centers * options.replicates * sum(1:component_cap);
end
