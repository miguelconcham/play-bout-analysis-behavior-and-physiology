function options = gaussian_mixture_options(options)
%GAUSSIAN_MIXTURE_OPTIONS  Defaults for the regional center mixtures.
%
% max_components is a hard cap. The regional search tests every count
% from 1 through 20. replicates is the restart count for every region.
% large_sample_replicates is used only when it is smaller than replicates
% and a region has more than large_sample_centers. The default keeps all
% 50 restarts so a large region gets the same search as a small one.
% slow_run_seconds is the clock budget. A projection past it is reported
% while the fit is running.

if nargin < 1 || isempty(options)
    options = struct();
end
options = fill_mixture_option(options, 'selection_criterion', 'BIC');
options = fill_mixture_option(options, 'max_components', 20);
options = fill_mixture_option(options, 'replicates', 50);
options = fill_mixture_option(options, 'large_sample_centers', 500);
options = fill_mixture_option(options, 'large_sample_replicates', 50);
options = fill_mixture_option(options, 'max_iter', 2000);
options = fill_mixture_option(options, 'regularization', 1e-4);
options = fill_mixture_option(options, ...
    'minimum_observations_per_component', 5);
options = fill_mixture_option(options, 'minimum_effective_membership', 3);
options = fill_mixture_option(options, 'random_seed', 1);
options = fill_mixture_option(options, 'slow_run_seconds', 10 * 3600);
options = fill_mixture_option(options, 'progress_label', '');
options = fill_mixture_option(options, 'other_regions_work_units', 0);

criterion = upper(strtrim(char(string(options.selection_criterion))));
if ~ismember(criterion, {'BIC', 'AIC'})
    error('regionGaussianCenters:criterion', ...
        'selection_criterion must be BIC or AIC.');
end
options.selection_criterion = criterion;
if options.max_components < 1 || options.replicates < 1 || ...
        options.large_sample_replicates < 1 || options.max_iter < 1 || ...
        options.minimum_observations_per_component < 1 || ...
        options.minimum_effective_membership < 1
    error('regionGaussianCenters:options', ...
        ['max_components, replicates, large_sample_replicates, ' ...
        'max_iter, and the membership minima must be at least 1.']);
end
if options.large_sample_centers < 0 || options.slow_run_seconds < 0 || ...
        options.other_regions_work_units < 0
    error('regionGaussianCenters:options', ...
        ['large_sample_centers, slow_run_seconds, and ' ...
        'other_regions_work_units must be non-negative.']);
end
end

function options = fill_mixture_option(options, name, default_value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = default_value;
end
end
