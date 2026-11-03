function options = actor_behavior_probability_options(overrides)
%ACTOR_BEHAVIOR_PROBABILITY_OPTIONS Settings for the actor-relative pipeline.
%
% This pipeline is deliberately independent of Figure 1/Behavior probability.
% target_behaviors define the events around which behavior is sampled.

if nargin < 1
    overrides = struct();
end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('actor_behavior_probability_options:Overrides', ...
        'Overrides must be a scalar structure.');
end

data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data'];

options = struct();
options.target_behaviors = {'Pounce', 'PounceI'};
options.all_behaviors = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', ...
    'Escape', 'CB', 'CD', 'Grooming', 'Rearing', 'Sniffing', ...
    'Scratching', 'PounceI', 'Bite'};
options.time_window = [-5 5];
options.bin_size = 0.01;
options.time = [];
options.comparison_window = [0 2];
options.flank_seconds = 2;
options.wrapped_bins = 20;
options.pre_lag = [];
options.post_lag = [];
options.wrap_fraction = [];
options.wrapped_plot_time = [];
options.behavior_folder = fullfile(data_root, 'Behavior backups');
options.output_folder = fullfile(data_root, 'Analysis results', ...
    'behavior probability by target actor');
options.output_filename = 'actor_behavior_probability_around_offset.mat';
options.alignment = 'target behavior offset';
options.self_definition = 'the animal that performed the target event';
options.other_definition = ...
    'logical OR across all other animals in the same Partners session';

names = fieldnames(overrides);
for field_n = 1:numel(names)
    if ~isfield(options, names{field_n})
        error('actor_behavior_probability_options:UnknownOption', ...
            'Unknown option "%s".', names{field_n});
    end
    options.(names{field_n}) = overrides.(names{field_n});
end

options.target_behaviors = normalize_names(options.target_behaviors);
options.all_behaviors = normalize_names(options.all_behaviors);
if isempty(options.target_behaviors)
    error('actor_behavior_probability_options:Targets', ...
        'target_behaviors cannot be empty.');
end
if numel(unique(options.all_behaviors)) ~= numel(options.all_behaviors)
    error('actor_behavior_probability_options:Behaviors', ...
        'all_behaviors contains duplicate names.');
end
unknown = setdiff(options.target_behaviors, options.all_behaviors);
if ~isempty(unknown)
    error('actor_behavior_probability_options:Targets', ...
        'Unknown target behavior(s): %s.', strjoin(unknown, ', '));
end

if ~isnumeric(options.time_window) || numel(options.time_window) ~= 2 || ...
        any(~isfinite(options.time_window)) || ...
        options.time_window(1) > 0 || options.time_window(2) < 0 || ...
        options.time_window(2) <= options.time_window(1)
    error('actor_behavior_probability_options:Window', ...
        'time_window must be a finite [negative_pre positive_post] pair.');
end
options.time_window = double(options.time_window(:).');
options.bin_size = positive_scalar(options.bin_size, 'bin_size');
n_steps = diff(options.time_window) / options.bin_size;
if abs(n_steps - round(n_steps)) > 1e-9 * max(1, abs(n_steps))
    error('actor_behavior_probability_options:WindowBins', ...
        'The time-window length must be an integer multiple of bin_size.');
end
options.time = options.time_window(1) + ...
    (0:round(n_steps)) * options.bin_size;
options.time(end) = options.time_window(2);
time_tolerance = max(eps(max(1, max(abs(options.time)))) * 16, ...
    options.bin_size * 1e-9);
options.time(abs(options.time) <= time_tolerance) = 0;

if ~isnumeric(options.comparison_window) || ...
        numel(options.comparison_window) ~= 2 || ...
        any(~isfinite(options.comparison_window)) || ...
        options.comparison_window(1) < 0 || ...
        options.comparison_window(2) <= options.comparison_window(1) || ...
        options.comparison_window(2) > options.time_window(2)
    error('actor_behavior_probability_options:ComparisonWindow', ...
        ['comparison_window must be an ascending post-offset range ' ...
        'inside time_window.']);
end
options.comparison_window = double(options.comparison_window(:).');
if ~any(options.time > time_tolerance & ...
        options.time >= options.comparison_window(1) - time_tolerance & ...
        options.time <= options.comparison_window(2) + time_tolerance)
    error('actor_behavior_probability_options:ComparisonWindowBins', ...
        'comparison_window contains no strictly post-offset sample.');
end

options.flank_seconds = positive_scalar( ...
    options.flank_seconds, 'flank_seconds');
n_flank = options.flank_seconds / options.bin_size;
if abs(n_flank - round(n_flank)) > 1e-9 * max(1, abs(n_flank))
    error('actor_behavior_probability_options:FlankBins', ...
        'flank_seconds must be an integer multiple of bin_size.');
end
n_flank = round(n_flank);
if ~isnumeric(options.wrapped_bins) || ...
        ~isscalar(options.wrapped_bins) || ...
        ~isfinite(options.wrapped_bins) || options.wrapped_bins < 2 || ...
        options.wrapped_bins ~= round(options.wrapped_bins)
    error('actor_behavior_probability_options:WrappedBins', ...
        'wrapped_bins must be an integer of at least 2.');
end
options.wrapped_bins = double(options.wrapped_bins);
options.pre_lag = -options.flank_seconds + ...
    (0:n_flank - 1) * options.bin_size;
options.post_lag = (1:n_flank) * options.bin_size;
options.wrap_fraction = linspace(0, 1, options.wrapped_bins);
options.wrapped_plot_time = [options.pre_lag, ...
    linspace(0, options.flank_seconds, options.wrapped_bins), ...
    options.flank_seconds + options.post_lag];
end

function value = positive_scalar(value, name)
if ~isnumeric(value) || ~isscalar(value) || ...
        ~isfinite(value) || value <= 0
    error('actor_behavior_probability_options:PositiveScalar', ...
        '%s must be a positive finite scalar.', name);
end
value = double(value);
end

function names = normalize_names(values)
if ischar(values)
    names = {values};
elseif isstring(values) || iscategorical(values)
    names = cellstr(values(:).');
elseif iscell(values)
    names = cell(size(values));
    for value_n = 1:numel(values)
        if ~(ischar(values{value_n}) || ...
                (isstring(values{value_n}) && isscalar(values{value_n})))
            error('actor_behavior_probability_options:BehaviorNames', ...
                'Behavior names must be text scalars.');
        end
        names{value_n} = char(values{value_n});
    end
    names = names(:).';
else
    error('actor_behavior_probability_options:BehaviorNames', ...
        'Behavior names must be text.');
end
names = strtrim(names);
if any(cellfun(@isempty, names))
    error('actor_behavior_probability_options:BehaviorNames', ...
        'Behavior names cannot be empty.');
end
end
