function options = behavior_probability_options(overrides)
%BEHAVIOR_PROBABILITY_OPTIONS Settings for behavior probability matrices.
%
% target_behaviors contains the Type2 behaviors whose offsets define time
% zero. time_window is [pre post] in seconds relative to that offset.

if nargin < 1
    overrides = struct();
end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('behavior_probability_options:Overrides', ...
        'Overrides must be a scalar structure.');
end

data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data'];

options = struct();
options.target_behaviors = {'Pounce'};
options.all_behaviors = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', ...
    'Escape', 'CB', 'CD', 'Grooming', 'Rearing', 'Sniffing', ...
    'Scratching', 'PounceI', 'Bite'};
options.time_window = [-5 5];
options.bin_size = 0.01;
options.time = [];
options.comparison_window = [];
options.data_root = data_root;
options.behavior_folder = fullfile(data_root, 'Behavior backups');
options.output_folder = fullfile(data_root, 'Analysis results', ...
    'behavior probability');
options.output_filename = 'behavior_probability_around_offset.mat';
options.alignment = 'offset';
options.self_definition = ...
    'animal that performed the target event on this matrix row';
options.partner_definition = ...
    'other animal(s) in the same partner-play session';

names = fieldnames(overrides);
for field_n = 1:numel(names)
    if ~isfield(options, names{field_n})
        error('behavior_probability_options:UnknownOption', ...
            'Unknown option "%s".', names{field_n});
    end
    options.(names{field_n}) = overrides.(names{field_n});
end

options.target_behaviors = behavior_names(options.target_behaviors);
options.all_behaviors = behavior_names(options.all_behaviors);

if isempty(options.target_behaviors)
    error('behavior_probability_options:Targets', ...
        'target_behaviors cannot be empty.');
end
if numel(unique(options.all_behaviors)) ~= numel(options.all_behaviors)
    error('behavior_probability_options:Behaviors', ...
        'all_behaviors contains duplicate names.');
end
unknown = setdiff(options.target_behaviors, options.all_behaviors);
if ~isempty(unknown)
    error('behavior_probability_options:Targets', ...
        'Unknown target behavior(s): %s. Valid Type2 names are: %s.', ...
        strjoin(unknown, ', '), strjoin(options.all_behaviors, ', '));
end
if ~isnumeric(options.time_window) || numel(options.time_window) ~= 2 || ...
        any(~isfinite(options.time_window)) || ...
        options.time_window(1) > 0 || options.time_window(2) < 0 || ...
        options.time_window(2) <= options.time_window(1)
    error('behavior_probability_options:Window', ...
        'time_window must be a finite [negative_pre positive_post] pair.');
end
options.time_window = double(options.time_window(:).');
if ~isscalar(options.bin_size) || ~isnumeric(options.bin_size) || ...
        ~isfinite(options.bin_size) || options.bin_size <= 0
    error('behavior_probability_options:BinSize', ...
        'bin_size must be a positive finite scalar.');
end
options.bin_size = double(options.bin_size);
if ~strcmp(options.alignment, 'offset')
    error('behavior_probability_options:Alignment', ...
        'This pipeline currently supports offset alignment only.');
end
if ~isempty(options.comparison_window)
    if ~isnumeric(options.comparison_window) || ...
            numel(options.comparison_window) ~= 2 || ...
            any(~isfinite(options.comparison_window)) || ...
            options.comparison_window(1) < 0 || ...
            options.comparison_window(2) <= options.comparison_window(1) || ...
            options.comparison_window(2) > options.time_window(2)
        error('behavior_probability_options:ComparisonWindow', ...
            ['comparison_window must be an ascending post-offset range ' ...
            'inside time_window.']);
    end
    options.comparison_window = double(options.comparison_window(:).');
end

n_steps = diff(options.time_window) / options.bin_size;
if abs(n_steps - round(n_steps)) > 1e-9 * max(1, abs(n_steps))
    error('behavior_probability_options:WindowBins', ...
        'The time-window length must be an integer multiple of bin_size.');
end
options.time = options.time_window(1) + ...
    (0:round(n_steps)) * options.bin_size;
options.time(end) = options.time_window(2);
time_tolerance = max(eps(max(1, max(abs(options.time)))) * 16, ...
    options.bin_size * 1e-9);
options.time(abs(options.time) <= time_tolerance) = 0;
if ~isempty(options.comparison_window) && ...
        ~any(options.time > time_tolerance & ...
        options.time >= options.comparison_window(1) - time_tolerance & ...
        options.time <= options.comparison_window(2) + time_tolerance)
    error('behavior_probability_options:ComparisonWindowBins', ...
        ['comparison_window does not contain a strictly post-offset ' ...
        'sample at the selected bin_size.']);
end
end

function names = behavior_names(values)
if ischar(values)
    names = {values};
elseif isstring(values) || iscategorical(values)
    names = cellstr(values(:).');
elseif iscell(values)
    names = cell(size(values));
    for value_n = 1:numel(values)
        if ~(ischar(values{value_n}) || ...
                (isstring(values{value_n}) && isscalar(values{value_n})))
            error('behavior_probability_options:BehaviorNames', ...
                'Behavior names must be text scalars.');
        end
        names{value_n} = char(values{value_n});
    end
    names = names(:).';
else
    error('behavior_probability_options:BehaviorNames', ...
        'Behavior names must be a cell array, string array, or character vector.');
end
names = strtrim(names);
if any(cellfun(@isempty, names))
    error('behavior_probability_options:BehaviorNames', ...
        'Behavior names cannot be empty.');
end
end
