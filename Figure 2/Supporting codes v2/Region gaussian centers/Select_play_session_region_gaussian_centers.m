%% Select_play_session_region_gaussian_centers
% Compares mixture-selection rules for the fits saved by
% Estimate_play_session_region_gaussian_centers.m.
%
% AIC, BIC, the slope heuristic, and the elbow rules need the likelihood
% columns stored by every fit. ICL, NEC, CLC, AWE, overlap, and
% goodness-of-fit run when Estimate stored the newer measurement columns.
% A region with no fit is reported and skipped. The chosen curves in the
% estimate file stay as they were. This script writes a separate file.

%% Regions and rule settings
saving_folder = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\' ...
    'Analysis results\psth power by frequency and behavior'];

% Empty uses every saved region. SupCol and SC are the same region.
regions_to_select = {};

% Fourier is the spectrum used for the mid-LPAG mixture.
%   spectral_method = 'Welch';
%   spectral_method = {'Welch', 'Fourier'};
spectral_method = 'Fourier';

% Successive elbow: stop when the next gain is below this fraction of
% the first gain. Cumulative elbow: first count that has reached
% (1 - elbow_fraction) of the total drop.
elbow_fraction = 0.05;

% Largest component count that still keeps neighboring Gaussians apart.
ashman_d_minimum = 2;
bhattacharyya_maximum = 0.5;

% Smallest count whose fit is within this fraction of the best
% goodness-of-fit value.
gof_relative_tolerance = 0.1;

% NEC above this value keeps the one-component model.
nec_one_component_threshold = 1;

%% Load the saved mixtures
source_file = fullfile(saving_folder, ...
    'play_session_region_gaussian_center_mixtures.mat');
if exist(source_file, 'file') ~= 2
    error(['No saved region mixtures at %s. Run ' ...
        'Estimate_play_session_region_gaussian_centers.m first.'], ...
        source_file);
end
loaded = load(source_file, 'region_gaussian_centers');
if ~isfield(loaded, 'region_gaussian_centers')
    error('regionGaussianCenters:savedFile', ...
        '%s has no region_gaussian_centers.', source_file);
end
region_gaussian_centers = loaded.region_gaussian_centers;

if ischar(spectral_method) || isstring(spectral_method)
    method_names = cellstr(spectral_method);
else
    method_names = cellstr(string(spectral_method));
end
selection_options = struct( ...
    'elbow_fraction', elbow_fraction, ...
    'ashman_d_minimum', ashman_d_minimum, ...
    'bhattacharyya_maximum', bhattacharyya_maximum, ...
    'gof_relative_tolerance', gof_relative_tolerance, ...
    'nec_one_component_threshold', nec_one_component_threshold);

region_gaussian_selection = struct();
region_gaussian_selection.source_file = source_file;
region_gaussian_selection.options = selection_options;
region_gaussian_selection.methods = struct();

for method_n = 1:numel(method_names)
    method_name = canonical_selection_method(method_names{method_n});
    [selected, missing] = select_saved_region_gaussian_centers( ...
        region_gaussian_centers, regions_to_select, method_name);
    report_missing_regions(method_name, missing);
    region_rows = repmat(empty_selection_row(), numel(selected), 1);
    for region_n = 1:numel(selected)
        region_rows(region_n) = selection_for_region( ...
            selected(region_n), selection_options);
        fprintf('\n%s %s (%s measurements):\n', method_name, ...
            region_rows(region_n).name, region_rows(region_n).diagnostics);
        disp(region_rows(region_n).selection)
    end
    method_result = struct();
    method_result.by_region = region_rows;
    method_result.missing = missing;
    method_result.summary = selection_summary(region_rows);
    if ~isempty(region_rows)
        fprintf('\n%s selected component counts:\n', method_name);
        disp(method_result.summary)
    end
    region_gaussian_selection.methods.(method_name) = method_result;
end

selection_file = fullfile(saving_folder, ...
    'play_session_region_gaussian_center_selection.mat');
save(selection_file, 'region_gaussian_selection', '-v7.3');
fprintf('Saved %s\n', selection_file);

function row = selection_for_region(region_result, selection_options)
row = empty_selection_row();
row.name = region_result.name;
fit = region_result.fit;
if ~isstruct(fit) || ~isfield(fit, 'model_selection') || ...
        ~istable(fit.model_selection) || height(fit.model_selection) == 0
    row.diagnostics = 'none';
    row.selection = select_region_gaussian_criteria(table());
    return
end
options = selection_options;
if isfield(fit, 'center_frequencies')
    options.n_observations = numel(fit.center_frequencies);
end
row.selection = select_region_gaussian_criteria( ...
    fit.model_selection, options);
names = fit.model_selection.Properties.VariableNames;
if ismember('Entropy', names) && ...
        any(isfinite(fit.model_selection.Entropy))
    row.diagnostics = 'full';
else
    row.diagnostics = 'likelihood';
end
end

function row = empty_selection_row()
row = struct( ...
    'name', '', ...
    'diagnostics', '', ...
    'selection', table());
end

function summary = selection_summary(region_rows)
if isempty(region_rows)
    summary = table();
    return
end
criteria = region_rows(1).selection.Criterion;
if isempty(criteria)
    summary = table(cell(numel(region_rows), 1), ...
        'VariableNames', {'Region'});
    summary.Region = {region_rows.name}.';
    return
end
if isstring(criteria)
    criteria = cellstr(criteria);
end
values = nan(numel(region_rows), numel(criteria));
names = cell(numel(region_rows), 1);
for region_n = 1:numel(region_rows)
    names{region_n} = region_rows(region_n).name;
    this_selection = region_rows(region_n).selection;
    for criterion_n = 1:numel(criteria)
        if isempty(this_selection) || ~istable(this_selection)
            continue
        end
        match = strcmp(this_selection.Criterion, criteria{criterion_n});
        if any(match)
            values(region_n, criterion_n) = ...
                this_selection.SelectedN(find(match, 1, 'first'));
        end
    end
end
summary = table(names, 'VariableNames', {'Region'});
for criterion_n = 1:numel(criteria)
    summary.(criteria{criterion_n}) = values(:, criterion_n);
end
end

function name = canonical_selection_method(token)
token = strtrim(char(string(token)));
if strcmpi(token, 'Fourier')
    name = 'Fourier';
elseif strcmpi(token, 'Welch')
    name = 'Welch';
else
    error('spectral_method must be Fourier, Welch, or both.');
end
end

function report_missing_regions(method_name, missing)
if isempty(missing)
    fprintf('%s: every requested region has a saved fit.\n', method_name);
    return
end
fprintf('%s regions skipped:\n', method_name);
for region_n = 1:numel(missing)
    fprintf('  %s: %s.\n', missing(region_n).name, missing(region_n).reason);
end
end
