%% Analyse_play_session_region_gaussian_centers
% Draws the Gaussian-center mixtures saved by
% Estimate_play_session_region_gaussian_centers.m.
%
% regions_to_plot chooses the rows. Leave it empty to draw every saved
% region. A region that has no fit yet is named in the command window and
% left out of the figure. SupCol and SC are the same region.
%
% selection_criterion chooses which saved rule supplies the mixture.
% Leave it empty to draw the criterion stored by Estimate, usually BIC.
% After Select_play_session_region_gaussian_centers.m, set it to one of
% the printed names, for example 'ICL' or 'SlopeHeuristicRegression'.
% The elbow and overlap settings below must match that script.

%% Regions to draw
saving_folder = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\' ...
    'Analysis results\psth power by frequency and behavior'];

% Example that draws the finished regions while SupCol is still missing:
%   regions_to_plot = {'DR', 'VLPAG', 'LPAG', 'DLPAG'};
regions_to_plot = {'DR', 'VLPAG', 'LPAG', 'DLPAG', 'SupCol'};

% Fourier is the spectrum used for the mid-LPAG mixture. Welch is the
% other saved parameterization:
%   spectral_method = 'Welch';
%   spectral_method = {'Welch', 'Fourier'};
spectral_method = 'Fourier';

% '' draws the criterion saved by Estimate. Otherwise use a name from
% Select_play_session_region_gaussian_centers.m:
%   selection_criterion = 'ICL';
%   selection_criterion = 'SlopeHeuristicRegression';
%   selection_criterion = 'ElbowBICSuccessive';
selection_criterion = '';

% Keep these equal to the selection script when the chosen rule uses them.
elbow_fraction = 0.05;
ashman_d_minimum = 2;
bhattacharyya_maximum = 0.5;
gof_relative_tolerance = 0.1;
nec_one_component_threshold = 1;

% NaN uses the value stored by Estimate. [] omits the reference line.
% A number such as 8 draws that line. It does not choose the mixture.
low_frequency_threshold_hz = NaN;

% NaN uses the bin width stored by Estimate.
histogram_bin_hz = NaN;

%% Load the saved mixtures
output_file = fullfile(saving_folder, ...
    'play_session_region_gaussian_center_mixtures.mat');
if exist(output_file, 'file') ~= 2
    error(['No saved region mixtures at %s. Run ' ...
        'Estimate_play_session_region_gaussian_centers.m first.'], ...
        output_file);
end
loaded = load(output_file, 'region_gaussian_centers');
if ~isfield(loaded, 'region_gaussian_centers')
    error('regionGaussianCenters:savedFile', ...
        '%s has no region_gaussian_centers.', output_file);
end
region_gaussian_centers = loaded.region_gaussian_centers;

if ischar(spectral_method) || isstring(spectral_method)
    method_names = cellstr(spectral_method);
else
    method_names = cellstr(string(spectral_method));
end
frequency = [];
if isfield(region_gaussian_centers, 'frequency')
    frequency = region_gaussian_centers.frequency;
end
if isempty(frequency)
    frequency_limits = [0.1 45];
else
    frequency_limits = frequency([1 end]);
end
if isempty(selection_criterion)
    if isfield(region_gaussian_centers, 'selection_criterion') && ...
            ~isempty(region_gaussian_centers.selection_criterion)
        selection_criterion = region_gaussian_centers.selection_criterion;
    else
        selection_criterion = 'BIC';
    end
end
selection_options = struct( ...
    'elbow_fraction', elbow_fraction, ...
    'ashman_d_minimum', ashman_d_minimum, ...
    'bhattacharyya_maximum', bhattacharyya_maximum, ...
    'gof_relative_tolerance', gof_relative_tolerance, ...
    'nec_one_component_threshold', nec_one_component_threshold);
threshold_hz = resolve_plot_threshold( ...
    low_frequency_threshold_hz, region_gaussian_centers);
bin_hz = resolve_histogram_bin(histogram_bin_hz, region_gaussian_centers);

%% Draw the requested regions that have a fit
for method_n = 1:numel(method_names)
    method_name = method_names{method_n};
    [selected, missing] = select_saved_region_gaussian_centers( ...
        region_gaussian_centers, regions_to_plot, method_name);
    report_missing_regions(method_name, missing);
    if isempty(selected)
        fprintf('%s: no saved region to plot.\n', method_name);
        continue
    end
    plotted = selected([]);
    for region_n = 1:numel(selected)
        [updated, outcome] = region_mixture_for_criterion( ...
            selected(region_n), selection_criterion, frequency_limits, ...
            selection_options);
        if ~outcome.applied
            fprintf('  %s left out: %s.\n', selected(region_n).name, ...
                outcome.note);
            continue
        end
        fprintf('  %s: %s selects %g.\n', updated.name, ...
            outcome.criterion, outcome.selected_n);
        if isempty(plotted)
            plotted = updated;
        else
            plotted(end + 1, 1) = updated; %#ok<AGROW>
        end
    end
    if isempty(plotted)
        fprintf('%s: no region can be drawn with %s.\n', ...
            method_name, char(string(selection_criterion)));
        continue
    end
    fprintf('%s figure: %s.\n', method_name, strjoin({plotted.name}, ', '));
    plot_region_gaussian_centers(plotted, frequency_limits, ...
        selection_criterion, threshold_hz, bin_hz, method_name);
    for region_n = 1:numel(plotted)
        this_region = plotted(region_n);
        best_n = NaN;
        if isfield(this_region.fit, 'best_n_components')
            best_n = this_region.fit.best_n_components;
        end
        fprintf('  %s: %s, %d centers', this_region.name, ...
            this_region.fit.status, this_region.n_peaks);
        if isfinite(best_n)
            fprintf(', %g Gaussians', best_n);
        end
        fprintf('.\n');
    end
    drawnow;
end

function threshold_hz = resolve_plot_threshold(setting, result)
if isnumeric(setting) && isempty(setting)
    threshold_hz = [];
    return
end
if isnumeric(setting) && all(isnan(setting(:)))
    if isfield(result, 'low_frequency_threshold_hz')
        threshold_hz = result.low_frequency_threshold_hz;
    else
        threshold_hz = [];
    end
    return
end
threshold_hz = setting;
end

function bin_hz = resolve_histogram_bin(setting, result)
if isempty(setting) || (isnumeric(setting) && all(isnan(setting(:))))
    if isfield(result, 'histogram_bin_hz') && ~isempty(result.histogram_bin_hz)
        bin_hz = result.histogram_bin_hz;
    else
        bin_hz = 1;
    end
    return
end
bin_hz = setting;
end

function report_missing_regions(method_name, missing)
if isempty(missing)
    fprintf('%s: every requested region has a saved fit.\n', method_name);
    return
end
fprintf('%s regions left out of the figure:\n', method_name);
for region_n = 1:numel(missing)
    fprintf('  %s: %s.\n', missing(region_n).name, missing(region_n).reason);
end
end
