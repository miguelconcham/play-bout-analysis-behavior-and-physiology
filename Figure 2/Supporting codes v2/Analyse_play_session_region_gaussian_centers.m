%% Analyse_play_session_region_gaussian_centers
% Extends the mid-LPAG mixture of FOOOF Gaussian centers to every
% recorded channel. The centers are the periodic power bumps fit above
% the aperiodic spectrum. Each listed region is pooled across channels,
% shanks, recordings, and partner sessions, then fit on its own.
%
% The figure has one row per region. BIC and AIC are both drawn. The
% mixture on the distribution is the one selected by selection_criterion,
% which is BIC unless that switch is set to AIC.
%
% A channel labeled with two listed regions contributes its bumps to both.
% Labels outside the list, including InfCol, are left out. SupCol and SC
% are the same region. VLPAG is not treated as LPAG.
%
% Requires the channel-wise spectra written by
% Estimate_play_session_probe_spectra.m.

%% Regions, spectrum, and model selection
saving_folder = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\' ...
    'Analysis results\psth power by frequency and behavior'];

% Order used for the figure and the saved summary.
regions_to_analyze = {'DR', 'VLPAG', 'LPAG', 'DLPAG', 'SupCol'};

% Fourier is the spectrum used for the mid-LPAG mixture. Welch is the
% other saved parameterization:
%   spectral_method = 'Welch';
%   spectral_method = {'Welch', 'Fourier'};
spectral_method = 'Fourier';

% Lower is better. AIC is plotted either way.
selection_criterion = 'BIC';

% Same component search as the mid-LPAG mixture.
mixture_options = struct( ...
    'max_components', 10, ...
    'replicates', 50, ...
    'max_iter', 2000, ...
    'regularization', 1e-4, ...
    'minimum_observations_per_component', 5, ...
    'minimum_effective_membership', 3, ...
    'random_seed', 1);

% Drawn on the distribution only. It does not choose the mixture.
% Set this empty to omit the line.
low_frequency_threshold_hz = 8;

histogram_bin_hz = 1;

% Empty keeps every saved recording and every partner session.
recording_index = [];
session_index = [];

mixture_options.selection_criterion = selection_criterion;

%% Load every selected channel-wise spectrum
file_patterns = probe_periodic_spectra_filename();
file_list = dir(fullfile(saving_folder, file_patterns{1}));
if ~isempty(file_list)
    keep = false(numel(file_list), 1);
    for file_n = 1:numel(file_list)
        [~, ~, extension] = fileparts(file_list(file_n).name);
        keep(file_n) = ~file_list(file_n).isdir && ...
            (isempty(extension) || strcmpi(extension, '.mat'));
    end
    file_list = file_list(keep);
end
if isempty(file_list)
    error(['No channel-wise spectra found in %s. Expected names such as ' ...
        'play_session_probe_periodic_spectra_B1D1_1007_Dual.mat.'], ...
        saving_folder);
end
[~, order] = sort({file_list.name});
file_list = file_list(order);

if isempty(recording_index)
    recording_index = 1:numel(file_list);
end
if any(recording_index < 1) || any(recording_index > numel(file_list))
    error('recording_index must be between 1 and %d.', numel(file_list));
end

fprintf('Saved recordings:\n');
for file_n = 1:numel(file_list)
    fprintf('  %d  %s\n', file_n, file_list(file_n).name);
end

frequency = [];
regions = {};
pieces = cell(numel(recording_index), 1);
piece_n = 0;
for file_n = recording_index(:).'
    input_file = fullfile(saving_folder, file_list(file_n).name);
    fprintf('Loading %s\n', input_file);
    loaded = load(input_file, 'probe_spectra');
    probe_spectra = loaded.probe_spectra;
    this_frequency = probe_spectra.frequency(:).';
    if isempty(frequency)
        frequency = this_frequency;
    elseif numel(frequency) ~= numel(this_frequency) || ...
            max(abs(frequency - this_frequency)) > 1e-9
        error('Frequency grids differ between saved recordings.');
    end
    piece_n = piece_n + 1;
    [pieces{piece_n}, regions] = region_gaussian_center_table( ...
        probe_spectra, regions_to_analyze, spectral_method, ...
        session_index, file_n);
end
pieces = pieces(1:piece_n);
pieces = pieces(~cellfun(@isempty, pieces));
if isempty(pieces)
    peak_table = region_gaussian_center_table(struct( ...
        'animal_name', '', 'n_sessions', 0), regions, spectral_method);
else
    peak_table = vertcat(pieces{:});
end

if ischar(spectral_method) || isstring(spectral_method)
    method_names = cellstr(spectral_method);
else
    method_names = cellstr(string(spectral_method));
end
if isempty(frequency)
    frequency_limits = [0.1 45];
else
    frequency_limits = frequency([1 end]);
end

%% Fit and draw one mixture per region
region_gaussian_centers = struct();
region_gaussian_centers.regions = regions;
region_gaussian_centers.frequency = frequency;
region_gaussian_centers.selection_criterion = upper(strtrim(char( ...
    string(selection_criterion))));
region_gaussian_centers.low_frequency_threshold_hz = ...
    low_frequency_threshold_hz;
region_gaussian_centers.histogram_bin_hz = histogram_bin_hz;
region_gaussian_centers.mixture_options = mixture_options;
region_gaussian_centers.peak_table = peak_table;
region_gaussian_centers.methods = struct();

for method_n = 1:numel(method_names)
    method_name = canonical_method_name(method_names{method_n});
    in_method = strcmp(peak_table.Method, method_name);
    region_results = repmat(empty_region_result(), numel(regions), 1);
    for region_n = 1:numel(regions)
        in_region = in_method & strcmp(peak_table.Region, regions{region_n});
        centers = peak_table.CenterFrequency(in_region);
        region_results(region_n).name = regions{region_n};
        region_results(region_n).n_peaks = numel(centers);
        region_results(region_n).n_channels = numel(unique([ ...
            peak_table.RecordingIndex(in_region), ...
            peak_table.ProbeNumber(in_region), ...
            peak_table.Channel(in_region)], 'rows'));
        region_results(region_n).n_sessions = numel(unique([ ...
            peak_table.RecordingIndex(in_region), ...
            peak_table.SessionIndex(in_region)], 'rows'));
        region_results(region_n).fit = fit_gaussian_center_mixture( ...
            centers, frequency_limits, mixture_options);
    end

    plot_region_gaussian_centers(region_results, frequency_limits, ...
        region_gaussian_centers.selection_criterion, ...
        low_frequency_threshold_hz, histogram_bin_hz, method_name);

    summary = region_summary_table(region_results);
    fprintf('\n%s Gaussian-center mixtures (%s):\n', method_name, ...
        region_gaussian_centers.selection_criterion);
    disp(summary)
    for region_n = 1:numel(region_results)
        if strcmp(region_results(region_n).fit.status, 'fitted')
            fprintf('\n%s %s components:\n', method_name, ...
                region_results(region_n).name);
            disp(region_results(region_n).fit.best_component_table)
        end
    end

    method_result = struct();
    method_result.summary = summary;
    method_result.by_region = region_results;
    region_gaussian_centers.methods.(method_name) = method_result;
end

output_file = fullfile(saving_folder, ...
    'play_session_region_gaussian_center_mixtures.mat');
save(output_file, 'region_gaussian_centers', '-v7.3');
disp(['Saved ', output_file])

function name = canonical_method_name(token)
token = strtrim(char(string(token)));
if strcmpi(token, 'Fourier')
    name = 'Fourier';
elseif strcmpi(token, 'Welch')
    name = 'Welch';
else
    error('spectral_method must be Fourier, Welch, or both.');
end
end

function region_result = empty_region_result()
region_result = struct( ...
    'name', '', ...
    'n_peaks', 0, ...
    'n_channels', 0, ...
    'n_sessions', 0, ...
    'fit', []);
end

function summary = region_summary_table(region_results)
n_regions = numel(region_results);
names = cell(n_regions, 1);
n_peaks = zeros(n_regions, 1);
n_channels = zeros(n_regions, 1);
n_sessions = zeros(n_regions, 1);
best_n = nan(n_regions, 1);
selected_bic = nan(n_regions, 1);
selected_aic = nan(n_regions, 1);
status = cell(n_regions, 1);
for region_n = 1:n_regions
    this_region = region_results(region_n);
    names{region_n} = this_region.name;
    n_peaks(region_n) = this_region.n_peaks;
    n_channels(region_n) = this_region.n_channels;
    n_sessions(region_n) = this_region.n_sessions;
    best_n(region_n) = this_region.fit.best_n_components;
    status{region_n} = this_region.fit.status;
    model_table = this_region.fit.model_selection;
    if isfinite(best_n(region_n)) && ~isempty(model_table)
        best_row = model_table.NComponents == best_n(region_n);
        selected_bic(region_n) = model_table.BIC(best_row);
        selected_aic(region_n) = model_table.AIC(best_row);
    end
end
summary = table(names, n_peaks, n_channels, n_sessions, best_n, ...
    selected_bic, selected_aic, status, ...
    'VariableNames', {'Region', 'NCenters', 'NChannels', 'NSessions', ...
    'BestNComponents', 'SelectedBIC', 'SelectedAIC', 'Status'});
end
