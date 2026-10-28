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

% Every count from 1 through max_components is fit, with 50 restarts.
% The previous cap of 10 was reached, so the cap is 20. Large regions
% keep all 50 restarts: large_sample_replicates matches replicates.
% Set large_sample_replicates lower only to shorten a region that has
% more than large_sample_centers. A projection past slow_run_seconds
% (10 h) is reported while the fit is running. Each finished region is
% written to the output file before the next region starts.
%
% resume_partial_fit keeps regions already stored in the output file and
% fits only regions whose fit is empty. The saved peak table is reused.
% Set it false to load the spectra and fit every region again.
mixture_options = struct( ...
    'max_components', 20, ...
    'replicates', 50, ...
    'large_sample_centers', 500, ...
    'large_sample_replicates', 50, ...
    'max_iter', 2000, ...
    'regularization', 1e-4, ...
    'minimum_observations_per_component', 5, ...
    'minimum_effective_membership', 3, ...
    'random_seed', 1, ...
    'slow_run_seconds', 10 * 3600);

% Drawn on the distribution only. It does not choose the mixture.
% Set this empty to omit the line.
low_frequency_threshold_hz = 8;

histogram_bin_hz = 1;

% Empty keeps every saved recording and every partner session.
recording_index = [];
session_index = [];

% true continues an interrupted run. Finished regions stay as they are.
resume_partial_fit = true;

mixture_options.selection_criterion = selection_criterion;
output_file = fullfile(saving_folder, ...
    'play_session_region_gaussian_center_mixtures.mat');
saved_for_resume = [];
resumed = false;
if resume_partial_fit && exist(output_file, 'file') == 2
    fprintf('Checking %s for finished regions.\n', output_file);
    drawnow;
    loaded_result = load(output_file, 'region_gaussian_centers');
    if isfield(loaded_result, 'region_gaussian_centers')
        saved_for_resume = loaded_result.region_gaussian_centers;
        if saved_mixture_has_peaks(saved_for_resume)
            region_gaussian_centers = saved_for_resume;
            resumed = true;
            fprintf(['Resuming. The saved peak table is reused, so the ' ...
                'spectra are not loaded again.\n']);
        else
            fprintf(['The saved file has no peak table. The spectra will ' ...
                'be loaded, and any finished region fits will still be kept.\n']);
        end
        drawnow;
    end
end

if resumed
    regions = region_gaussian_centers.regions;
    if isstring(regions)
        regions = cellstr(regions);
    end
    regions = regions(:).';
    peak_table = region_gaussian_centers.peak_table;
    frequency = region_gaussian_centers.frequency(:).';
    mixture_options = region_gaussian_centers.mixture_options;
    if isfield(region_gaussian_centers, 'selection_criterion') && ...
            ~isempty(region_gaussian_centers.selection_criterion)
        selection_criterion = region_gaussian_centers.selection_criterion;
    end
    if isfield(region_gaussian_centers, 'low_frequency_threshold_hz')
        low_frequency_threshold_hz = ...
            region_gaussian_centers.low_frequency_threshold_hz;
    end
    if isfield(region_gaussian_centers, 'histogram_bin_hz') && ...
            ~isempty(region_gaussian_centers.histogram_bin_hz)
        histogram_bin_hz = region_gaussian_centers.histogram_bin_hz;
    end
else
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
selected_files = recording_index(:).';
total_bytes = sum([file_list(selected_files).bytes]);
bytes_done = 0;
load_clock = tic;
fprintf(['Loading %d channel-wise spectra (%.0f MB). Each file holds ' ...
    'every channel, so this step can take a few minutes.\n'], ...
    numel(selected_files), total_bytes / 1e6);
for file_n = selected_files
    input_file = fullfile(saving_folder, file_list(file_n).name);
    this_bytes = file_list(file_n).bytes;
    fprintf('Loading spectra %d/%d: %s (%.0f MB).\n', ...
        piece_n + 1, numel(selected_files), file_list(file_n).name, ...
        this_bytes / 1e6);
    if bytes_done > 0
        load_left = mixture_remaining_seconds( ...
            toc(load_clock), bytes_done, total_bytes - bytes_done);
        fprintf('  About %s left for loading.\n', ...
            format_search_duration(load_left));
    end
    drawnow;
    file_clock = tic;
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
    clear loaded probe_spectra
    bytes_done = bytes_done + this_bytes;
    fprintf('  Peaks collected in %s. Loading elapsed %s.\n', ...
        format_search_duration(toc(file_clock)), ...
        format_search_duration(toc(load_clock)));
    drawnow;
end
pieces = pieces(1:piece_n);
pieces = pieces(~cellfun(@isempty, pieces));
if isempty(pieces)
    peak_table = region_gaussian_center_table(struct( ...
        'animal_name', '', 'n_sessions', 0), regions, spectral_method);
else
    peak_table = vertcat(pieces{:});
end
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
if ~resumed
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
    if ~isempty(saved_for_resume) && isfield(saved_for_resume, 'methods')
        region_gaussian_centers.methods = saved_for_resume.methods;
    end
end

for method_n = 1:numel(method_names)
    method_name = canonical_method_name(method_names{method_n});
    prior_results = saved_method_regions(region_gaussian_centers, method_name);
    if method_is_complete(region_gaussian_centers, method_name)
        fprintf(['%s mixtures are already complete. ' ...
            'Drawing the saved figure.\n'], method_name);
        plot_region_gaussian_centers(prior_results, frequency_limits, ...
            region_gaussian_centers.selection_criterion, ...
            low_frequency_threshold_hz, histogram_bin_hz, method_name);
        drawnow;
        continue
    end
    in_method = strcmp(peak_table.Method, method_name);
    region_results = repmat(empty_region_result(), numel(regions), 1);
    region_centers = cell(numel(regions), 1);
    search_plans = cell(numel(regions), 1);
    region_kept = false(numel(regions), 1);
    for region_n = 1:numel(regions)
        in_region = in_method & strcmp(peak_table.Region, regions{region_n});
        centers = peak_table.CenterFrequency(in_region);
        centers = centers(isfinite(centers));
        region_centers{region_n} = centers;
        region_results(region_n).name = regions{region_n};
        region_results(region_n).n_peaks = numel(centers);
        region_results(region_n).n_channels = numel(unique([ ...
            peak_table.RecordingIndex(in_region), ...
            peak_table.ProbeNumber(in_region), ...
            peak_table.Channel(in_region)], 'rows'));
        region_results(region_n).n_sessions = numel(unique([ ...
            peak_table.RecordingIndex(in_region), ...
            peak_table.SessionIndex(in_region)], 'rows'));
        search_plans{region_n} = gaussian_mixture_search_plan( ...
            numel(centers), numel(unique(centers)), mixture_options);
        kept_fit = saved_fit_for_region(prior_results, regions{region_n});
        region_kept(region_n) = ~isempty(kept_fit);
        if region_kept(region_n)
            region_results(region_n).fit = kept_fit;
            search_plans{region_n}.already_saved = true;
            search_plans{region_n}.work_units = 0;
            search_plans{region_n}.runtime_note = ...
                'already saved; this fit is kept';
        else
            search_plans{region_n}.already_saved = false;
        end
    end
    print_mixture_plan(method_name, regions, search_plans, mixture_options);

    for region_n = 1:numel(regions)
        if region_kept(region_n)
            fprintf('%s %s (%d/%d): saved fit kept.\n', ...
                method_name, regions{region_n}, region_n, numel(regions));
            drawnow;
            continue
        end
        fit_options = mixture_options;
        fit_options.progress_label = sprintf('%s %s (%d/%d)', ...
            method_name, regions{region_n}, region_n, numel(regions));
        later_work = 0;
        for later_n = (region_n + 1):numel(regions)
            later_work = later_work + search_plans{later_n}.work_units;
        end
        fit_options.other_regions_work_units = later_work;
        region_results(region_n).fit = fit_gaussian_center_mixture( ...
            region_centers{region_n}, frequency_limits, fit_options);
        partial_result = struct();
        partial_result.by_region = region_results;
        partial_result.complete = false;
        region_gaussian_centers.methods.(method_name) = partial_result;
        save(output_file, 'region_gaussian_centers', '-v7.3');
        fprintf('Saved %s through %s %s.\n', output_file, ...
            method_name, regions{region_n});
        drawnow;
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
    method_result.complete = true;
    region_gaussian_centers.methods.(method_name) = method_result;
    save(output_file, 'region_gaussian_centers', '-v7.3');
    fprintf('Saved completed %s mixtures to %s.\n', method_name, output_file);
end
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
restarts = nan(n_regions, 1);
fit_seconds = nan(n_regions, 1);
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
    restarts(region_n) = this_region.fit.replicates_used;
    these_seconds = this_region.fit.component_fit_seconds;
    if isempty(these_seconds)
        fit_seconds(region_n) = NaN;
    else
        fit_seconds(region_n) = sum(these_seconds, 'omitnan');
    end
    status{region_n} = this_region.fit.status;
    model_table = this_region.fit.model_selection;
    if isfinite(best_n(region_n)) && ~isempty(model_table)
        best_row = model_table.NComponents == best_n(region_n);
        selected_bic(region_n) = model_table.BIC(best_row);
        selected_aic(region_n) = model_table.AIC(best_row);
    end
end
summary = table(names, n_peaks, n_channels, n_sessions, best_n, ...
    restarts, fit_seconds, selected_bic, selected_aic, status, ...
    'VariableNames', {'Region', 'NCenters', 'NChannels', 'NSessions', ...
    'BestNComponents', 'Restarts', 'FitSeconds', 'SelectedBIC', ...
    'SelectedAIC', 'Status'});
end

function print_mixture_plan(method_name, regions, search_plans, options)
options = gaussian_mixture_options(options);
fprintf(['\n%s mixture plan. At most %d Gaussians are tested. ' ...
    'Every count from 1 through that cap is fit.\n'], ...
    method_name, options.max_components);
if options.large_sample_replicates < options.replicates
    fprintf(['A region with more than %d centers uses %d restarts ' ...
        'instead of %d.\n'], options.large_sample_centers, ...
        options.large_sample_replicates, options.replicates);
else
    fprintf('Every region uses %d restarts.\n', options.replicates);
end
total_work = 0;
for region_n = 1:numel(regions)
    plan = search_plans{region_n};
    total_work = total_work + plan.work_units;
    if isfield(plan, 'already_saved') && plan.already_saved
        fprintf('  %s: saved fit kept.\n', regions{region_n});
        continue
    end
    if plan.restarts_reduced
        restart_text = sprintf('%d restarts, reduced from %d', ...
            plan.replicates, plan.replicates_requested);
    else
        restart_text = sprintf('%d restarts', plan.replicates);
    end
    if plan.max_components_tested < 1
        component_text = 'no fit';
    else
        component_text = sprintf('components 1:%d', ...
            plan.max_components_tested);
    end
    fprintf('  %s: %d centers, %s, %s. %s.\n', ...
        regions{region_n}, plan.n_centers, component_text, ...
        restart_text, plan.runtime_note);
end
if total_work <= 0
    kept = false(numel(search_plans), 1);
    for region_n = 1:numel(search_plans)
        plan = search_plans{region_n};
        kept(region_n) = isfield(plan, 'already_saved') && plan.already_saved;
    end
    if all(kept)
        fprintf('Every listed region already has a saved fit.\n');
    else
        fprintf('No region has enough centers for a mixture.\n');
    end
    drawnow;
    return
end
reference = search_plans{1}.reference_work_units;
if isempty(reference) || reference <= 0
    scale = NaN;
else
    scale = total_work / reference;
end
fprintf(['Combined work is %.1f times one mid-LPAG-sized search ' ...
    '(300 centers, %d restarts, up to %d Gaussians).\n'], ...
    scale, options.replicates, options.max_components);
if isfinite(scale) && scale > 10
    fprintf(['This is the full restart search and can take hours. ' ...
        'The first component will turn that into a clock-time estimate. ' ...
        'A projection past %s is reported while the fit is running. ' ...
        'Each finished region is saved before the next one starts.\n'], ...
        format_search_duration(options.slow_run_seconds));
else
    fprintf(['The first component will turn this into a clock-time ' ...
        'estimate.\n']);
end
drawnow;
end

function ready = saved_mixture_has_peaks(result)
ready = isstruct(result) && isfield(result, 'peak_table') && ...
    istable(result.peak_table) && ...
    all(ismember({'Region', 'Method', 'CenterFrequency'}, ...
    result.peak_table.Properties.VariableNames)) && ...
    isfield(result, 'frequency') && isfield(result, 'regions') && ...
    isfield(result, 'methods') && isfield(result, 'mixture_options');
end

function complete = method_is_complete(result, method_name)
complete = false;
method_result = saved_method_result(result, method_name);
if isempty(method_result) || ~isfield(method_result, 'complete')
    return
end
complete = isequal(method_result.complete, true);
end

function prior = saved_method_regions(result, method_name)
prior = [];
method_result = saved_method_result(result, method_name);
if isempty(method_result) || ~isfield(method_result, 'by_region')
    return
end
prior = method_result.by_region;
end

function method_result = saved_method_result(result, method_name)
method_result = [];
if ~isstruct(result) || ~isfield(result, 'methods') || ...
        ~isstruct(result.methods) || ~isfield(result.methods, method_name)
    return
end
method_result = result.methods.(method_name);
end

function kept_fit = saved_fit_for_region(prior, region_name)
kept_fit = [];
if isempty(prior)
    return
end
for region_n = 1:numel(prior)
    if ~isfield(prior(region_n), 'name') || ...
            ~strcmp(prior(region_n).name, region_name)
        continue
    end
    if region_mixture_already_saved(prior(region_n))
        kept_fit = prior(region_n).fit;
        return
    end
end
end
