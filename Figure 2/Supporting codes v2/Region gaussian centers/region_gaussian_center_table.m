function [peak_table, regions] = region_gaussian_center_table( ...
        probe_spectra, regions, spectral_method, session_index, recording_index)
%REGION_GAUSSIAN_CENTER_TABLE  FOOOF bump centers on one recording.
%
% regions is the requested list after SupCol/SC canonicalization, in the
% same order. peak_table uses those names.
%
% Each finite Gaussian center is one power bump. gaussian_params columns
% are center frequency (Hz), height in log10 power, and standard deviation
% (Hz). A channel whose area label names more than one requested region
% contributes the same bumps to each of those regions. VLPAG does not
% match LPAG, and DLPAG does not match LPAG. SupCol and SC are one region.

if nargin < 4
    session_index = [];
end
if nargin < 5 || isempty(recording_index)
    recording_index = NaN;
end

regions = normalize_region_list(regions);
method_names = normalize_method_list(spectral_method);
peak_table = empty_peak_table();

if isempty(session_index)
    sessions_in_file = 1:probe_spectra.n_sessions;
else
    sessions_in_file = session_index(:).';
end
if any(sessions_in_file < 1) || ...
        any(sessions_in_file > probe_spectra.n_sessions)
    error('regionGaussianCenters:session', ...
        '%s has %d partner sessions.', probe_spectra.animal_name, ...
        probe_spectra.n_sessions);
end
if isempty(sessions_in_file)
    return
end

channel_table = probe_spectra.channel_table;
n_channels = height(channel_table);
area_labels = cellstr(string(channel_table.Area));
pieces = cell(numel(sessions_in_file) * numel(method_names), 1);
piece_n = 0;

for session_n = sessions_in_file
    this_session = probe_spectra.sessions(session_n);
    for method_n = 1:numel(method_names)
        method_name = method_names{method_n};
        fitted = this_session.(method_field_name(method_name));
        gaussians = gaussian_cells(fitted, n_channels, ...
            probe_spectra.animal_name, session_n, method_name);
        piece_n = piece_n + 1;
        pieces{piece_n} = peaks_for_session( ...
            gaussians, area_labels, channel_table, regions, ...
            probe_spectra.animal_name, recording_index, session_n, ...
            method_name);
    end
end

pieces = pieces(1:piece_n);
pieces = pieces(~cellfun(@isempty, pieces));
if ~isempty(pieces)
    peak_table = vertcat(pieces{:});
end
end

function names = normalize_region_list(regions)
if ischar(regions) || isstring(regions)
    requested = cellstr(regions);
else
    requested = cellstr(string(regions));
end
requested = requested(:).';
if isempty(requested) || any(cellfun(@(name) isempty(strtrim(name)), requested))
    error('regionGaussianCenters:regions', ...
        'regions_to_analyze must name at least one region.');
end
names = cell(1, numel(requested));
for region_n = 1:numel(requested)
    names{region_n} = canonical_analysis_region(requested{region_n});
end
if numel(unique(names, 'stable')) ~= numel(names)
    error('regionGaussianCenters:regions', ...
        'regions_to_analyze lists the same region twice.');
end
end

function name = canonical_analysis_region(token)
token = strtrim(char(string(token)));
known = {'DR', 'VLPAG', 'LPAG', 'DLPAG', 'DMPAG', 'PAG'};
if strcmpi(token, 'SupCol') || strcmpi(token, 'SC')
    name = 'SupCol';
    return
end
match = strcmp(known, token);
if ~any(match)
    error('regionGaussianCenters:regions', ...
        ['Unknown region "%s". Use DR, VLPAG, LPAG, DLPAG, SupCol, ' ...
        'DMPAG, or PAG.'], token);
end
name = known{match};
end

function method_names = normalize_method_list(spectral_method)
if ischar(spectral_method) || isstring(spectral_method)
    requested = cellstr(spectral_method);
else
    requested = cellstr(string(spectral_method));
end
requested = requested(:).';
method_names = cell(1, numel(requested));
for method_n = 1:numel(requested)
    token = strtrim(requested{method_n});
    if strcmpi(token, 'Fourier')
        method_names{method_n} = 'Fourier';
    elseif strcmpi(token, 'Welch')
        method_names{method_n} = 'Welch';
    else
        error('regionGaussianCenters:method', ...
            'spectral_method must be Fourier, Welch, or both.');
    end
end
if isempty(method_names) || numel(unique(method_names, 'stable')) ~= numel(method_names)
    error('regionGaussianCenters:method', ...
        'spectral_method must list Fourier, Welch, or both, once each.');
end
end

function field_name = method_field_name(method_name)
if strcmp(method_name, 'Fourier')
    field_name = 'fourier';
else
    field_name = 'welch';
end
end

function gaussians = gaussian_cells(fitted, n_channels, animal_name, ...
        session_n, method_name)
if isfield(fitted, 'gaussian_params') && ~isempty(fitted.gaussian_params)
    gaussians = fitted.gaussian_params;
elseif isfield(fitted, 'peak_params') && ~isempty(fitted.peak_params)
    gaussians = fitted.peak_params;
    for channel_n = 1:numel(gaussians)
        if ~isempty(gaussians{channel_n})
            gaussians{channel_n}(:, 3) = gaussians{channel_n}(:, 3) / 2;
        end
    end
else
    error('regionGaussianCenters:gaussians', ...
        '%s session %d has no %s gaussian_params.', ...
        animal_name, session_n, method_name);
end
if ~iscell(gaussians) || numel(gaussians) ~= n_channels
    error('regionGaussianCenters:gaussians', ...
        ['%s session %d %s gaussian_params has %d entries for ' ...
        '%d channels.'], animal_name, session_n, method_name, ...
        numel(gaussians), n_channels);
end
end

function peak_table = peaks_for_session(gaussians, area_labels, ...
        channel_table, regions, animal_name, recording_index, ...
        session_n, method_name)
pieces = cell(numel(gaussians), 1);
for channel_n = 1:numel(gaussians)
    region_names = regions_in_area_label(area_labels{channel_n}, regions);
    if isempty(region_names)
        continue
    end
    parameters = gaussian_matrix(gaussians{channel_n});
    if isempty(parameters)
        continue
    end
    pieces{channel_n} = rows_for_channel(parameters, region_names, ...
        animal_name, recording_index, session_n, method_name, ...
        channel_table.Channel(channel_n), ...
        channel_table.ProbeNumber(channel_n), ...
        area_labels{channel_n}, channel_table.Yum(channel_n));
end
pieces = pieces(~cellfun(@isempty, pieces));
if isempty(pieces)
    peak_table = empty_peak_table();
else
    peak_table = vertcat(pieces{:});
end
end

function parameters = gaussian_matrix(parameters)
if isempty(parameters)
    parameters = zeros(0, 3);
    return
end
if size(parameters, 2) < 3
    error('regionGaussianCenters:gaussians', ...
        'Gaussian parameters must have center, height, and SD columns.');
end
parameters = parameters(:, 1:3);
parameters = parameters(isfinite(parameters(:, 1)), :);
end

function peak_table = rows_for_channel(parameters, region_names, ...
        animal_name, recording_index, session_n, method_name, ...
        channel, probe_number, area_label, y_um)
n_peaks = size(parameters, 1);
n_regions = numel(region_names);
n_rows = n_peaks * n_regions;
region_grid = repmat(region_names(:).', n_peaks, 1);
region_column = reshape(region_grid, n_rows, 1);
center = repmat(parameters(:, 1), n_regions, 1);
height = repmat(parameters(:, 2), n_regions, 1);
sd = repmat(parameters(:, 3), n_regions, 1);
peak_table = table( ...
    region_column, ...
    repmat({animal_name}, n_rows, 1), ...
    repmat(recording_index, n_rows, 1), ...
    repmat(session_n, n_rows, 1), ...
    repmat({method_name}, n_rows, 1), ...
    repmat(channel, n_rows, 1), ...
    repmat(probe_number, n_rows, 1), ...
    repmat({area_label}, n_rows, 1), ...
    repmat(y_um, n_rows, 1), ...
    center, height, sd, ...
    'VariableNames', {'Region', 'Animal', 'RecordingIndex', ...
    'SessionIndex', 'Method', 'Channel', 'ProbeNumber', 'Area', 'Yum', ...
    'CenterFrequency', 'Height', 'SD'});
end

function names = regions_in_area_label(label, requested_names)
tokens = strtrim(strsplit(char(string(label)), ','));
tokens = tokens(~cellfun('isempty', tokens));
names = cell(1, numel(tokens));
n_names = 0;
for token_n = 1:numel(tokens)
    canon = canonical_region_token(tokens{token_n});
    if ~isempty(canon) && any(strcmp(requested_names, canon))
        n_names = n_names + 1;
        names{n_names} = canon;
    end
end
names = unique(names(1:n_names), 'stable');
end

function name = canonical_region_token(token)
token = strtrim(char(string(token)));
if strcmpi(token, 'SupCol') || strcmpi(token, 'SC')
    name = 'SupCol';
    return
end
known = {'DR', 'VLPAG', 'LPAG', 'DLPAG', 'DMPAG', 'PAG'};
match = strcmp(known, token);
if any(match)
    name = known{match};
else
    name = '';
end
end

function peak_table = empty_peak_table()
peak_table = table(cell(0, 1), cell(0, 1), zeros(0, 1), zeros(0, 1), ...
    cell(0, 1), zeros(0, 1), zeros(0, 1), cell(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', {'Region', 'Animal', 'RecordingIndex', ...
    'SessionIndex', 'Method', 'Channel', 'ProbeNumber', 'Area', 'Yum', ...
    'CenterFrequency', 'Height', 'SD'});
end
