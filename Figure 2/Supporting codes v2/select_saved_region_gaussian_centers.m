function [selected, missing] = select_saved_region_gaussian_centers( ...
        result, regions_to_plot, spectral_method)
%SELECT_SAVED_REGION_GAUSSIAN_CENTERS  Regions to draw from a saved fit.
%
% regions_to_plot names the rows of the figure. An empty list uses every
% region stored in the file, in that saved order. A region with an empty
% fit, or a name that is not in the file, is returned in missing and left
% out of selected. SupCol and SC are the same region. A repeated region
% is kept once, at its first request.

method_name = canonical_saved_method_name(spectral_method);
missing = struct('name', {}, 'reason', {});
selected = empty_selected_region();
requested = requested_region_tokens(regions_to_plot, result);
if isempty(requested)
    return
end

if ~isstruct(result) || ~isfield(result, 'methods') || ...
        ~isstruct(result.methods) || ~isfield(result.methods, method_name)
    missing = missing_for_tokens(requested, sprintf( ...
        '%s results were not saved', method_name));
    return
end

method_result = result.methods.(method_name);
prior = [];
if isstruct(method_result) && isfield(method_result, 'by_region')
    prior = method_result.by_region;
end

seen = {};
for token_n = 1:numel(requested)
    [region_name, known] = canonical_saved_region_name(requested{token_n});
    if ~known
        missing(end + 1, 1) = struct( ... %#ok<AGROW>
            'name', strtrim(char(string(requested{token_n}))), ...
            'reason', 'unknown region');
        continue
    end
    if any(strcmp(seen, region_name))
        continue
    end
    seen{end + 1} = region_name; %#ok<AGROW>
    match = find_saved_region(prior, region_name);
    if isempty(match)
        missing(end + 1, 1) = struct( ... %#ok<AGROW>
            'name', region_name, ...
            'reason', 'not in the saved regions');
        continue
    end
    if ~region_mixture_already_saved(match)
        missing(end + 1, 1) = struct( ... %#ok<AGROW>
            'name', region_name, ...
            'reason', 'fit not saved yet');
        continue
    end
    if isempty(selected)
        selected = match;
    else
        selected(end + 1, 1) = match; %#ok<AGROW>
    end
end
end

function tokens = requested_region_tokens(regions_to_plot, result)
if nargin < 1 || isempty(regions_to_plot)
    if isstruct(result) && isfield(result, 'regions') && ...
            ~isempty(result.regions)
        tokens = cellstr(string(result.regions));
    else
        tokens = {};
    end
elseif ischar(regions_to_plot) || isstring(regions_to_plot)
    tokens = cellstr(regions_to_plot);
else
    tokens = cellstr(string(regions_to_plot));
end
tokens = tokens(:).';
keep = true(1, numel(tokens));
for token_n = 1:numel(tokens)
    keep(token_n) = ~isempty(strtrim(char(string(tokens{token_n}))));
end
tokens = tokens(keep);
end

function missing = missing_for_tokens(tokens, reason)
missing = struct('name', {}, 'reason', {});
for token_n = 1:numel(tokens)
    [region_name, known] = canonical_saved_region_name(tokens{token_n});
    if ~known
        region_name = strtrim(char(string(tokens{token_n})));
        this_reason = 'unknown region';
    else
        this_reason = reason;
    end
    if any(strcmp({missing.name}, region_name))
        continue
    end
    missing(end + 1, 1) = struct( ... %#ok<AGROW>
        'name', region_name, ...
        'reason', this_reason);
end
end

function [name, known] = canonical_saved_region_name(token)
token = strtrim(char(string(token)));
known_names = {'DR', 'VLPAG', 'LPAG', 'DLPAG', 'DMPAG', 'PAG'};
if strcmpi(token, 'SupCol') || strcmpi(token, 'SC')
    name = 'SupCol';
    known = true;
    return
end
match = strcmp(known_names, token);
known = any(match);
if known
    name = known_names{match};
else
    name = '';
end
end

function name = canonical_saved_method_name(spectral_method)
if ~(ischar(spectral_method) || isstring(spectral_method))
    error('regionGaussianCenters:method', ...
        'spectral_method must be Fourier or Welch.');
end
token = strtrim(char(string(spectral_method)));
if strcmpi(token, 'Fourier')
    name = 'Fourier';
elseif strcmpi(token, 'Welch')
    name = 'Welch';
else
    error('regionGaussianCenters:method', ...
        'spectral_method must be Fourier or Welch.');
end
end

function match = find_saved_region(prior, region_name)
match = [];
if isempty(prior)
    return
end
for region_n = 1:numel(prior)
    if isfield(prior(region_n), 'name') && ...
            strcmp(prior(region_n).name, region_name)
        match = prior(region_n);
        return
    end
end
end

function selected = empty_selected_region()
selected = struct( ...
    'name', {}, ...
    'n_peaks', {}, ...
    'n_channels', {}, ...
    'n_sessions', {}, ...
    'fit', {});
end
