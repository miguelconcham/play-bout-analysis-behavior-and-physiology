function [heights, note] = region_bump_heights( ...
        region_result, peak_table, method_name)
%REGION_BUMP_HEIGHTS  FOOOF heights aligned with the stored bump centers.
%
% Heights come from the peak table saved by
% Estimate_play_session_region_gaussian_centers.m. They stay in the same
% order as fit.center_frequencies: finite centers only. heights is empty
% and note explains why when the column is missing or the centers do not
% match, so a height is never paired with the wrong frequency.

heights = zeros(0, 1);
note = 'Bump height was not stored';
if ~isstruct(region_result) || ~isfield(region_result, 'fit') || ...
        ~isstruct(region_result.fit) || ...
        ~isfield(region_result.fit, 'center_frequencies')
    note = 'No bump centers were stored';
    return
end
stored_centers = double(region_result.fit.center_frequencies(:));
if ~istable(peak_table) || ~all(ismember( ...
        {'Region', 'Method', 'CenterFrequency', 'Height'}, ...
        peak_table.Properties.VariableNames))
    return
end
if ~isfield(region_result, 'name') || isempty(region_result.name)
    note = 'The region name was not stored';
    return
end

method_name = char(string(method_name));
in_region = strcmp(peak_table.Region, region_result.name) & ...
    strcmp(peak_table.Method, method_name);
table_centers = double(peak_table.CenterFrequency(in_region));
table_heights = double(peak_table.Height(in_region));
finite_center = isfinite(table_centers);
table_centers = table_centers(finite_center);
table_heights = table_heights(finite_center);
if numel(table_centers) ~= numel(stored_centers) || ...
        any(abs(table_centers(:) - stored_centers(:)) > 1e-6)
    note = ['Stored centers do not match the peak table, so heights ' ...
        'were left out'];
    return
end
heights = table_heights(:);
note = '';
end
