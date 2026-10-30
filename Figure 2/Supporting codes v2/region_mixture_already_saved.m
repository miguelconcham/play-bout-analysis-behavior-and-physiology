function already_saved = region_mixture_already_saved(region_result)
%REGION_MIXTURE_ALREADY_SAVED  True when a region has a stored fit.
%
% An interrupted run leaves the unfinished region with an empty fit.
% A finished region has a fit struct whose status is non-empty text,
% including a sample that was too short to fit.

already_saved = false;
if ~isstruct(region_result) || numel(region_result) ~= 1 || ...
        ~isfield(region_result, 'fit') || ~isstruct(region_result.fit) || ...
        numel(region_result.fit) ~= 1 || ~isfield(region_result.fit, 'status')
    return
end
status = region_result.fit.status;
if isstring(status)
    if numel(status) ~= 1
        return
    end
    status = char(status);
end
already_saved = ischar(status) && ~isempty(strtrim(status));
end
