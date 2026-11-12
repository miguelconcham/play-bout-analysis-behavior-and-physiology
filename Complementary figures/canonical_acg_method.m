function method = canonical_acg_method(method)
%CANONICAL_ACG_METHOD  'eligible' or 'overlap'.
%
% 'eligible' keeps one observation per spike that can see the whole lag
% bin. 'overlap' is pair counts divided by T-|lag|. An empty value is
% 'eligible'.

if nargin < 1 || isempty(method)
    method = 'eligible';
end
method = lower(strtrim(char(method)));
if ~any(strcmp(method, {'eligible', 'overlap'}))
    error('canonical_acg_method:name', ...
        'acg_method must be ''eligible'' or ''overlap''.');
end
end
