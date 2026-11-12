function name = canonical_kept_area(label, areas_to_keep)
%CANONICAL_KEPT_AREA  Map a depth-table label onto the areas used here.
%
% Superior colliculus is stored as SupCol and returned as SC. Matching
% ignores case and spaces. VLPAG is not treated as LPAG. A label that is
% not in areas_to_keep, including a bare PAG, returns ''.

if nargin < 2 || isempty(areas_to_keep)
    areas_to_keep = {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'};
end
name = '';
while iscell(label)
    if isempty(label)
        label = '';
        break
    end
    label = label{1};
end
if isstring(label)
    if numel(label) ~= 1
        return
    end
    label = char(label);
end
if isnumeric(label) || isempty(label)
    return
end
text = strtrim(char(label));
if isempty(text)
    return
end
tokens = regexp(text, '[,;/]+', 'split');
for k = 1:numel(tokens)
    token = upper(regexprep(strtrim(tokens{k}), '\s+', ''));
    if any(strcmp(token, {'SC', 'SUPCOL', 'SUPERIORCOLLICULUS'}))
        token = 'SC';
    end
    hit = strcmpi(areas_to_keep, token);
    if any(hit)
        name = areas_to_keep{find(hit, 1)};
        return
    end
end
end
