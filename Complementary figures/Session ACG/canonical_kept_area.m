function name = canonical_kept_area(label, areas_to_keep)
%CANONICAL_KEPT_AREA  Map a depth-table label onto the areas kept for ACGs.
%
% SupCol and SC are the superior colliculus. Matching ignores case and
% spaces. VLPAG does not match LPAG. A label that is not in areas_to_keep
% returns ''.

name = '';
text = strtrim(char(string(label)));
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
