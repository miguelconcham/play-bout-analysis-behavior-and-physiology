function text = session_acg_char(value)
%SESSION_ACG_CHAR  One trimmed character row from a saved name.

if iscell(value)
    if isempty(value)
        text = '';
        return
    end
    value = value{1};
end
if isempty(value)
    text = '';
    return
end
text = char(value);
text = strtrim(text(:).');
end
