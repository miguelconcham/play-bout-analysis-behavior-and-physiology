function text = ccg_plain_text(value)
%CCG_PLAIN_TEXT  One cell, string, or char value as a trimmed character vector.

while iscell(value)
    if isempty(value)
        text = '';
        return
    end
    value = value{1};
end
if is_string(value)
    if numel(value) ~= 1
        text = '';
        return
    end
    value = char(value);
end
if ischar(value)
    text = strtrim(reshape(value, 1, []));
    return
end
if isempty(value) || isnumeric(value)
    text = '';
    return
end
text = strtrim(char(value));
end


function tf = is_string(value)
tf = false;
if exist('isstring', 'builtin') == 5 || exist('isstring', 'file') == 2
    tf = isstring(value);
end
end
