function matches = session_name_matches(session_column, animal_code)
%SESSION_NAME_MATCHES  Rows of a session column equal to one animal folder.

animal_code = char(string(animal_code));
if iscell(session_column)
    matches = strcmp(session_column, animal_code);
    empty = cellfun(@isempty, session_column);
    matches(empty) = false;
    nested = cellfun(@iscell, session_column);
    if any(nested)
        matches(nested) = cellfun(@(value) ~isempty(value) && strcmp(value{1}, animal_code), ...
            session_column(nested));
    end
else
    matches = strcmp(string(session_column), animal_code);
end
matches = matches(:);
end
