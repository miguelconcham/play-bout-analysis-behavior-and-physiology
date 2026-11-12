function text = table_text(cluster_table, column_name, row)
%TABLE_TEXT  One table cell as a character vector.

if ~ismember(column_name, cluster_table.Properties.VariableNames)
    text = '';
    return
end
value = cluster_table.(column_name)(row, :);
if iscell(value)
    value = value{1};
    if iscell(value)
        value = value{1};
    end
end
text = char(string(value));
end
