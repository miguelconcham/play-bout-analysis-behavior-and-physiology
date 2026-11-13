function id_column = table_cluster_id(cluster_table)
%TABLE_CLUSTER_ID  Numeric cluster id from a cluster_id or id column.

if ismember('cluster_id', cluster_table.Properties.VariableNames)
    id_column = cluster_table.cluster_id;
elseif ismember('id', cluster_table.Properties.VariableNames)
    id_column = cluster_table.id;
else
    error('table_cluster_id:id', 'The table has neither cluster_id nor id.');
end
id_column = double(id_column);
end
