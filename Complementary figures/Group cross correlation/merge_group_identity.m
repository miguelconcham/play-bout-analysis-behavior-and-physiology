function identity = merge_group_identity(primary, extra)
%MERGE_GROUP_IDENTITY  Keep the population row when a neuron is in both.
%
% extra supplies session and neuron ids that primary does not already
% contain. A neuron present in the population prep keeps that lock type
% and area, including an empty lock, so the delta table cannot reclassify
% a neuron the population figures left unclassified.

if isempty(primary) || ~isfield(primary, 'neuron_id') || isempty(primary.neuron_id)
    identity = extra;
    return
end
if isempty(extra) || ~isfield(extra, 'neuron_id') || isempty(extra.neuron_id)
    identity = primary;
    return
end
primary_keys = neuron_keys(primary.session, primary.neuron_id);
extra_keys = neuron_keys(extra.session, extra.neuron_id);
[~, first_primary] = unique(primary_keys, 'stable');
if numel(first_primary) ~= numel(primary_keys)
    error('group_ccg:identity', ...
        'The population identity has more than one row for a session and neuron id.');
end
add = ~ismember(extra_keys, primary_keys);
identity.session = [primary.session(:); extra.session(add)];
identity.area = [primary.area(:); extra.area(add)];
identity.neuron_id = [double(primary.neuron_id(:)); double(extra.neuron_id(add))];
identity.lock_group = [primary.lock_group(:); extra.lock_group(add)];
end


function keys = neuron_keys(session, neuron_id)
keys = cell(numel(neuron_id), 1);
for i = 1:numel(neuron_id)
    keys{i} = [ccg_plain_text(session{i}), '|', sprintf('%.0f', double(neuron_id(i)))];
end
end
