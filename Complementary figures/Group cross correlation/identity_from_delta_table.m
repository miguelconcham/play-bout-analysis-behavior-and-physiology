function identity = identity_from_delta_table(all_neurons, alpha, non_entrained_p)
%IDENTITY_FROM_DELTA_TABLE  Lock type from EntireSession PPC and angle.
%
% all_neurons is the table in delta_all_neurons_v2.mat, or a struct with
% session, cluster_id, area, and EntireSession.PPCPval and PreferedAngle.
% The split is delta_lock_group: locked cells are trough or peak, and a
% PPC p-value at or above non_entrained_p is unlocked.

if is_table(all_neurons)
    if ~ismember('EntireSession', all_neurons.Properties.VariableNames)
        error('group_ccg:identity', ...
            'The delta table has no EntireSession columns.');
    end
    n = height(all_neurons);
    session = cell(n, 1);
    area = cell(n, 1);
    for i = 1:n
        session{i} = table_text(all_neurons, 'session', i);
        area{i} = table_text(all_neurons, 'area', i);
    end
    neuron_id = table_cluster_id(all_neurons);
    ppc = double(all_neurons.EntireSession.PPCPval(:));
    angle = double(all_neurons.EntireSession.PreferedAngle(:));
else
    session = all_neurons.session(:);
    area = all_neurons.area(:);
    neuron_id = double(all_neurons.cluster_id(:));
    ppc = double(all_neurons.EntireSession.PPCPval(:));
    angle = double(all_neurons.EntireSession.PreferedAngle(:));
    n = numel(neuron_id);
    if numel(session) ~= n || numel(area) ~= n || numel(ppc) ~= n || numel(angle) ~= n
        error('group_ccg:identity', ...
            'The delta struct does not have one session, area, PPC, and angle per neuron.');
    end
    for i = 1:n
        session{i} = ccg_plain_text(session{i});
        area{i} = ccg_plain_text(area{i});
    end
end
if numel(ppc) ~= numel(neuron_id) || numel(angle) ~= numel(neuron_id)
    error('group_ccg:identity', ...
        'EntireSession PPC and preferred angle must have one value per neuron.');
end
lock_group = cell(numel(neuron_id), 1);
for i = 1:numel(neuron_id)
    lock_group{i} = delta_lock_group(ppc(i), angle(i), alpha, non_entrained_p);
end
identity.session = session;
identity.area = area;
identity.neuron_id = neuron_id;
identity.lock_group = lock_group;
end


function tf = is_table(value)
tf = false;
if exist('istable', 'builtin') == 5 || exist('istable', 'file') == 2
    tf = istable(value);
end
end
