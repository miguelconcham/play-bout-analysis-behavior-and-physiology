function identity = identity_from_population(session, area, cluster_id, neuron_type)
%IDENTITY_FROM_POPULATION  Lock type and area in the population-logit order.
%
% session, area, cluster_id, and neuron_type are the fields saved in
% acg_population_logit_prep.mat. neuron_type is 'trough', 'peak',
% 'nonentrained', or '' when that neuron was not classified.

session = session(:);
area = area(:);
cluster_id = double(cluster_id(:));
neuron_type = neuron_type(:);
n = numel(neuron_type);
if numel(session) ~= n || numel(area) ~= n || numel(cluster_id) ~= n
    error('group_ccg:identity', ...
        'Session, area, cluster id, and lock type must have one entry per neuron.');
end
identity.session = cell(n, 1);
identity.area = cell(n, 1);
identity.lock_group = cell(n, 1);
for i = 1:n
    identity.session{i} = ccg_plain_text(session{i});
    identity.area{i} = ccg_plain_text(area{i});
    identity.lock_group{i} = ccg_plain_text(neuron_type{i});
end
identity.neuron_id = cluster_id;
end
