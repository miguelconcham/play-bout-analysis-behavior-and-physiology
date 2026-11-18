function members = join_group_membership(labels, identity, areas_to_keep)
%JOIN_GROUP_MEMBERSHIP  Neurons that land in one of the six groups.
%
% A row of the label table is kept when its classification agrees with the
% oscillatory flag, the same session and neuron id has a peak, trough, or
% unlocked label, and its area is one of areas_to_keep. 'not classified'
% rows stay out. SupCol is kept as SC.

if nargin < 3 || isempty(areas_to_keep)
    areas_to_keep = {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'};
end
n = numel(labels.neuron_id);
if numel(labels.session) ~= n || numel(labels.classification) ~= n ...
        || numel(labels.oscillatory) ~= n
    error('group_ccg:members', ...
        'The label table does not have one session, id, and class per row.');
end
names = group_ccg_names();
label_keys = neuron_keys(labels.session, labels.neuron_id);
[~, first_label] = unique(label_keys, 'stable');
if numel(first_label) ~= n
    error('group_ccg:members', ...
        'The label table repeats a session and neuron id.');
end
identity_keys = neuron_keys(identity.session, identity.neuron_id);
[~, first_identity] = unique(identity_keys, 'stable');
if numel(first_identity) ~= numel(identity_keys)
    error('group_ccg:members', ...
        'The lock-and-area table repeats a session and neuron id.');
end

session = cell(n, 1);
area = cell(n, 1);
lock_group = cell(n, 1);
group_id = cell(n, 1);
neuron_id = nan(n, 1);
group_index = nan(n, 1);
oscillatory = false(n, 1);
keep = false(n, 1);
n_unclassified = 0;
n_inconsistent = 0;
n_unmatched = 0;
n_no_lock = 0;
n_outside_area = 0;

for i = 1:n
    [is_oscillatory, state] = oscillation_state( ...
        labels.classification{i}, labels.oscillatory(i));
    if strcmp(state, 'unclassified')
        n_unclassified = n_unclassified + 1;
        continue
    end
    if strcmp(state, 'inconsistent')
        n_inconsistent = n_inconsistent + 1;
        continue
    end
    hit = find(strcmp(identity_keys, label_keys{i}), 1);
    if isempty(hit)
        n_unmatched = n_unmatched + 1;
        continue
    end
    this_lock = identity.lock_group{hit};
    [this_group, this_index] = assign_ccg_group(this_lock, is_oscillatory);
    if isempty(this_group)
        n_no_lock = n_no_lock + 1;
        continue
    end
    this_area = canonical_kept_area(identity.area{hit}, areas_to_keep);
    if isempty(this_area)
        n_outside_area = n_outside_area + 1;
        continue
    end
    keep(i) = true;
    session{i} = labels.session{i};
    neuron_id(i) = labels.neuron_id(i);
    area{i} = this_area;
    lock_group{i} = names.lock{this_index};
    group_id{i} = this_group;
    group_index(i) = this_index;
    oscillatory(i) = is_oscillatory;
end

members.session = session(keep);
members.neuron_id = neuron_id(keep);
members.area = area(keep);
members.lock_group = lock_group(keep);
members.group_id = group_id(keep);
members.group_index = group_index(keep);
members.oscillatory = oscillatory(keep);
members.n_labels = n;
members.n_kept = sum(keep);
members.n_unclassified = n_unclassified;
members.n_inconsistent = n_inconsistent;
members.n_unmatched = n_unmatched;
members.n_no_lock = n_no_lock;
members.n_outside_area = n_outside_area;
end


function keys = neuron_keys(session, neuron_id)
keys = cell(numel(neuron_id), 1);
for i = 1:numel(neuron_id)
    keys{i} = [ccg_plain_text(session{i}), '|', sprintf('%.0f', double(neuron_id(i)))];
end
end


function [is_oscillatory, state] = oscillation_state(classification, flag)
is_oscillatory = false;
name = lower(ccg_plain_text(classification));
flag = double(flag);
if strcmp(name, 'oscillatory')
    if isfinite(flag) && flag == 1
        is_oscillatory = true;
        state = 'oscillatory';
    else
        state = 'inconsistent';
    end
elseif strcmp(name, 'non-oscillatory')
    if isfinite(flag) && flag == 0
        state = 'non-oscillatory';
    else
        state = 'inconsistent';
    end
else
    state = 'unclassified';
end
end
