function names = group_ccg_names()
%GROUP_CCG_NAMES  The six lock-by-oscillation groups and their pair order.
%
% Lock type is peak, trough, or unlocked. Oscillation is the label in
% oscillatory_neuron_labels.csv. A group is not paired with itself.
% pair_reference is the train whose spikes are the time origin.
% pair_target is the train counted around that origin. Lag is target time
% minus reference time, so a negative excess at a positive lag means the
% target fires less than its own rate after the reference group fires.

names.group_id = { ...
    'trough_oscillatory', 'trough_nonoscillatory', ...
    'peak_oscillatory', 'peak_nonoscillatory', ...
    'unlocked_oscillatory', 'unlocked_nonoscillatory'};
names.lock = {'trough', 'trough', 'peak', 'peak', 'unlocked', 'unlocked'};
names.oscillatory = [true, false, true, false, true, false];
names.title = { ...
    'Trough, oscillatory', 'Trough, non-oscillatory', ...
    'Peak, oscillatory', 'Peak, non-oscillatory', ...
    'Unlocked, oscillatory', 'Unlocked, non-oscillatory'};
names.color = [ ...
    0.10 0.35 0.75; ...
    0.55 0.75 0.90; ...
    0.75 0.15 0.15; ...
    0.90 0.55 0.50; ...
    0.20 0.20 0.20; ...
    0.62 0.62 0.62];
n_groups = numel(names.group_id);
pair_reference = zeros(n_groups * (n_groups - 1), 1);
pair_target = zeros(n_groups * (n_groups - 1), 1);
k = 0;
for i = 1:n_groups
    for j = 1:n_groups
        if i == j
            continue
        end
        k = k + 1;
        pair_reference(k) = i;
        pair_target(k) = j;
    end
end
names.pair_reference = pair_reference;
names.pair_target = pair_target;
end
