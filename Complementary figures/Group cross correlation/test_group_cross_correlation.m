function test_group_cross_correlation()
%TEST_GROUP_CROSS_CORRELATION  Groups, pooled correlograms, and session means.
%
% These checks do not need the experiment files. Run this file by name
% after add_repo_paths.

this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_dir));
addpath(fullfile(repo_root, 'Complementary figures'));
addpath(fullfile(repo_root, 'Complementary figures', 'Session ACG'));
addpath(this_dir);

root = group_ccg_data_root();
assert(~isempty(strfind(root, 'experimentfs.home.bccn-berlin.de')), ...
    'The data root is not the shared experiment host.'); %#ok<STREMP>

names = group_ccg_names();
assert(numel(names.group_id) == 6, 'The lock-by-oscillation split is not six groups.');
assert(numel(names.pair_reference) == 30, 'A group was paired with itself or a pair is missing.');
assert(all(names.pair_reference ~= names.pair_target), ...
    'A group was cross-correlated with itself.');
[group_id, group_index] = assign_ccg_group('nonentrained', true);
assert(strcmp(group_id, 'unlocked_oscillatory') && group_index == 5, ...
    'nonentrained oscillatory did not become the unlocked oscillatory group.');
[group_id, group_index] = assign_ccg_group('non-entrained', 0);
assert(strcmp(group_id, 'unlocked_nonoscillatory') && group_index == 6, ...
    'A non-oscillatory unlocked neuron was given another group.');
[group_id, group_index] = assign_ccg_group('trough', true);
assert(strcmp(group_id, 'trough_oscillatory') && group_index == 1, ...
    'Trough oscillatory was not the first group.');
[group_id, group_index] = assign_ccg_group('', false);
assert(isempty(group_id) && group_index == 0, ...
    'An unclassified lock was given a group.');

settings = group_ccg_settings(struct());
assert(numel(settings.centers) == 200, 'A +/- 1 s axis at 10 ms does not have 200 bins.');
assert(abs(settings.centers(1) + 0.995) < 1e-12, 'The first bin is not centered at -995 ms.');
assert(max(abs(settings.centers + fliplr(settings.centers))) < 1e-12, ...
    'The lag axis is not symmetric.');
empty_type = group_ccg_settings(struct('session_type', ''));
assert(isempty(empty_type.session_type), ...
    'An empty session type was replaced by the Partners session.');
failed = false;
try
    group_ccg_settings(struct('acg_lag', 0.015, 'acg_bin', 0.01));
catch
    failed = true;
end
assert(failed, 'A lag that is not a whole number of bins was accepted.');

labels = label_rows( ...
    {'A', 'A', 'A', 'A', 'A', 'A', 'A', 'B'}, ...
    1:8, ...
    [1 0 0 1 1 1 1 1], ...
    {'oscillatory', 'non-oscillatory', 'not classified', 'oscillatory', ...
    'oscillatory', 'oscillatory', 'oscillatory', 'non-oscillatory'});
identity = identity_rows( ...
    {'A', 'A', 'A', 'A', 'A', 'A'}, ...
    1:6, ...
    {'trough', 'peak', 'trough', 'nonentrained', 'trough', ''}, ...
    {'SC', 'VLPAG', 'SC', 'SupCol', 'DMPAG', 'SC'});
members = join_group_membership(labels, identity, {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'});
assert(members.n_labels == 8 && members.n_kept == 3, ...
    'The six-group filter did not keep the three eligible neurons.');
assert(members.n_unclassified == 1 && members.n_inconsistent == 1 ...
    && members.n_unmatched == 1 && members.n_no_lock == 1 ...
    && members.n_outside_area == 1, ...
    'A rejected label was not counted in its own reason.');
assert(isequal(members.group_id, ...
    {'trough_oscillatory'; 'peak_nonoscillatory'; 'unlocked_oscillatory'}), ...
    'The kept neurons are not in label order.');
assert(strcmp(members.area{3}, 'SC'), 'SupCol was not returned as SC.');

primary = identity_rows({'A'}, 1, {'trough'}, {'SC'});
extra = identity_rows({'A', 'A'}, [1 2], {'peak', 'nonentrained'}, {'DR', 'LPAG'});
merged = merge_group_identity(primary, extra);
assert(numel(merged.neuron_id) == 2, 'A neuron missing from the population was not added.');
assert(strcmp(merged.lock_group{1}, 'trough') && strcmp(merged.area{1}, 'SC'), ...
    'The delta table replaced the population lock or area.');
assert(merged.neuron_id(2) == 2 && strcmp(merged.lock_group{2}, 'nonentrained'), ...
    'The added neuron did not keep its delta lock.');

delta_only = identity_from_delta_table(struct( ...
    'session', {{'S', 'S'}}, ...
    'cluster_id', [7 8], ...
    'area', {{'LPAG', 'SC'}}, ...
    'EntireSession', struct('PPCPval', [0.001; 0.2], 'PreferedAngle', [pi; 0])), ...
    0.01, 0.1);
assert(strcmp(delta_only.lock_group{1}, 'trough'), ...
    'A locked cell at pi was not trough.');
assert(strcmp(delta_only.lock_group{2}, 'nonentrained'), ...
    'A PPC p-value of 0.2 was not unlocked.');

if exist('writetable', 'file') == 2 || exist('writetable', 'builtin') == 5
    label_file = fullfile(tempdir, 'oscillatory_neuron_labels_test.csv');
    writetable(table(1, {'B1 D1'}, 70, 0, {'non-oscillatory'}, ...
        'VariableNames', {'population_row', 'session', 'neuron_id', 'oscillatory', 'classification'}), ...
        label_file);
    read_back = read_oscillatory_labels(label_file);
    assert(strcmp(read_back.session{1}, 'B1 D1') && read_back.neuron_id(1) == 70 ...
        && read_back.oscillatory(1) == 0, ...
        'The label file did not keep the session and neuron id.');
end

edges = (-2:2) * 0.01;
centers = edges(1:end - 1) + 0.005;
one = group_pair_ccg(1, 1 + centers(4), [0 2], edges);
assert(isequal(one.counts, [0 0 0 1]), 'A positive lag was not placed in its bin.');
assert(one.n_reference == 1 && abs(one.conditional_rate(4) - 100) < 1e-9, ...
    'The conditional rate is not one count divided by bin width.');
assert(abs(one.target_rate - 0.5) < 1e-9, ...
    'The target rate is not one spike in a two-second window.');
assert(abs(one.excess_rate(4) - 99.5) < 1e-9, ...
    'Excess rate did not subtract the target rate.');
flipped = group_pair_ccg(1 + centers(4), 1, [0 2], edges);
assert(isequal(flipped.counts, fliplr(one.counts)), ...
    'Swapping the trains did not reverse the lag.');
on_edge = group_pair_ccg(1, [1 + edges(1), 1 + edges(end)], [0 2], edges);
assert(isequal(on_edge.counts, [1 0 0 1]), ...
    'The first and last edges were not included.');
outside = group_pair_ccg(1, 1 + edges(end) + 0.01, [0 2], edges);
assert(isequal(outside.counts, [0 0 0 0]) && outside.n_reference == 1, ...
    'A lag beyond the axis was counted.');
too_short = group_pair_ccg(0.02, 0.02, [0 0.04], edges);
assert(too_short.n_reference == 0 && all(isnan(too_short.excess_rate)), ...
    'A window equal to two lag margins was used.');

windows = [0 10; 20 30];
pooled = group_pair_ccg([5; 25], [5 + centers(4); 25 + centers(2)], windows, edges);
assert(isequal(pooled.counts, [0 1 0 1]), ...
    'Counts from two windows were not added.');
assert(pooled.n_reference == 2 && pooled.n_windows == 2 && abs(pooled.duration - 20) < 1e-12, ...
    'The two windows did not keep their own durations.');
assert(abs(pooled.target_rate - 0.1) < 1e-9 && abs(pooled.excess_rate(4) - 49.9) < 1e-9 ...
    && abs(pooled.excess_rate(1) + 0.1) < 1e-9, ...
    'The expected rate was not weighted by the reference spikes in each window.');
assert(abs(pooled.rate_ratio(4) - 500) < 1e-6, ...
    'The rate ratio is not the conditional rate divided by the target rate.');
overlap = group_pair_ccg(0.6, 0.6, [0 1; 0.5 1.5], edges);
assert(overlap.n_reference == 1 && overlap.n_windows == 1, ...
    'Overlapping windows counted the same spike twice.');

spikes = cell(1, 6);
spikes{1} = 1;
spikes{2} = 1 + centers(4);
for g = 3:6
    spikes{g} = [];
end
pairs = session_group_ccgs(spikes, [1 1 0 0 0 0], [0 2], edges, 1);
assert(pairs(1).usable && strcmp(pairs(1).reference, 'trough_oscillatory') ...
    && strcmp(pairs(1).target, 'trough_nonoscillatory'), ...
    'The first directed pair is not trough oscillatory against trough non-oscillatory.');
assert(isequal(pairs(1).counts, [0 0 0 1]), ...
    'The session pair did not use the pooled spike times.');
assert(~pairs(2).usable, 'A pair with an empty target group was marked usable.');
strict = session_group_ccgs(spikes, [1 1 0 0 0 0], [0 2], edges, 2);
assert(~strict(1).usable && strict(1).n_reference == 1, ...
    'One reference spike passed a minimum of two.');

centers = [-0.03 -0.01 0.01 0.03];
first = session_with_excess('A', centers, [-1 -1], true);
second = session_with_excess('B', centers, [-3 -3], true);
summary = pool_group_ccg_sessions([first, second], 0.02);
assert(abs(summary.central_mean(1) + 2) < 1e-12, ...
    'The session mean was not the equal-weight mean of the two sessions.');
assert(isequal(summary.mean_excess(1, 2:3), [-2 -2]), ...
    'Bins inside the central window were not averaged across sessions.');
assert(summary.n_sessions(1) == 2 && summary.n_sessions(2) == 0, ...
    'An unusable pair was included in the session count.');
assert(abs(summary.holm_p(1) - summary.raw_p(1)) < 1e-12 && isfinite(summary.raw_p(1)), ...
    'A single finite p-value was changed by the Holm-Sidak correction.');
assert(all(isnan(summary.holm_p(2:end))), ...
    'A pair without two sessions received a p-value.');

adjusted = holm_sidak_adjust([0.01 0.04]);
assert(abs(adjusted(1) - (1 - 0.99 ^ 2)) < 1e-12, ...
    'The smaller p-value was not adjusted for both tests.');
assert(abs(adjusted(2) - 0.04) < 1e-12, ...
    'The larger p-value fell below its raw value.');
with_missing = holm_sidak_adjust([NaN 0.2]);
assert(isnan(with_missing(1)) && abs(with_missing(2) - 0.2) < 1e-12, ...
    'A missing p-value entered the family.');

disp('Group cross-correlation checks passed.')
end


function labels = label_rows(session, neuron_id, oscillatory, classification)
labels.session = session(:);
labels.neuron_id = neuron_id(:);
labels.oscillatory = oscillatory(:);
labels.classification = classification(:);
end


function identity = identity_rows(session, neuron_id, lock_group, area)
identity.session = session(:);
identity.neuron_id = neuron_id(:);
identity.lock_group = lock_group(:);
identity.area = area(:);
end


function session = session_with_excess(animal, centers, central_values, usable)
names = group_ccg_names();
n_bins = numel(centers);
n_pairs = numel(names.pair_reference);
excess = zeros(1, n_bins);
inside = abs(centers) <= 0.02;
excess(inside) = central_values;
pair = struct( ...
    'reference', '', ...
    'target', '', ...
    'usable', false, ...
    'excess_rate', nan(1, n_bins), ...
    'rate_ratio', nan(1, n_bins), ...
    'n_neurons_reference', 0, ...
    'n_neurons_target', 0);
session.animal_code = animal;
session.centers = centers;
session.version = 1;
session.pairs = repmat(pair, n_pairs, 1);
for k = 1:n_pairs
    session.pairs(k).reference = names.group_id{names.pair_reference(k)};
    session.pairs(k).target = names.group_id{names.pair_target(k)};
    session.pairs(k).usable = k == 1 && usable;
    if session.pairs(k).usable
        session.pairs(k).excess_rate = excess;
        session.pairs(k).rate_ratio = excess + 1;
        session.pairs(k).n_neurons_reference = 2;
        session.pairs(k).n_neurons_target = 4;
    end
end
end
