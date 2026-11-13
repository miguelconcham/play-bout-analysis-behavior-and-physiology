function test_neuron_pair_cross_correlation()
%TEST_NEURON_PAIR_CROSS_CORRELATION  Pair matrices from crosscorrelogram.
%
% These checks do not need the experiment files. Run this file by name
% after add_repo_paths.

this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_dir));
addpath(fullfile(repo_root, 'Custom functions'));
addpath(fullfile(repo_root, 'Complementary figures'));
addpath(this_dir);

edges = (-50:50) * 0.001;
both = crosscorrelogram_bin_counts([1; 1.2], [1.015; 1.215], edges);
assert(sum(both) == 2, 'A reference spike was dropped by crosscorrelogram.');
[~, peak] = max(both);
assert(edges(peak) <= 0.015 && edges(peak + 1) > 0.015, ...
    'The 15 ms lag was not counted in the 1 ms bin that contains it.');
one = crosscorrelogram_bin_counts(1, 1.015, edges);
assert(sum(one) == 1 && one(peak) == 1, ...
    'The single reference spike did not keep its 1 ms bin.');

wide = (-1000:1000) * 0.001;
trains = cell(1, 6);
ids = cell(1, 6);
for g = 1:6
    trains{g} = {};
    ids{g} = zeros(0, 1);
end
trains{1} = {5, 5.2};
ids{1} = [10; 11];
trains{2} = {5.015};
ids{2} = 20;
combinations = session_neuron_pair_ccgs(trains, ids, [0 10], wide, 1);
assert(combinations(1).n_pairs == 2, ...
    'The two reference neurons did not each meet the target neuron.');
assert(combinations(1).reference_id(1) == 10 && combinations(1).target_id(2) == 20, ...
    'The pair matrix did not keep the neuron ids.');
assert(combinations(2).n_pairs == 0, ...
    'An empty target group still produced pairs.');
assert(sum(combinations(1).counts(1, :)) == 1, ...
    'The first neuron pair did not keep one lag from crosscorrelogram.');
second_lag = 5.015 - 5.2;
[~, second_peak] = max(combinations(1).counts(2, :));
second_center = (wide(second_peak) + wide(second_peak + 1)) / 2;
assert(abs(second_center - second_lag) <= 0.001, ...
    'The second pair lag was not placed in its 1 ms bin.');

centers = wide(1:end - 1) + 0.0005;
first = pair_session('A', centers, wide, 2, 100, 10);
second = pair_session('B', centers, wide, 0.2, 50, 2);
summary = pool_neuron_pair_sessions([first, second], 0.05);
assert(summary.n_pairs(1) == 2, 'Both sessions'' pairs were not kept.');
assert(abs(summary.central_mean(1) - 1) < 1e-9, ...
    'Pairs were averaged before each one was divided by its own target rate.');
assert(max(abs(summary.ci_low(1, :) - summary.mean_excess(1, :))) < 1e-8, ...
    'The confidence interval was not the spread of the normalized session means.');
assert(size(summary.pair_excess{1}, 1) == 2, ...
    'The neuron-pair image does not have one row per pair.');

disp('Neuron-pair cross-correlation checks passed.')
end


function session = pair_session(animal, centers, edges, count_value, n_reference, target_rate)
trains = cell(1, 6);
ids = cell(1, 6);
for g = 1:6
    trains{g} = {};
    ids{g} = zeros(0, 1);
end
session.animal_code = animal;
session.centers = centers;
session.version = 2;
session.combinations = session_neuron_pair_ccgs(trains, ids, [0 10], edges, 1);
n_bins = numel(centers);
session.combinations(1).counts = count_value * ones(1, n_bins);
session.combinations(1).n_reference = n_reference;
session.combinations(1).n_target = 1;
session.combinations(1).duration = 10;
session.combinations(1).target_rate = target_rate;
session.combinations(1).reference_id = 1;
session.combinations(1).target_id = 2;
session.combinations(1).usable = true;
session.combinations(1).n_pairs = 1;
session.combinations(1).n_usable = 1;
end
