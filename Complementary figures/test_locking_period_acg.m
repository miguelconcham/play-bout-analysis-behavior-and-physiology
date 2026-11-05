function test_locking_period_acg()
%TEST_LOCKING_PERIOD_ACG  Checks that do not need the experiment files.
%
% Delta-lock labels, within-neuron locking quartiles, overlap-weighted
% pooling, and per-neuron ACG normalization. Run this file by name after
% add_repo_paths.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
    addpath(fullfile(fileparts(this_dir), 'Figure 2', 'Supporting codes', 'Delta episodes'));
end

assert(strcmp(delta_lock_group(0.001, pi, 0.01, 0.1), 'trough'), ...
    'A locked cell at pi was not trough-locked.');
assert(strcmp(delta_lock_group(0.001, -pi, 0.01, 0.1), 'trough'), ...
    'A locked cell at -pi was not trough-locked.');
assert(strcmp(delta_lock_group(0.001, 0, 0.01, 0.1), 'peak'), ...
    'A locked cell at 0 was not peak-locked.');
assert(strcmp(delta_lock_group(0.001, pi / 2, 0.01, 0.1), 'peak'), ...
    'The closed edge of [-pi/2, pi/2] was treated as a trough.');
assert(strcmp(delta_lock_group(0.001, -pi / 2, 0.01, 0.1), 'peak'), ...
    'The closed edge -pi/2 was treated as a trough.');
assert(strcmp(delta_lock_group(0.2, pi, 0.01, 0.1), 'nonentrained'), ...
    'A cell above the non-entrained threshold stayed locked.');
assert(strcmp(delta_lock_group(0.05, 0, 0.01, 0.1), ''), ...
    'A cell between 0.01 and 0.1 was assigned to a group.');
assert(strcmp(delta_lock_group(0.001, NaN, 0.01, 0.1), ''), ...
    'A locked cell with no preferred angle was assigned.');
assert(strcmp(delta_lock_group(NaN, 0, 0.01, 0.1), ''), ...
    'A missing p-value was assigned.');

ppc = (0.1:0.1:0.8)';
pooled_probe = pool_locking_period_acgs(zeros(8, 1), [0 1], ones(8, 1), ppc, false);
assert(isequal(pooled_probe.low_index, [true; true; false; false; false; false; false; false]), ...
    'The low 25% of locking was not the two weakest runs.');
assert(isequal(pooled_probe.high_index, [false; false; false; false; false; false; true; true]), ...
    'The high 25% of locking was not the two strongest runs.');

edges = [0.05 0.15];
counts = [90; 10];
run_length = [1; 0.3];
rates = [100; 50];
pooled = pool_locking_period_acgs(counts, edges, run_length, [0.1; 0.9], false);
assert(all(isnan(pooled.low)) && all(isnan(pooled.high)), ...
    'Two runs still produced a low/high split.');
assert(abs(pooled.all - (100 / 1.1)) < 1e-9, ...
    'The pooled ACG was not total counts over total overlap.');
assert(abs(pooled.all - mean(rates)) > 1, ...
    'The pooled ACG collapsed to the mean of the per-run rates.');
pooled_rates = pool_locking_period_acgs(rates, edges, run_length, [0.1; 0.9], true);
assert(abs(pooled_rates.all - pooled.all) < 1e-9, ...
    'Edge-corrected rows were not converted back to counts before pooling.');

n_bins = 1;
acg = zeros(9, n_bins);
acg(1:2, 1) = 10;
acg(7:8, 1) = 40;
acg(9, 1) = 1000;
run_length = ones(9, 1);
ppc = [ppc; NaN];
% One bin centered at 0, so overlap equals the run length.
pooled = pool_locking_period_acgs(acg, [-0.5 0.5], run_length, ppc, false);
assert(pooled.n_low == 2 && pooled.n_high == 2, ...
    'A run with no PPC entered a locking quartile.');
assert(abs(pooled.low - 10) < 1e-9, 'Low-quartile pool did not keep the weak runs.');
assert(abs(pooled.high - 40) < 1e-9, 'High-quartile pool did not keep the strong runs.');
assert(pooled.n_all == 9, 'The run with no PPC was dropped from the all-runs pool.');
assert(pooled.all > 40, 'The unlabeled run did not enter the all-runs pool.');

power = [1 2 3 4 40 50 60 70]';
power_pool = pool_locking_period_acgs(zeros(8, 1), [-0.5 0.5], ones(8, 1), power, false);
assert(isequal(power_pool.low_index, [true; true; false; false; false; false; false; false]), ...
    'The low 25% of delta power was not the two weakest runs.');
assert(isequal(power_pool.high_index, [false; false; false; false; false; false; true; true]), ...
    'The high 25% of delta power was not the two strongest runs.');

tied = pool_locking_period_acgs(ones(6, 1), [-0.5 0.5], ones(6, 1), ones(6, 1), false);
assert(all(isnan(tied.low)) && all(isnan(tied.high)) && tied.n_low == 0, ...
    'Identical PPC values still produced low and high tails.');
assert(abs(tied.all - 1) < 1e-9, 'Identical runs did not pool to their common rate.');

raw = [100 0 0; 1 1 NaN];
normalized = normalize_acg_rows(raw);
expected = [1 0 0; 0.5 0.5 NaN];
assert(max(abs(normalized(1, :) - expected(1, :))) < 1e-12, ...
    'A row was not scaled by its own pair count.');
assert(max(abs(normalized(2, 1:2) - expected(2, 1:2))) < 1e-12, ...
    'NaN bins were included in the normalizing sum.');
assert(isnan(normalized(2, 3)), 'A NaN lag bin became a number.');
mean_acg = mean(normalized, 1, 'omitnan');
assert(abs(mean_acg(1) - 0.75) < 1e-12, ...
    'The neuron with more pairs dominated the mean autocorrelogram.');

disp('test_locking_period_acg passed')
end
