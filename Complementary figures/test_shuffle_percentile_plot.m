function test_shuffle_percentile_plot()
%TEST_SHUFFLE_PERCENTILE_PLOT  Percentile rows, tail rates, and finite logits.
%
% These checks do not need the experiment files. Run this file by name
% after add_repo_paths.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

assert(strcmp(shuffle_checkpoint_filename(500), 'locking_period_acg_shuffle_n500.mat'), ...
    'A missing method did not keep the overlap checkpoint name.');
assert(strcmp(shuffle_checkpoint_filename(500, 'overlap'), 'locking_period_acg_shuffle_n500.mat'), ...
    'Overlap did not keep locking_period_acg_shuffle_nN.mat.');
assert(strcmp(shuffle_checkpoint_filename(500, 'eligible'), ...
    'locking_period_acg_shuffle_eligible_n500.mat'), ...
    'Eligible did not use its own checkpoint name.');

n_shuffles = 500;
logit = correct_percentile_logit([0 50 100 NaN], n_shuffles);
assert(all(isfinite(logit(1:3))), 'A percentile of 0 or 100 made the logit infinite.');
assert(isnan(logit(4)), 'A missing percentile became a number.');
assert(abs(logit(2)) < 1e-12, 'The logit of the 50th percentile was not 0.');
edge = 50 / (n_shuffles + 1);
assert(abs(logit(1) - log10((edge / 100) / (1 - edge / 100))) < 1e-12, ...
    'A percentile of 0 was not clipped to half a count.');
assert(abs(logit(3) + logit(1)) < 1e-12, ...
    'Percentiles of 0 and 100 were not clipped to opposite logits.');
assert(logit(1) < 0 && logit(3) > 0, 'The clipped logits had the wrong sign.');

stack = [0; 90];
clipped = correct_percentile_logit(stack, 9);
assert(abs(mean(clipped) - correct_percentile_logit(mean(stack), 9)) > 0.05, ...
    'The mean logit collapsed to the logit of the mean percentile.');

fractions = [96 10; 94 4; 4 96; NaN 5];
[high_percent, low_percent] = percentile_tail_fractions(fractions);
assert(abs(high_percent(1) - 100 / 3) < 1e-12, ...
    'The share at or above 95 ignored a finite neuron.');
assert(abs(low_percent(1) - 100 / 3) < 1e-12, ...
    'The share at or below 5 ignored a finite neuron.');
assert(abs(high_percent(2) - 25) < 1e-12, ...
    'One neuron at 96 was not 25 percent of four.');
assert(abs(low_percent(2) - 50) < 1e-12, ...
    'Percentile 5 was left out of the low tail.');

lag = [0.1 0.2 0.3];
early = make_record('trough', lag, [10 20 5]);
late = make_record('trough', lag, [5 15 40]);
weaker = make_record('trough', lag, [1 2 3]);
[got_lag, got_stack] = percentile_curve_stack({late, weaker, early}, 'locking', 'low', 'sum_pctl');
assert(isequal(got_lag, lag), 'The percentile lag axis was not kept.');
assert(isequal(size(got_stack), [3 3]), 'A finite percentile curve was dropped.');
assert(abs(got_stack(1, 2) - 20) < 1e-12, ...
    'The earliest peak was not the first row.');
assert(abs(got_stack(2, 3) - 40) < 1e-12 && abs(got_stack(3, 3) - 3) < 1e-12, ...
    'A shared peak lag did not put the higher percentile first.');
mismatch = make_record('trough', [0.1 0.2], [8 9]);
[~, matched] = percentile_curve_stack({early, mismatch}, 'locking', 'low', 'sum_pctl');
assert(isequal(size(matched), [1 3]), 'A curve on a different lag axis was kept.');

lag = [0.05 0.25 0.35 0.60];
unsorted = [90 10 10 90; 10 80 70 10; NaN NaN NaN 40];
[sorted_pctl, sorted_logit, score, order] = order_percentile_rows(unsorted, lag, [0.2 0.4], 500);
assert(abs(score(1) - mean([80 70])) < 1e-12, ...
    'The sort score was not the mean percentile from 200 to 400 ms.');
assert(abs(score(2) - mean([10 10])) < 1e-12, ...
    'Lags outside 200 to 400 ms changed the sort score.');
assert(isequal(sorted_pctl(1, :), unsorted(2, :)), ...
    'The higher 200-400 ms percentile was not the first row.');
assert(order(1) == 2, ...
    'The sort index did not point at the higher 200-400 ms row.');
assert(all(isnan(sorted_pctl(3, 1:3))), ...
    'A neuron with no percentile in the window was not placed last.');
assert(all(isfinite(sorted_logit(1:2, :))), ...
    'A sorted percentile of 0 or 100 made the logit infinite.');

lag = [0.10 0.25 0.35];
kept_neuron = population_record('trough', lag, [1 2 3]);
missing_low = population_record('peak', lag, [4 5 6]);
missing_low.power = rmfield(missing_low.power, 'low');
unlocked = population_record('nonentrained', lag, [7 8 9]);
unclassified = population_record('', lag, [0 1 1]);
other_axis = population_record('peak', [0.10 0.20], [1 2]);
[got_lag, got_percentile, got_type, got_index] = align_population_percentiles( ...
    {kept_neuron, missing_low, unlocked, unclassified, other_axis});
assert(isequal(got_lag, lag), 'The shared lag axis was not kept.');
assert(isequal(got_index, [1; 3; 4]), ...
    'A neuron missing one curve, or on another lag axis, was kept.');
assert(isequal(got_type, {'trough'; 'nonentrained'; ''}), ...
    'Lock group was not kept in the same row order as the curves.');
assert(isequal(got_percentile.sum.all(2, :), [7 8 9]), ...
    'The all-runs row was not the same neuron as the lock-group row.');
assert(isequal(got_percentile.sum.power_high(2, :), [9 10 11]), ...
    'The high-power row was not aligned with the all-runs row.');
assert(isequal(got_percentile.ratio.power_low(1, :), [6 7 8]), ...
    'The rate-ratio low-power row was not aligned with the sum-to-1 row.');

osc_lag = [0.10 0.30 0.40 0.50];
osc_logit = [0 1.5 0 0; 0 0 0 0; 0 0.4 0 0; NaN NaN 0 0];
[osc_rows, osc_score] = oscillatory_row_index( ...
    osc_logit, osc_lag, [0.25 0.35], [0.05 0.15], 0);
assert(isequal(osc_rows, [1; 3]), ...
    'Oscillatory rows were not the positive peak-minus-trough scores, strongest first.');
assert(abs(osc_score(1) - 1.5) < 1e-12 && abs(osc_score(3) - 0.4) < 1e-12, ...
    'The oscillation score was not the peak-window mean minus the trough-window mean.');
assert(abs(osc_score(2)) < 1e-12, ...
    'A flat logit did not score zero.');

[type_count, type_proportion, type_label] = lock_type_proportion( ...
    {'Trough', 'peak', 'nonentrained', '', 'delta'});
assert(isequal(type_count, [1 1 1 2]), ...
    'Trough, peak, unlocked, and non-classified were not counted separately.');
assert(abs(type_proportion(4) - 0.4) < 1e-12, ...
    'An unlabeled neuron and an unknown label were not both non-classified.');
assert(isequal(type_label, {'Trough', 'Peak', 'Unlocked', 'Non-classified'}), ...
    'The lock-type labels were renamed.');

population_type = {'trough', 'trough', 'peak', 'nonentrained', ''};
[type_percent, type_observed, type_total] = lock_type_osc_rate(population_type, [1; 3]);
assert(isequal(type_observed, [1 1 0 0]), ...
    'The oscillatory count was not the number of that lock type inside the group.');
assert(isequal(type_total, [2 1 1 1]), ...
    'The lock-type total was not the number of that type in the population.');
assert(abs(type_percent(1) - 50) < 1e-12 && abs(type_percent(2) - 100) < 1e-12, ...
    'The bar height was not 100 times oscillatory cells over cells of that type.');
assert(abs(type_percent(3)) < 1e-12 && abs(type_percent(4)) < 1e-12, ...
    'A lock type with no oscillatory cell was not zero percent.');

member = [true; false; true; true];
[member_rows] = oscillatory_member_rows(member, [1; 9; 3; NaN]);
assert(isequal(member_rows, [3; 1; 4]), ...
    'Oscillatory rows were not strongest modulation first, with missing scores last.');

positive_lag = (0.01:0.01:1)';
rhythm = 2 * ones(size(positive_lag));
rhythm(10) = 0;
rhythm(30) = 5;
rhythm(40) = 1;
rhythm(60) = 4;
R = quantify_acg_rhythmicity(rhythm, positive_lag);
assert(abs(R.FirstModulation - 5) < 1e-12, ...
    'The first modulation was not the first peak minus the first trough.');
assert(abs(R.SecondModulation - 3) < 1e-12, ...
    'The second modulation was not the second peak minus the second trough.');
assert(abs(R.Persistence - 0.6) < 1e-12, ...
    'Persistence was not the second modulation divided by the first.');
assert(abs(R.PeakInterval - (positive_lag(60) - positive_lag(30))) < 1e-12 ...
    && abs(R.TroughInterval - (positive_lag(40) - positive_lag(10))) < 1e-12, ...
    'The peak and trough intervals did not use the lag times.');
assert(abs(R.Baseline - 2) < 1e-12, ...
    'Baseline was not the median of the 0.8 to 1 s tail.');
assert(isfinite(R.OscillationScore) && R.OscillationScore > 1, ...
    'The spectral oscillation score was not above the mean spectrum.');

full_lag = [-flipud(positive_lag); positive_lag];
full_acg = [zeros(numel(positive_lag), 1); rhythm];
if exist('array2table', 'file') == 2 && exist('movevars', 'file') == 2
    feature_table = acg2features(full_acg', full_lag);
    assert(feature_table.NeuronIndex == 1, ...
        'The feature table did not keep the neuron row.');
    assert(abs(feature_table.Modulation1 - R.FirstModulation) < 1e-12, ...
        'Negative lags changed the positive-lag modulation.');
    assert(abs(feature_table.OscillationScore - R.OscillationScore) < 1e-9, ...
        'Negative lags changed the oscillation score.');
end

disp('test_shuffle_percentile_plot passed')
end


function rec = make_record(group_name, lag, curve)
rec.lock_group = group_name;
rec.locking.low.pctl_lag = lag;
rec.locking.low.sum_pctl = curve;
end


function rec = population_record(group_name, lag, curve)
rec.lock_group = group_name;
rec.all.pctl_lag = lag;
rec.all.sum_pctl = curve;
rec.all.ratio_pctl = curve + 1;
rec.power.high.pctl_lag = lag;
rec.power.high.sum_pctl = curve + 2;
rec.power.high.ratio_pctl = curve + 3;
rec.power.low.pctl_lag = lag;
rec.power.low.sum_pctl = curve + 4;
rec.power.low.ratio_pctl = curve + 5;
end
