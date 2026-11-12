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

disp('test_shuffle_percentile_plot passed')
end


function rec = make_record(group_name, lag, curve)
rec.lock_group = group_name;
rec.locking.low.pctl_lag = lag;
rec.locking.low.sum_pctl = curve;
end
