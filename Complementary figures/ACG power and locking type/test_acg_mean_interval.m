function test_acg_mean_interval()
%TEST_ACG_MEAN_INTERVAL  Mean and 95% interval across autocorrelogram rows.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

matrix = [1 2; 3 4; 5 6];
interval = acg_mean_interval(matrix);
assert(abs(interval.mean(1) - 3) < 1e-12, 'The first lag mean is not 3.');
assert(abs(interval.mean(2) - 4) < 1e-12, 'The second lag mean is not 4.');
expected_half = mean_t_critical_for_test(0.025, 2) * std([1; 3; 5], 0) / sqrt(3);
assert(abs(interval.half(1) - expected_half) < 1e-12, ...
    'The confidence interval is not the 95% t interval of the mean.');
assert(abs(interval.low(1) - (3 - expected_half)) < 1e-12, ...
    'The lower edge is not the mean minus the half-width.');
assert(abs(interval.high(2) - (4 + interval.half(2))) < 1e-12, ...
    'The upper edge is not the mean plus the half-width.');

with_gap = [1 nan; 3 4; 5 8];
gapped = acg_mean_interval(with_gap);
assert(gapped.n(1) == 3 && gapped.n(2) == 2, ...
    'A missing value was still counted at that lag.');
assert(abs(gapped.mean(2) - 6) < 1e-12, ...
    'The mean did not skip the missing value.');

one_row = acg_mean_interval([2 5]);
assert(abs(one_row.mean(1) - 2) < 1e-12, 'One neuron did not keep its value as the mean.');
assert(~isfinite(one_row.half(1)), 'One neuron still received a confidence interval.');
assert(abs(one_row.low(1) - one_row.mean(1)) < 1e-12, ...
    'A missing interval did not stay on the mean.');

empty_interval = acg_mean_interval([]);
assert(isempty(empty_interval.mean), 'An empty autocorrelogram still returned a mean.');

script_text = fileread(fullfile(this_dir, ...
    'Plot_acg_population_nonoscillatory_description.m'));
assert(~isempty(strfind(script_text, ...
    'draw_lock_type_acg(imagesc_logit.sum.all, lag, neuron_type, osc_mask,')), ...
    'Section 4 does not draw the oscillatory lock-type autocorrelogram.');
assert(~isempty(strfind(script_text, ...
    'draw_lock_type_acg(imagesc_logit.sum.all, lag, neuron_type, ~osc_mask,')), ...
    'Section 4 does not draw the non-oscillatory lock-type autocorrelogram.');
assert(~isempty(strfind(script_text, 'imagesc_logit.sum.all')), ...
    'The lock-type figure does not use the all-runs autocorrelogram.');

disp('test_acg_mean_interval passed')
end


function tc = mean_t_critical_for_test(tail_p, df)
x = betaincinv(2 * tail_p, df / 2, 0.5);
tc = sqrt(df * (1 - x) / x);
end
