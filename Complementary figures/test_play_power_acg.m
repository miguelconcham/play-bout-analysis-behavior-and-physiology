function test_play_power_acg()
%TEST_PLAY_POWER_ACG  Checks that do not need the experiment files.
%
% Epochs inside and outside play, the 25/75 power tails, one-to-one power
% matching, overlap-weighted autocorrelograms, and the paired property and
% cluster tests. Run this file by name after add_repo_paths.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
    addpath(fullfile(fileparts(this_dir), 'Figure 2', 'Supporting codes', 'Delta episodes'));
    addpath(fullfile(fileparts(this_dir), 'Custom functions'));
end

assert(isequal(merge_time_intervals([0 1; 1 2]), [0 2]), ...
    'Touching intervals were not joined.');
assert(isequal(merge_time_intervals([0 1.5; 1 2; 4 5]), [0 2; 4 5]), ...
    'Overlapping intervals were not joined.');
assert(isequal(intersect_time_intervals([0 10], [2 5; 8 12]), [2 5; 8 10]), ...
    'The time-limit intersection was wrong.');
assert(isequal(outside_play_intervals([0 10], [2 4; 7 8]), [0 2; 4 7; 8 10]), ...
    'Gaps outside play bouts were wrong.');
assert(isequal(outside_play_intervals([0 10], [-1 1]), [1 10]), ...
    'A play bout that starts before the window was not clipped.');
assert(isequal(fixed_duration_epochs([0 2.5], 1), [0 1; 1 2]), ...
    'A short leftover became an epoch.');
epochs = play_and_outside_epochs([0 10], [2 4; 7 8], 1, []);
assert(isequal(epochs.play, [2 3; 3 4; 7 8]), ...
    'Play epochs crossed a bout or kept a short leftover.');
assert(isequal(epochs.outside, [0 1; 1 2; 4 5; 5 6; 6 7; 8 9; 9 10]), ...
    'Outside epochs were not the one-second blocks between play bouts.');
clipped = play_and_outside_epochs([0 10], [0 10], 1, [2 5.5]);
assert(isequal(clipped.play, [2 3; 3 4; 4 5]), ...
    'Epochs were not clipped to the LFP time limit.');
assert(isempty(clipped.outside), ...
    'Play that fills the clipped window still left outside epochs.');

power = [1 2 3 4 40 50 60 70]';
tails = play_power_tails(power, 25, 75);
assert(isequal(tails.low, [true; true; false; false; false; false; false; false]), ...
    'The low 25% of play power was not the two weakest epochs.');
assert(isequal(tails.high, [false; false; false; false; false; false; true; true]), ...
    'The high 25% of play power was not the two strongest epochs.');
tied = play_power_tails(ones(6, 1), 25, 75);
assert(~any(tied.low) && ~any(tied.high), ...
    'Equal powers still produced a low/high split.');
short = play_power_tails([1 2 3]', 25, 75);
assert(~any(short.low) && ~any(short.high), ...
    'Three values still produced a low/high split.');

pairs = match_power_periods([0; 0.3], [0.2; 0.5], 0.25);
assert(isequal(pairs, [1 1; 2 2]), ...
    'The caliper match kept the closest pair and dropped the second valid pair.');
pairs = match_power_periods([0; 10], [0.1; 10.2; 100], 1);
assert(isequal(pairs, [1 1; 2 2]), ...
    'A candidate outside the caliper was paired.');
pairs = match_power_periods([0; 0], [0.1], 0.2);
assert(size(pairs, 1) == 1, 'One outside period was paired twice.');
pairs = match_power_periods([0; 1], 10, 0.5);
assert(isempty(pairs), 'A period beyond the caliper was paired.');
assert_error(@() match_power_periods(0, 0, -1), ...
    'A negative caliper was accepted.');

play_db = [0; 0; 10; 10];
outside_db = [0.5; 9.5; 3];
assigned = assign_play_power_conditions(10 .^ (play_db / 10), 10 .^ (outside_db / 10), ...
    struct('caliper_sd', 0.25));
assert(isequal(assigned.play_low, [1; 2]) && isequal(assigned.play_high, [3; 4]), ...
    'The outer play tails were not the two low and two high epochs.');
assert(numel(assigned.play_matched) == 2, ...
    'The overall match did not keep the two outside epochs inside the caliper.');
assert(numel(assigned.play_low_matched) == 1 && numel(assigned.play_high_matched) == 1, ...
    'Each power tail did not keep its one matched outside epoch.');
assert(~ismember(3, [assigned.outside_matched; assigned.outside_low_matched; assigned.outside_high_matched]), ...
    'The outside epoch beyond the caliper was used.');
assert(abs(assigned.gap_db.play_matched - 0.5) < 1e-9, ...
    'The matched decibel gap was not 0.5.');
assert(max(abs(play_db(assigned.play_matched) - outside_db(assigned.outside_matched))) <= assigned.caliper_db + 1e-9, ...
    'A matched pair exceeded the caliper.');

edges = [0.05 0.15];
pooled_counts = pool_selected_window_acgs([90; 10], edges, [1; 0.3], false);
assert(abs(pooled_counts - (100 / 1.1)) < 1e-9, ...
    'The pooled autocorrelogram was not total counts over total overlap.');
pooled_rates = pool_selected_window_acgs([100; 50], edges, [1; 0.3], true);
assert(abs(pooled_rates - pooled_counts) < 1e-9, ...
    'Edge-corrected rows were not converted back to counts before pooling.');

bin_edges = -1:0.1:1;
one = window_condition_acg([0.05; 0.20], [0 1], bin_edges, 4);
two = window_condition_acg([0.05; 0.20; 2.05; 2.20], [0 1; 2 3], bin_edges, [4; 8]);
centers = bin_edges(1:end - 1) + 0.05;
hit = find(one.acg > 0);
assert(~isempty(hit), 'A two-spike epoch produced no autocorrelogram pairs.');
assert(all(abs(two.acg(hit) - one.acg(hit)) < 1e-9), ...
    'Two equal epochs did not pool to the same overlap-weighted rate.');
positive = hit(centers(hit) > 0);
assert(numel(positive) == 1, 'The spike pair did not land in one positive lag bin.');
assert(abs(one.acg(positive) - 1 / (1 - abs(centers(positive)))) < 1e-9, ...
    'The pair count was not divided by the overlap at that lag.');
assert(one.n_spikes == 2 && abs(one.duration - 1) < 1e-12 && abs(one.rate - 2) < 1e-12, ...
    'Spike count, duration, or rate of one epoch was wrong.');
assert(abs(one.mean_power - 4) < 1e-12, 'A single epoch did not keep its own power.');
assert(two.n_spikes == 4 && abs(two.duration - 2) < 1e-12 && abs(two.rate - 2) < 1e-12, ...
    'Spike count, duration, or rate of the two epochs was wrong.');
assert(abs(two.mean_power - 6) < 1e-12, 'Epoch power was not the mean of the two epochs.');
assert(abs(two.mean_power_db - mean(10 * log10([4 8]))) < 1e-12, ...
    'Epoch decibel power was not the mean of the epoch decibel values.');
boundary = window_condition_acg(1, [0 1; 1 2], bin_edges, []);
assert(boundary.n_spikes == 1, 'A spike on the shared epoch edge was counted twice or not at all.');
squared = epoch_mean_square(2 * ones(1, 5), 10, [0.1 0.3; 0.15 0.35]);
assert(abs(squared(1) - 4) < 1e-12 && abs(squared(2) - 4) < 1e-12, ...
    'Epoch mean square did not use the half-open LFP sample clock.');

centers = -0.02:0.01:0.02;
blanked = blank_acg_center(ones(1, numel(centers)), centers, 0.005);
assert(isnan(blanked(centers == 0)) && all(blanked(abs(centers) >= 0.01) == 1), ...
    'The refractory blank did not cover only the lag-0 bin.');
assert(strcmp(canonical_kept_area('SupCol'), 'SC'), 'SupCol was not mapped to SC.');
assert(strcmp(canonical_kept_area('superior colliculus'), 'SC'), ...
    'The superior-colliculus label was not mapped to SC.');
assert(strcmp(canonical_kept_area({'SC'}), 'SC'), 'A cell label was not unwrapped.');
assert(strcmp(canonical_kept_area('lpag'), 'LPAG'), 'LPAG was not recognized.');
assert(strcmp(canonical_kept_area('VLPAG'), 'VLPAG'), 'VLPAG was rewritten as LPAG.');
assert(isempty(canonical_kept_area('PAG')) && isempty(canonical_kept_area('isRt')), ...
    'PAG or isRt was kept.');

scores_a = fake_scores([1; 2], [0.2; 0.3], [true; true], [false; true]);
scores_b = fake_scores([1; 2], [0.4; nan], [true; false], [true; true]);
assert(isequal(acg_measure_mask(scores_a, scores_b, 'first_prominence'), [true; false]), ...
    'A first peak that was not searched on both sides was paired.');
assert(isequal(acg_measure_mask(scores_a, scores_b, 'first_lag'), [false; false]), ...
    'A missing first peak still contributed a latency.');
assert(~any(acg_measure_mask(scores_a, scores_b, 'second_prominence')), ...
    'A second cycle that was not observed was paired.');
p_same = paired_signrank_p([1; 2; 3], [1; 2; 3], 'both');
assert(isnan(p_same), 'Identical pairs produced a p-value.');

if exist('normcdf', 'file') ~= 2
    pkg load statistics
end
a = (1:12)';
p_increase = paired_signrank_p(a + 1, a, 'right');
assert(p_increase < 0.01, 'A consistent increase was not detected.');

centers = -1:0.01:1;
n = 10;
osc = zeros(1, numel(centers));
far = abs(centers) >= 0.1;
osc(far) = cos(2 * pi * 3 * centers(far));
side_a = ones(n, numel(centers));
side_b = side_a + 0.8 * osc + 0.01 * (1:n)';
options = struct('cluster_perm', 99);
compared = compare_play_power_acgs(side_a, side_b, centers, ones(n, 1), ones(n, 1), 'sum', options);
assert(compared.n == n, 'A finite autocorrelogram was dropped.');
assert(all(compared.scores_b.acg_delta > compared.scores_a.acg_delta), ...
    'A 3 Hz side peak did not raise folded autocorrelogram delta power.');
assert(compared.properties.acg_delta.median_delta > 0, ...
    'The paired delta-power difference was not positive.');
assert(compared.properties.acg_delta.p_increase < 0.05, ...
    'The paired delta-power increase was not detected.');
assert(any([compared.clusters.clusters.survives]), ...
    'The time-wise cluster test did not keep the 3 Hz difference.');

disp('test_play_power_acg passed')
end


function scores = fake_scores(prominence, lag, measurable, present)
n = numel(prominence);
scores.burst = ones(n, 1);
scores.first_measurable = measurable(:);
scores.first_present = present(:);
scores.first_lag = lag(:);
scores.first_prominence = prominence(:);
scores.second_measurable = false(n, 1);
scores.second_present = false(n, 1);
scores.second_lag = nan(n, 1);
scores.second_prominence = nan(n, 1);
scores.n_cycles = zeros(n, 1);
scores.acg_delta = ones(n, 1);
end


function assert_error(fn, message)
threw = false;
try
    fn();
catch
    threw = true;
end
assert(threw, message);
end
