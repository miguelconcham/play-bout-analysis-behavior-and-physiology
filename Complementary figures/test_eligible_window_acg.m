function test_eligible_window_acg()
%TEST_ELIGIBLE_WINDOW_ACG  Eligible-spike counts without experiment files.
%
% One spike that can see the bin and has one partner, one spike that can
% see it and has none, and two spikes that cannot see it. The mean partner
% count is 1/2. Dividing by the same mean rate times the bin width is 1.25.
% A window whose spikes all sit too late is missing, not zero.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
    addpath(fullfile(fileparts(this_dir), 'Figure 2', 'Supporting codes', 'Delta episodes'));
    addpath(fullfile(fileparts(this_dir), 'Custom functions'));
end

assert(strcmp(canonical_acg_method(''), 'eligible'), ...
    'An empty method was not the eligible-spike estimator.');
assert(acg_rate_exponent('eligible') == 1, ...
    'The eligible-spike Poisson ratio did not use rate to the first power.');
assert(acg_rate_exponent('overlap') == 2 && acg_rate_exponent('') == 2, ...
    'The overlap ratio did not keep rate squared when the method is absent.');
assert_error(@() canonical_acg_method('counts'), ...
    'An unknown autocorrelogram method was accepted.');

bin_edges = [0.45 0.55];
spikes = [0.10; 0.20; 0.60; 0.80];
[partner_counts, n_eligible] = eligible_window_counts(spikes, 0, 1, bin_edges, false);
assert(partner_counts == 1 && n_eligible == 2, ...
    'The two early spikes did not give one partner count and two eligible spikes.');
assert(abs((partner_counts / n_eligible) / (4 * 0.10) - 1.25) < 1e-12, ...
    'The mean partner count was not 1.25 after dividing by rate times bin width.');

second = [1.05; 1.55];
[both_counts, both_n] = eligible_window_counts( ...
    [spikes; second], [0; 1], [1; 2], bin_edges, false);
pooled = pool_eligible_period_acgs(both_counts, both_n, [1; 1], [1; 1]);
assert(abs(pooled.all - (2 / 3)) < 1e-12, ...
    'Two epochs did not pool to two partner counts over three eligible spikes.');
assert(abs((pooled.all) / (3 * 0.10) - (2 / 3) / 0.3) < 1e-12, ...
    'The pooled ratio did not use the mean rate of both epochs.');

[late_counts, late_n] = eligible_window_counts([0.70; 0.90], 0, 1, bin_edges, false);
assert(late_counts == 0 && late_n == 0, ...
    'Spikes too late to see the bin were not left out.');
late_pool = pool_eligible_period_acgs(late_counts, late_n, 1, 1);
assert(isnan(late_pool.all), ...
    'A bin with no eligible spike was stored as zero.');

[closed_counts, closed_n] = eligible_window_counts(1, 0, 1, [-0.2 -0.1], true);
[open_counts, open_n] = eligible_window_counts(1, 0, 1, [-0.2 -0.1], false);
assert(closed_n == 1 && closed_counts == 0, ...
    'A spike on the closed window end was not eligible for a bin it can see.');
assert(open_n == 0 && open_counts == 0, ...
    'A spike on the half-open window end was counted.');

scores = (1:8)';
eligible_pool = pool_eligible_period_acgs(zeros(8, 1), ones(8, 1), ones(8, 1), scores);
assert(isequal(eligible_pool.low_index, [true; true; false; false; false; false; false; false]), ...
    'The low 25% was not the two weakest windows.');
assert(isequal(eligible_pool.high_index, [false; false; false; false; false; false; true; true]), ...
    'The high 25% was not the two strongest windows.');
assert(all(eligible_pool.all == 0), ...
    'Windows with eligible spikes and no partners were not a partner count of zero.');

if exist('normcdf', 'file') ~= 2
    pkg load statistics
end
centers = [0.40 0.50];
curve = [0.5 0.5];
options = struct('refractory_s', 0, 'smooth_half', 0.06, 'skip_abs_lag', 0.02, ...
    'delta_lag', [0.2 0.5], 'cluster_alpha', 0.05, 'cluster_perm', 20, ...
    'acg_method', 'eligible');
compared = compare_play_power_acgs(curve, curve, centers, 4, 4, 'ratio', options);
assert(max(abs(compared.curve_a - 1.25)) < 1e-9, ...
    'The eligible rate ratio was not the mean partner count over rate times bin width.');
options = rmfield(options, 'acg_method');
overlap_compared = compare_play_power_acgs(curve, curve, centers, 4, 4, 'ratio', options);
assert(max(abs(overlap_compared.curve_a - 0.3125)) < 1e-9, ...
    'A ratio without a method did not keep the overlap rate-squared scale.');

block = condition_window_acg(spikes, [0 1], bin_edges, [], 'eligible');
assert(abs(block.acg - 0.5) < 1e-12 && abs(block.rate - 4) < 1e-12, ...
    'The play-condition curve was not the mean partner count, or the rate changed.');
wide = -1:0.1:1;
overlap_block = condition_window_acg([0.05; 0.20], [0 1], wide, 4, 'overlap');
direct = window_condition_acg([0.05; 0.20], [0 1], wide, 4);
assert(max(abs(overlap_block.acg - direct.acg)) < 1e-12, ...
    'The overlap selector changed the play-condition autocorrelogram.');

assert(strcmp(locking_summary_filename(false, 'overlap'), 'locking_period_acg_summary.mat'), ...
    'The overlap summary name changed.');
assert(strcmp(locking_summary_filename(true, 'overlap'), 'locking_period_acg_summary_max_strata.mat'), ...
    'The maximized overlap summary name changed.');
assert(strcmp(locking_summary_filename(false, 'eligible'), 'locking_period_acg_summary_eligible.mat'), ...
    'The eligible summary did not get its own file name.');
assert(strcmp(locking_summary_filename(true, 'eligible'), ...
    'locking_period_acg_summary_eligible_max_strata.mat'), ...
    'The maximized eligible summary did not get its own file name.');
assert(strcmp(locking_tail_filename(25, 'overlap'), 'locking_period_acg_tails_25.mat'), ...
    'The overlap tail file name changed.');
assert(strcmp(locking_tail_filename(10, 'eligible'), 'locking_period_acg_tails_eligible_10.mat'), ...
    'The eligible tail file did not get its own file name.');
assert(strcmp(locking_shuffle_filename(500, 'overlap'), 'locking_period_acg_shuffle_n500.mat'), ...
    'The overlap shuffle checkpoint name changed.');
assert(strcmp(locking_shuffle_filename(500, 'eligible'), ...
    'locking_period_acg_shuffle_eligible_n500.mat'), ...
    'The eligible shuffle checkpoint did not get its own file name.');
assert(strcmp(play_power_acg_filename('A', 'overlap'), 'A_play_power_acg.mat'), ...
    'The overlap play file name changed.');
assert(strcmp(play_power_acg_filename('A', 'eligible'), 'A_play_power_acg_eligible.mat'), ...
    'The eligible play file did not get its own file name.');
assert(strcmp(play_power_acg_filename('*', 'eligible'), '*_play_power_acg_eligible.mat'), ...
    'The eligible play pattern would also match overlap files.');
assert(strcmp(play_power_retry_filename('overlap'), 'play_power_acg_retry_sessions.mat'), ...
    'The overlap retry log name changed.');
assert(strcmp(play_power_retry_filename('eligible'), 'play_power_acg_eligible_retry_sessions.mat'), ...
    'The eligible retry log did not get its own file name.');
locking_text = fileread(fullfile(this_dir, 'Analyze_locking_period_acg.m'));
locking_method = regexp(locking_text, '(?m)^acg_method\s*=\s*''([^'']+)''', 'tokens', 'once');
assert(strcmp(locking_method{1}, 'eligible'), ...
    'The locking script does not default to the eligible-spike method.');
tail_token = regexp(locking_text, '(?m)^tail_percent\s*=\s*([^;%]+)', 'tokens', 'once');
assert(~isempty(tail_token) && str2num(tail_token{1}) == 25, ...
    'tail_percent was not a number the shuffle can read.');
low_token = regexp(locking_text, '(?m)^low_percent\s*=\s*([^;%]+)', 'tokens', 'once');
assert(~isempty(low_token) && isempty(str2num(low_token{1})), ...
    'low_percent is still a literal, so the shuffle would not follow tail_percent.');
estimate_text = fileread(fullfile(this_dir, 'Estimate_play_power_acg.m'));
estimate_method = regexp(estimate_text, '(?m)^acg_method\s*=\s*''([^'']+)''', 'tokens', 'once');
assert(strcmp(estimate_method{1}, 'eligible'), ...
    'The play estimate does not default to the eligible-spike method.');
analyze_text = fileread(fullfile(this_dir, 'Analyze_play_power_acg.m'));
analyze_method = regexp(analyze_text, '(?m)^acg_method\s*=\s*''([^'']+)''', 'tokens', 'once');
assert(strcmp(analyze_method{1}, 'eligible'), ...
    'The play analysis does not default to the eligible-spike method.');

disp('test_eligible_window_acg passed')
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
