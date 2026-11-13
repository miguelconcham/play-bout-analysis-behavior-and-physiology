function test_window_acg_methods()
%TEST_WINDOW_ACG_METHODS  Pair, FFT, and crosscorr autocorrelogram counts.
%
% The crosscorr count uses crosscorrelogram. On a run from 0 to 10 s with
% lags of +/- 1 s, reference spikes are those in [1 9] and partners are
% those in [0 10]. Play epochs use the same count on [start, end).

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
    addpath(fullfile(fileparts(this_dir), 'Custom functions'));
end

assert(strcmp(canonical_count_method(''), 'pairs'), ...
    'An empty count was not the pair-count autocorrelogram.');
assert(strcmp(canonical_count_method('crosscorrelogram'), 'crosscorr'), ...
    'The crosscorrelogram alias was not accepted.');

edges = linspace(-1, 1, 21);
centers = (edges(1:end - 1) + edges(2:end)) / 2;
spikes = [0; 4; 4.3; 10];
acg = run_window_acgs(spikes, 0, 10, edges, 'crosscorr');
near_plus = centers > 0.2 & centers < 0.4;
near_minus = centers < -0.2 & centers > -0.4;
assert(abs(sum(acg(near_plus)) - 1 / 8) < 1e-12, ...
    'The +0.3 s pair was not one count over the [1 9] reference window.');
assert(abs(sum(acg(near_minus)) - 1 / 8) < 1e-12, ...
    'The -0.3 s pair was not counted from the later spike.');
zero_bin = abs(edges(1:end - 1)) < 1e-9;
assert(any(zero_bin) && acg(find(zero_bin, 1)) == 0, ...
    'A spike paired with itself stayed in the zero-lag bin.');

outside = run_window_acgs([-2; 4; 4.3; 12], 0, 10, edges, 'crosscorr');
assert(max(abs(outside - acg)) < 1e-12, ...
    'Spikes outside [0 10] entered the crosscorr window.');

short_acg = run_window_acgs([0; 0.2], 0, 1.5, edges, 'crosscorr');
assert(all(isnan(short_acg)), ...
    'A run shorter than the two lag margins still returned a rate.');

rhythm = (0:0.3:9.9)';
rhythm_acg = run_window_acgs(rhythm, 0, 10, edges, 'crosscorr');
inner = rhythm(rhythm >= 1 & rhythm <= 9);
assert(abs(sum(rhythm_acg(near_plus)) - numel(inner) / 8) < 1e-9, ...
    'Each reference spike in [1 9] did not contribute its +0.3 s partner.');

pair_edges = -1:0.5:1;
plain = run_window_acgs([0; 1], 0, 2, pair_edges);
named = run_window_acgs([0; 1], 0, 2, pair_edges, 'pairs');
assert(max(abs(plain - named)) < 1e-12, ...
    'The default autocorrelogram was not the pair count.');

fft_acg = run_window_acgs([0; 0.5], 0, 2, pair_edges, 'fft');
assert(abs(fft_acg(4) - 1 / 1.25) < 1e-8, ...
    'The FFT autocorrelogram missed the pair half a second apart.');
assert(abs(fft_acg(2) - 1 / 1.75) < 1e-8, ...
    'The FFT autocorrelogram missed the mirrored pair.');

open_spikes = [9; 10];
closed_row = run_window_acgs(open_spikes, 0, 10, edges, 'crosscorr', true);
open_row = run_window_acgs(open_spikes, 0, 10, edges, 'crosscorr', false);
assert(max(abs(closed_row - open_row)) > 1e-12, ...
    'A spike on the closed window end was counted in a half-open epoch.');
epoch = window_condition_acg(open_spikes, [0 10], edges, [], 'crosscorr');
assert(max(abs(epoch.acg - open_row)) < 1e-12, ...
    'One play epoch did not use the half-open crosscorr count.');

rates = [1; 3];
lengths = [10; 6];
pooled = pool_crosscorr_rows(rates, lengths, edges);
assert(abs(pooled(1) - (1 * 8 + 3 * 4) / 12) < 1e-12, ...
    'Crosscorr rows were not pooled by inner duration.');
short_pool = pool_crosscorr_rows([1; 9], [10; 1.5], edges);
assert(abs(short_pool(1) - 1) < 1e-12, ...
    'A window shorter than two lag margins entered the crosscorr pool.');

[rebuilt, ~, corrected] = firing_run_acgs([0; 1], 0, 2, pair_edges, 'overlap', [], [], 'fft');
direct = run_window_acgs([0; 1], 0, 2, pair_edges, 'fft', true);
assert(corrected && max(abs(rebuilt - direct)) < 1e-12, ...
    'The firing-run estimator did not rebuild the FFT count.');
assert_error(@() firing_run_acgs([0; 1], 0, 2, pair_edges, 'eligible', [], [], 'fft'), ...
    'An eligible estimate accepted an overlap count.');

assert(strcmp(locking_summary_filename(false, 'overlap', 'fft'), ...
    'locking_period_acg_summary_fft.mat'), ...
    'The FFT summary did not get its own file name.');
assert(strcmp(locking_summary_filename(false, 'overlap', 'pairs'), ...
    'locking_period_acg_summary.mat'), ...
    'Pair counts changed the overlap summary name.');
assert(strcmp(locking_tail_filename(25, 'overlap', 'crosscorr'), ...
    'locking_period_acg_tails_crosscorr_25.mat'), ...
    'The crosscorr tail file did not get its own file name.');
assert(strcmp(locking_shuffle_filename(500, 'overlap', 'fft'), ...
    'locking_period_acg_shuffle_fft_n500.mat'), ...
    'The FFT shuffle checkpoint did not get its own file name.');
assert(strcmp(play_power_acg_filename('A', 'overlap', 'crosscorr'), ...
    'A_play_power_acg_crosscorr.mat'), ...
    'The crosscorr play file did not get its own file name.');
assert(strcmp(play_power_acg_filename('*', 'overlap', 'fft'), ...
    '*_play_power_acg_fft.mat'), ...
    'The FFT play pattern would also match pair-count files.');
assert(strcmp(play_power_retry_filename('overlap', 'fft'), ...
    'play_power_acg_fft_retry_sessions.mat'), ...
    'The FFT retry log did not get its own file name.');
assert_error(@() acg_count_token('eligible', 'fft'), ...
    'An eligible file name accepted an overlap count.');

locking_text = fileread(fullfile(this_dir, 'ACG power and locking type', ...
    'Analyze_locking_period_acg.m'));
count_token = regexp(locking_text, '(?m)^acg_count\s*=\s*''([^'']+)''', 'tokens', 'once');
assert(strcmp(count_token{1}, 'pairs'), ...
    'The locking script does not default to pair counts.');
estimate_text = fileread(fullfile(this_dir, 'Estimate_play_power_acg.m'));
estimate_count = regexp(estimate_text, '(?m)^acg_count\s*=\s*''([^'']+)''', 'tokens', 'once');
assert(strcmp(estimate_count{1}, 'pairs'), ...
    'The play estimate does not default to pair counts.');

threw = false;
try
    run_window_acgs([0; 1], 0, 2, edges, 'wavelet');
catch ME
    threw = strcmp(ME.identifier, 'canonical_count_method:name');
end
assert(threw, 'An unknown autocorrelogram count was accepted.');

disp('test_window_acg_methods passed')
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
