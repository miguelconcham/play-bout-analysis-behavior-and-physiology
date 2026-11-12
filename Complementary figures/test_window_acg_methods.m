function test_window_acg_methods()
%TEST_WINDOW_ACG_METHODS  The three autocorrelogram choices.
%
% The crosscorr choice uses crosscorrelogram. On a run from 0 to 10 s
% with lags of +/- 1 s, reference spikes are those in [1 9] and partners
% are those in [0 10].

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
    addpath(fullfile(fileparts(this_dir), 'Custom functions'));
end

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

fft_edges = -1:0.5:1;
fft_acg = run_window_acgs([0; 0.5], 0, 2, fft_edges, 'fft');
assert(abs(fft_acg(4) - 1 / 1.25) < 1e-8, ...
    'The FFT autocorrelogram missed the pair half a second apart.');
assert(abs(fft_acg(2) - 1 / 1.75) < 1e-8, ...
    'The FFT autocorrelogram missed the mirrored pair.');

threw = false;
try
    run_window_acgs([0; 1], 0, 2, edges, 'wavelet');
catch ME
    threw = strcmp(ME.identifier, 'run_window_acgs:method');
end
assert(threw, 'An unknown autocorrelogram method was accepted.');

disp('test_window_acg_methods passed')
end
