function test_analyze_session_acg()
%TEST_ANALYZE_SESSION_ACG  Percentile rows and the shuffle band.
%
% These checks do not need the experiment files or a figure window. Run
% this file by name after add_repo_paths, or run it directly.

this_dir = fileparts(mfilename('fullpath'));
complementary_dir = fileparts(this_dir);
repo_root = fileparts(complementary_dir);
addpath(this_dir);
addpath(complementary_dir);
addpath(fullfile(complementary_dir, 'ACG power and locking type'));

column = (1:100)';
band = session_acg_null_band([column, nan(100, 1)]);
assert(abs(band.median(1) - median(column)) < 1e-12, ...
    'The shuffle median was not the middle of the finite rows.');
assert(abs(band.low(1) - prctile(column, 2.5)) < 1e-12, ...
    'The lower edge was not the 2.5th percentile.');
assert(abs(band.high(1) - prctile(column, 97.5)) < 1e-12, ...
    'The upper edge was not the 97.5th percentile.');
assert(band.low(1) < band.median(1) && band.median(1) < band.high(1), ...
    'The shuffle band was not ordered from low to high.');
assert(isnan(band.low(2)) && isnan(band.median(2)) && isnan(band.high(2)), ...
    'A bin with no finite shuffle became a number.');
with_gap = session_acg_null_band([1 NaN 3; 3 NaN 9]);
assert(abs(with_gap.median(1) - 2) < 1e-12, ...
    'A finite column did not ignore the neighboring NaN.');
assert(isnan(with_gap.median(2)), ...
    'A NaN column was given a median.');
assert(abs(with_gap.median(3) - 6) < 1e-12, ...
    'The last shuffle column was not scored on its own.');
empty_band = session_acg_null_band([]);
assert(isempty(empty_band.low) && isempty(empty_band.high), ...
    'An empty shuffle became a band.');

root = [tempname '_session_acg'];
mkdir(root);
cleanup = onCleanup(@() rmdir(root, 's')); %#ok<NASGU>
lag = [0.1 0.3 0.6];
save_neuron(fullfile(root, 'SessionA'), 'S', 1, 'trough', [10 90 10], lag);
save_neuron(fullfile(root, 'SessionA'), 'S', 2, 'peak', [10 40 10], lag);
save_neuron(fullfile(root, 'SessionA'), 'S', 3, 'nonentrained', [5 5 5], lag);
save_neuron(fullfile(root, 'SessionB'), 'S', 1, 'trough', [1 1 1], lag);
save_neuron(fullfile(root, 'SessionC'), 'S', 4, 'trough', [80 80], [0.2 0.4]);
notes = fullfile(root, 'SessionA', 'neuron_notes.mat');
save(notes, 'lag');

records = collect_session_acg_percentiles(root);
assert(numel(records) == 4, 'A duplicate neuron was kept, or a neuron was dropped.');
kept = records([records.cluster_id] == 1);
assert(numel(kept) == 1 && abs(kept.curves.session.pctl(2) - 90) < 1e-12, ...
    'The first saved duplicate was not the curve that was kept.');
assert(kept.curves.session.oscillates, ...
    'The saved oscillation call was dropped.');
assert(isempty(kept.curves.fooof_high.pctl), ...
    'A tail that was not estimated gained a percentile curve.');

packed = session_acg_percentile_matrix(records, 'delta', 'session', 500, [0.2 0.5]);
assert(packed.n == 2, 'Delta-locked did not keep the trough and peak neurons.');
assert(packed.n_dropped == 1, 'The other lag axis stayed in the delta-locked matrix.');
assert(abs(packed.percentile(1, 2) - 90) < 1e-12, ...
    'The higher 200-500 ms percentile was not the first row.');
assert(abs(packed.percentile(2, 2) - 40) < 1e-12, ...
    'The peak-locked neuron was not the second delta-locked row.');
assert(all(isfinite(packed.logit(1, :))), ...
    'A finite percentile made the logit infinite.');
assert(packed.logit(1, 2) > 0, ...
    'A percentile above 50 had a negative logit.');

trough = session_acg_percentile_matrix(records, 'trough', 'session', 500, [0.2 0.5]);
assert(trough.n == 1 && trough.n_dropped == 1, ...
    'The trough-locked matrix mixed two lag axes.');
nonentrained = session_acg_percentile_matrix(records, 'nonentrained', 'session', 500, [0.2 0.5]);
assert(nonentrained.n == 1 && abs(nonentrained.percentile(1, 1) - 5) < 1e-12, ...
    'The non-entrained neuron was left out or moved.');
empty_tail = session_acg_percentile_matrix(records, 'all', 'fooof_high', 500, [0.2 0.5]);
assert(empty_tail.n == 0 && empty_tail.n_dropped == 0, ...
    'An unestimated tail still produced rows.');
all_session = session_acg_percentile_matrix(records, 'all', 'session', 500, [0.2 0.5]);
assert(all_session.n == 3 && all_session.n_dropped == 1, ...
    'The shared lag axis was not the one drawn for all neurons.');

disp('test_analyze_session_acg passed')
end


function save_neuron(folder, session_name, cluster_id, lock_group, pctl, lag)
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
neuron.session = session_name;
neuron.cluster_id = cluster_id;
neuron.lock_group = lock_group;
neuron.scores.session.pctl = pctl;
neuron.scores.session.pctl_lag = lag;
neuron.scores.session.oscillates = cluster_id == 1;
neuron.scores.session.first_lag = 0.3;
neuron.scores.session.first_p = 0.01;
save(fullfile(folder, sprintf('neuron_%d.mat', cluster_id)), 'neuron');
end
