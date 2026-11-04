%% Frequency-resolved theta phase-amplitude coupling
% This is a replacement for the broad 0.1-6 Hz PAC analysis.
%
% Objectives
%   1. Test whether mid-LPAG theta power is coupled to any slower LFP phase
%      without defining a delta band in advance.
%   2. If coupling is present, test coupling separately during annotated
%      play and non-play and directly test the play-minus-non-play
%      interaction. "Significant in play but not in non-play" is not treated
%      as evidence for a difference unless the interaction is significant.
%   3. Report a hierarchy of evidence:
%        a) session x LPAG-channel units;
%        b) an intermediate pooled-unit test (sensitive, not population
%           inference because units from one animal are not independent);
%        c) each animal after averaging its units;
%        d) an equal-animal group test (most conservative/generalizable).
%
% Method
%   - Theta power is abs(hilbert(6-12 Hz LFP)).^2 after downsampling.
%   - Candidate phase frequencies are scanned with seven-cycle Morlet
%     wavelets from 0.5:0.25:6 Hz. Frequencies whose upper two-sigma
%     wavelet support reaches the 6-12 Hz theta band are explicitly marked
%     as overlapping/exploratory; they cannot establish independent
%     phase-amplitude coupling from this same signal.
%   - Two maximum-statistic families are reported: primary inference over
%     0.5-4.5 Hz and an exploratory full-scan inference over 0.5-6 Hz.
%   - PAC is the magnitude of the first circular harmonic of z-scored
%     log-theta power. Its null independently circularly shifts theta power
%     inside contiguous 20 s blocks. Each block's autocorrelation and
%     spectrum are preserved, while independent block offsets prevent a
%     stationary coupling from surviving as only a changed preferred phase.
%   - Pointwise and maximum-statistic family-wise corrected p-values are
%     saved. All headline p-values are maximum-statistic corrected over the
%     complete modulation-frequency scan.
%   - Play and non-play coupling use the same circular null. Their direct
%     difference pairs play/non-play bouts by duration and session time,
%     truncates each pair to equal duration, and swaps labels within pairs.
%     This preserves within-bout structure and exactly matches condition
%     duration in every permutation.
%   - A harmonic control tests n:m phase locking between candidate slow
%     phase and theta carrier phase. Significant PAC together with
%     significant harmonic locking remains ambiguous with a non-sinusoidal
%     slow waveform; strong visible theta alone does not remove that issue.
%
% The script checkpoints after every recording and resumes by default.
% It is intended to be left running overnight. Progress is printed for every
% recording, LPAG channel, partner session, and fifth scan frequency.

%% Options
options = struct();
options.data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data'];
options.sampling_rate = 2500;
options.analysis_sampling_rate = 250;
options.n_lfp_channels = 384;
options.area_name = 'LPAG';

options.theta_band_hz = [6 12];
options.theta_filter_seconds = 4;
options.theta_measure = 'log_power'; % 'log_power' or 'linear_power'
options.modulation_frequencies_hz = 0.5:0.25:6;
options.primary_inference_max_hz = 4.5;
options.morlet_cycles = 7;
options.morlet_half_width_sd = 3.5;
options.edge_trim_seconds = 12;

options.n_permutations = 1000;
options.n_group_permutations = 5000;
options.surrogate_block_seconds = 20;
options.minimum_surrogate_blocks = 3;
options.minimum_within_block_shift_seconds = 1;
options.minimum_condition_seconds = 20;
options.minimum_condition_bouts = 5;
options.minimum_animals_group = 3;
options.maximum_harmonic = 24;
options.random_seed = 1907;

options.play_behaviors = ...
    {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD'};
options.nonplay_behaviors = ...
    {'Grooming', 'Rearing', 'Sniffing', 'Scratching', 'PounceI', 'Bite'};
options.bout_merge_gap_seconds = 1;

options.resume = true;
options.overwrite = false;
options.save_figures = true;
options.output_subfolder = fullfile( ...
    'Analysis results', 'frequency resolved theta PAC');
options.output_filename = 'frequency_resolved_theta_pac.mat';
% Add a recording folder name here only after documenting why it should not
% enter the population analysis. Exclusions are saved in the signature.
options.recordings_to_exclude = {};

% For a short pipeline check, set true. These results are deliberately saved
% under another filename and must not be used for inference.
options.quick_test = false;
if options.quick_test
    options.n_permutations = 20;
    options.n_group_permutations = 100;
    options.modulation_frequencies_hz = 0.5:0.5:6;
    options.output_filename = 'frequency_resolved_theta_pac_QUICK_TEST.mat';
end

%% Paths and checkpoint
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
supporting_code_dir = fileparts(script_dir);
repo_root = fileparts(fileparts(supporting_code_dir));
addpath(script_dir);
addpath(supporting_code_dir);
addpath(fullfile(repo_root, 'Custom functions'));

raw_data_dir = fullfile(options.data_root, 'NPX data', 'NPX raw data');
output_dir = fullfile(options.data_root, options.output_subfolder);
if ~isfolder(output_dir)
    mkdir(output_dir);
end
output_file = fullfile(output_dir, options.output_filename);
log_file = fullfile(output_dir, ...
    strrep(options.output_filename, '.mat', '_progress.txt'));

validate_options(options);
signature = analysis_signature(options);
pac_units = empty_pac_units();
processed_recordings = {};
failures = empty_failures();
hierarchy = struct();
unit_table = table();
animal_table = table();
group_table = table();

if options.resume && ~options.overwrite && isfile(output_file)
    previous = load(output_file, 'signature', 'pac_units', ...
        'processed_recordings', 'failures');
    if ~isfield(previous, 'signature') || ~isequaln(previous.signature, signature)
        error(['The checkpoint was created with different analysis settings. ' ...
            'Rename it, set options.overwrite=true, or restore the old settings.']);
    end
    if isfield(previous, 'pac_units')
        pac_units = previous.pac_units;
    end
    if isfield(previous, 'processed_recordings')
        processed_recordings = previous.processed_recordings;
    end
    if isfield(previous, 'failures')
        failures = previous.failures;
    end
    fprintf('Resuming %s with %d completed recording(s) and %d PAC unit(s).\n', ...
        output_file, numel(processed_recordings), numel(pac_units));
elseif options.overwrite && isfile(output_file)
    delete(output_file);
end

diary(log_file);
diary on
cleanup_diary = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('\nFrequency-resolved theta PAC started %s\n', datestr(now, 31));
fprintf('Checkpoint: %s\n', output_file);

%% Analyze every recording
recording_list = dir(raw_data_dir);
recording_list = recording_list([recording_list.isdir]);
recording_list = recording_list(~ismember({recording_list.name}, {'.', '..'}));
has_lfp = arrayfun(@(x) ...
    isfile(fullfile(raw_data_dir, x.name, 'LFP_PAG.mat')) || ...
    isfile(fullfile(raw_data_dir, x.name, 'LFP_PAG.dat')), recording_list);
recording_list = recording_list(has_lfp);
recording_list = recording_list(~ismember( ...
    {recording_list.name}, options.recordings_to_exclude));
[~, recording_order] = sort(lower(string({recording_list.name})));
recording_list = recording_list(recording_order);

for recording_n = 1:numel(recording_list)
    recording_name = recording_list(recording_n).name;
    if ismember(recording_name, processed_recordings)
        fprintf('Skipping completed recording %s (%d/%d).\n', ...
            recording_name, recording_n, numel(recording_list));
        continue
    end

    fprintf('\n=== Recording %s (%d/%d) ===\n', ...
        recording_name, recording_n, numel(recording_list));
    recording_dir = fullfile(raw_data_dir, recording_name);
    if ~isempty(pac_units)
        pac_units(ismember({pac_units.recording}, recording_name)) = [];
    end
    if ~isempty(failures)
        failures(ismember({failures.recording}, recording_name)) = [];
    end
    try
        [new_units, new_failures] = analyze_recording( ...
            recording_dir, options);
        if isempty(pac_units)
            pac_units = new_units(:);
        elseif ~isempty(new_units)
            pac_units = [pac_units(:); new_units(:)]; %#ok<AGROW>
        end
        if ~isempty(new_failures)
            failures = [failures(:); new_failures(:)]; %#ok<AGROW>
            fprintf(['Recording %s remains incomplete because %d unit(s) ' ...
                'failed; a resumed run will retry it.\n'], ...
                recording_name, numel(new_failures));
        else
            processed_recordings{end + 1, 1} = recording_name; %#ok<SAGROW>
        end
        fprintf('Finished %s: %d valid session-channel unit(s).\n', ...
            recording_name, numel(new_units));
    catch exception
        warning('FrequencyPAC:RecordingFailed', ...
            'Recording %s failed: %s', recording_name, exception.message);
        failure = struct( ...
            'recording', recording_name, ...
            'channel', NaN, ...
            'session', NaN, ...
            'identifier', exception.identifier, ...
            'message', exception.message, ...
            'report', getReport(exception, 'extended', 'hyperlinks', 'off'), ...
            'time', datestr(now, 31));
        if isempty(failures)
            failures = failure;
        else
            failures(end + 1, 1) = failure; %#ok<SAGROW>
        end
    end

    checkpoint = struct( ...
        'signature', signature, ...
        'options', options, ...
        'pac_units', pac_units, ...
        'processed_recordings', {processed_recordings}, ...
        'failures', failures, ...
        'hierarchy', hierarchy, ...
        'unit_table', unit_table, ...
        'animal_table', animal_table, ...
        'group_table', group_table);
    save_checkpoint(output_file, checkpoint);
end

all_recordings = {recording_list.name};
incomplete_recordings = setdiff( ...
    all_recordings, processed_recordings, 'stable');
if ~isempty(incomplete_recordings)
    error(['The overnight estimate is incomplete for: %s. The checkpoint ' ...
        'is safe; inspect failures and rerun after correcting transient or ' ...
        'data-specific errors. Group inference was intentionally not run.'], ...
        strjoin(incomplete_recordings, ', '));
end

if isempty(pac_units)
    error('No valid PAC units were generated. See %s.', log_file);
end

%% Hierarchical inference and output
fprintf('\nCombining %d PAC units through the inference hierarchy...\n', ...
    numel(pac_units));
hierarchy = build_hierarchy(pac_units, options);
[unit_table, animal_table, group_table] = result_tables( ...
    pac_units, hierarchy);

writetable(unit_table, fullfile(output_dir, ...
    'frequency_resolved_theta_pac_units.csv'));
writetable(animal_table, fullfile(output_dir, ...
    'frequency_resolved_theta_pac_animals.csv'));
writetable(group_table, fullfile(output_dir, ...
    'frequency_resolved_theta_pac_group.csv'));

print_summary(pac_units, hierarchy);
if options.save_figures
    make_summary_figures(pac_units, hierarchy, output_dir);
end

checkpoint = struct( ...
    'signature', signature, ...
    'options', options, ...
    'pac_units', pac_units, ...
    'processed_recordings', {processed_recordings}, ...
    'failures', failures, ...
    'hierarchy', hierarchy, ...
    'unit_table', unit_table, ...
    'animal_table', animal_table, ...
    'group_table', group_table);
save_checkpoint(output_file, checkpoint);
fprintf('Complete. Results saved to %s\n', output_file);
fprintf('Finished %s\n', datestr(now, 31));

%% Local functions
function validate_options(options)
if options.analysis_sampling_rate >= options.sampling_rate
    error('analysis_sampling_rate must be below the acquisition sampling rate.');
end
if options.theta_band_hz(2) >= options.analysis_sampling_rate / 2
    error('The theta band must lie below the downsampled Nyquist frequency.');
end
if any(options.modulation_frequencies_hz <= 0) || ...
        any(options.modulation_frequencies_hz > options.theta_band_hz(1))
    error('Modulation frequencies must lie above zero and at/below theta onset.');
end
if options.primary_inference_max_hz < ...
        min(options.modulation_frequencies_hz) || ...
        options.primary_inference_max_hz > ...
        max(options.modulation_frequencies_hz)
    error('primary_inference_max_hz must lie inside the scan.');
end
upper_two_sigma = max(options.modulation_frequencies_hz) * ...
    (1 + 2 / options.morlet_cycles);
if upper_two_sigma >= options.theta_band_hz(1)
    overlap_start = options.theta_band_hz(1) / ...
        (1 + 2 / options.morlet_cycles);
    warning('FrequencyPAC:PhaseThetaOverlap', ...
        ['Phase frequencies from %.2f Hz upward overlap theta within two ' ...
        'wavelet spectral SD. They will be computed but flagged as ' ...
        'exploratory/ambiguous.'], overlap_start);
end
if options.n_permutations < 19 || options.n_group_permutations < 99
    error('Use at least 19 unit and 99 group permutations.');
end
if options.surrogate_block_seconds <= ...
        2 * options.minimum_within_block_shift_seconds
    error('Surrogate blocks must exceed twice the minimum shift.');
end
if ~ismember(options.theta_measure, {'log_power', 'linear_power'})
    error('theta_measure must be log_power or linear_power.');
end
end

function signature = analysis_signature(options)
signature = struct( ...
    'version', 5, ...
    'code_version', 'frequency-resolved-pac-2026-11-04-v5', ...
    'sampling_rate', options.sampling_rate, ...
    'analysis_sampling_rate', options.analysis_sampling_rate, ...
    'n_lfp_channels', options.n_lfp_channels, ...
    'area_name', options.area_name, ...
    'theta_band_hz', options.theta_band_hz, ...
    'theta_filter_seconds', options.theta_filter_seconds, ...
    'theta_measure', options.theta_measure, ...
    'modulation_frequencies_hz', options.modulation_frequencies_hz, ...
    'primary_inference_max_hz', options.primary_inference_max_hz, ...
    'morlet_cycles', options.morlet_cycles, ...
    'morlet_half_width_sd', options.morlet_half_width_sd, ...
    'edge_trim_seconds', options.edge_trim_seconds, ...
    'n_permutations', options.n_permutations, ...
    'n_group_permutations', options.n_group_permutations, ...
    'surrogate_block_seconds', options.surrogate_block_seconds, ...
    'minimum_surrogate_blocks', options.minimum_surrogate_blocks, ...
    'minimum_within_block_shift_seconds', ...
    options.minimum_within_block_shift_seconds, ...
    'minimum_condition_seconds', options.minimum_condition_seconds, ...
    'minimum_condition_bouts', options.minimum_condition_bouts, ...
    'minimum_animals_group', options.minimum_animals_group, ...
    'maximum_harmonic', options.maximum_harmonic, ...
    'random_seed', options.random_seed, ...
    'recordings_to_exclude', {options.recordings_to_exclude}, ...
    'play_behaviors', {options.play_behaviors}, ...
    'nonplay_behaviors', {options.nonplay_behaviors}, ...
    'bout_merge_gap_seconds', options.bout_merge_gap_seconds);
end

function frequency = phase_theta_overlap_start(options)
frequency = options.theta_band_hz(1) / ...
    (1 + 2 / options.morlet_cycles);
end

function [units, local_failures] = analyze_recording(recording_dir, options)
[recording_name, animal_id, repeated_animal] = recording_names(recording_dir);
local_failures = empty_failures();
[behavior, partner_intervals] = load_behavior(recording_name, options);
if isempty(partner_intervals)
    warning('FrequencyPAC:NoPartnerSession', ...
        '%s has no Partners session intervals.', recording_name);
    units = empty_pac_units();
    return
end

play_intervals = behavior_intervals(behavior, ...
    options.play_behaviors, options.bout_merge_gap_seconds);
nonplay_intervals = behavior_intervals(behavior, ...
    options.nonplay_behaviors, options.bout_merge_gap_seconds);

[xcoords, ycoords, chan_map] = load_channel_map(recording_dir);
npx_type = neuropixels_type(xcoords);
area_table = readtable(fullfile(options.data_root, ...
    'Area_limits_GoodLooking.xlsx'));
area_table = area_table(ismember(area_table.AnimalName, animal_id), :);
lpag = area_table(ismember(area_table.area, options.area_name), :);
if isempty(lpag)
    error('No %s area limits found for %s.', options.area_name, animal_id);
end
limits = lpag{:, {'ProbeNum', 'depth_start', 'depth_end'}};
selected = select_mid_lpag_channels( ...
    xcoords, ycoords, chan_map, limits, npx_type);
[lfp, n_samples] = read_selected_lfp(recording_dir, ...
    [selected.lfp_row], options.n_lfp_channels);

recording_end = n_samples / options.sampling_rate;
partner_intervals(:, 1) = max(partner_intervals(:, 1), 0);
partner_intervals(:, 2) = min(partner_intervals(:, 2), recording_end);
partner_intervals = partner_intervals( ...
    partner_intervals(:, 2) > partner_intervals(:, 1), :);

units = empty_pac_units();
for channel_n = 1:numel(selected)
    for session_n = 1:size(partner_intervals, 1)
        bounds = partner_intervals(session_n, :);
        fprintf('%s | LPAG channel %g (%d/%d) | partner session %d/%d\n', ...
            recording_name, selected(channel_n).lfp_row, ...
            channel_n, numel(selected), session_n, size(partner_intervals, 1));
        try
            unit = analyze_session( ...
                lfp(channel_n, :), bounds, play_intervals, ...
                nonplay_intervals, options, ...
                deterministic_seed(options.random_seed, ...
                recording_name, session_n));
            unit.animal_id = animal_id;
            unit.recording = recording_name;
            unit.repeated_animal = repeated_animal;
            unit.lpag_channel = selected(channel_n).channel;
            unit.lpag_lfp_row = selected(channel_n).lfp_row;
            unit.probe_number = selected(channel_n).probe_number;
            unit.session_index = session_n;
            unit.session_bounds_seconds = bounds;
            if isempty(units)
                units = unit;
            else
                units(end + 1, 1) = unit; %#ok<AGROW>
            end
        catch exception
            warning('FrequencyPAC:SessionSkipped', ...
                '%s channel %g session %d skipped: %s', ...
                recording_name, selected(channel_n).lfp_row, ...
                session_n, exception.message);
            failure = struct( ...
                'recording', recording_name, ...
                'channel', selected(channel_n).lfp_row, ...
                'session', session_n, ...
                'identifier', exception.identifier, ...
                'message', exception.message, ...
                'report', getReport(exception, 'extended', ...
                'hyperlinks', 'off'), ...
                'time', datestr(now, 31));
            if isempty(local_failures)
                local_failures = failure;
            else
                local_failures(end + 1, 1) = failure; %#ok<AGROW>
            end
        end
    end
end
end

function unit = analyze_session(raw_lfp, bounds, play_intervals, ...
        nonplay_intervals, options, unit_seed)
rng(unit_seed, 'twister');
fs = options.sampling_rate;
target_fs = options.analysis_sampling_rate;
first_sample = max(1, floor(bounds(1) * fs) + 1);
last_sample = min(numel(raw_lfp), ceil(bounds(2) * fs));
segment = double(raw_lfp(first_sample:last_sample));
segment = detrend(segment(:), 'linear');

[resample_num, resample_den] = rat(target_fs / fs);
segment = resample(segment, resample_num, resample_den);
segment = segment(:);
time = (first_sample - 1) / fs + (0:numel(segment) - 1)' / target_fs;
wavelet_signal = segment;

theta_order = round(options.theta_filter_seconds * target_fs);
if mod(theta_order, 2) == 1
    theta_order = theta_order + 1;
end
if numel(segment) <= 3 * theta_order + 1
    error('Partner session is too short for the theta filter.');
end
theta_filter = fir1(theta_order, ...
    options.theta_band_hz / (target_fs / 2), 'bandpass', ...
    hamming(theta_order + 1));
theta_lfp = filtfilt(theta_filter, 1, segment);
theta_analytic = hilbert(theta_lfp);
theta_power = abs(theta_analytic) .^ 2;
if strcmp(options.theta_measure, 'log_power')
    positive_power = theta_power(theta_power > 0 & isfinite(theta_power));
    if isempty(positive_power)
        error('Theta power contains no positive finite values.');
    end
    theta_measure = log10(theta_power + median(positive_power) * eps);
else
    theta_measure = theta_power;
end

lowest_frequency = min(options.modulation_frequencies_hz);
max_wave_half_seconds = options.morlet_half_width_sd * ...
    options.morlet_cycles / (2 * pi * lowest_frequency);
edge_seconds = max(options.edge_trim_seconds, max_wave_half_seconds);
valid = time >= time(1) + edge_seconds & time <= time(end) - edge_seconds;
minimum_required = options.minimum_surrogate_blocks * ...
    options.surrogate_block_seconds * target_fs;
if sum(valid) < minimum_required
    error('Too little data remains after the common filter/wavelet edge trim.');
end
surrogate_block_samples = round( ...
    options.surrogate_block_seconds * target_fs);
n_complete_blocks = floor(sum(valid) / surrogate_block_samples);
n_usable = n_complete_blocks * surrogate_block_samples;
wavelet_valid_index = find(valid);
wavelet_valid_index = wavelet_valid_index(1:n_usable);
time = time(wavelet_valid_index);
theta_measure = theta_measure(wavelet_valid_index);

play_mask = interval_mask(time, play_intervals);
nonplay_mask = interval_mask(time, nonplay_intervals);
overlap = play_mask & nonplay_mask;
play_mask(overlap) = false;
nonplay_mask(overlap) = false;
all_mask = true(size(time));

n_samples = numel(time);
surrogate_blocks = surrogate_block_partition(n_samples, ...
    surrogate_block_samples, ...
    options.minimum_surrogate_blocks);
block_lags = independent_block_lags(surrogate_blocks, ...
    options.n_permutations, round( ...
    options.minimum_within_block_shift_seconds * target_fs));
all_cache = pac_normalization_cache( ...
    theta_measure, all_mask, surrogate_blocks, block_lags);
play_valid = sum(play_mask) / target_fs >= ...
    options.minimum_condition_seconds;
nonplay_valid = sum(nonplay_mask) / target_fs >= ...
    options.minimum_condition_seconds;
if play_valid
    play_cache = pac_normalization_cache( ...
        theta_measure, play_mask, surrogate_blocks, block_lags);
else
    play_cache = struct();
end
if nonplay_valid
    nonplay_cache = pac_normalization_cache( ...
        theta_measure, nonplay_mask, surrogate_blocks, block_lags);
else
    nonplay_cache = struct();
end

[blocks, block_labels] = matched_condition_bouts( ...
    play_mask, nonplay_mask);
contrast_valid = play_valid && nonplay_valid && ...
    sum(block_labels) >= options.minimum_condition_bouts && ...
    sum(~block_labels) >= options.minimum_condition_bouts;
if contrast_valid
    block_permutations = permuted_paired_labels( ...
        sum(block_labels), options.n_permutations);
    block_cache = block_normalization_cache(theta_measure, blocks);
else
    block_permutations = false(0, 0);
    block_cache = struct();
end

frequencies = options.modulation_frequencies_hz(:)';
overlap_start_hz = phase_theta_overlap_start(options);
n_frequency = numel(frequencies);
n_draw = options.n_permutations + 1;
whole_values = nan(n_draw, n_frequency);
play_values = nan(n_draw, n_frequency);
nonplay_values = nan(n_draw, n_frequency);
contrast_values = nan(n_draw, n_frequency);
if contrast_valid
    contrast_values = nan(size(block_permutations, 1), n_frequency);
end
harmonic_values = nan(n_draw, n_frequency);
preferred_phase_whole = nan(1, n_frequency);
preferred_phase_play = nan(1, n_frequency);
preferred_phase_nonplay = nan(1, n_frequency);
slow_wavelet_power = nan(1, n_frequency);

max_wave_half_samples = ceil(max_wave_half_seconds * target_fs);
maximum_wave_length = 2 * max_wave_half_samples + 1;
n_fft_convolution = 2 ^ nextpow2( ...
    numel(wavelet_signal) + maximum_wave_length - 1);
segment_fft = fft(wavelet_signal, n_fft_convolution);
[theta_harmonic_frequency, theta_harmonic_fft_blocks] = ...
    harmonic_phase_transforms(segment_fft, numel(wavelet_signal), ...
    n_fft_convolution, target_fs, wavelet_valid_index, ...
    surrogate_blocks, frequencies, options.theta_band_hz, ...
    options.maximum_harmonic, options.morlet_cycles, ...
    options.morlet_half_width_sd);

for frequency_n = 1:n_frequency
    frequency = frequencies(frequency_n);
    [phase_full, wavelet_power_full] = morlet_phase_from_fft( ...
        segment_fft, numel(wavelet_signal), n_fft_convolution, ...
        target_fs, frequency, options.morlet_cycles, ...
        options.morlet_half_width_sd);
    phase = phase_full(wavelet_valid_index);
    slow_wavelet_power(frequency_n) = ...
        mean(wavelet_power_full(wavelet_valid_index), 'omitnan');
    phase_phasor = exp(1i * phase);

    [whole_values(:, frequency_n), preferred_phase_whole(frequency_n)] = ...
        shifted_pac(phase_phasor, all_cache);
    if play_valid
        [play_values(:, frequency_n), preferred_phase_play(frequency_n)] = ...
            shifted_pac(phase_phasor, play_cache);
    end
    if nonplay_valid
        [nonplay_values(:, frequency_n), preferred_phase_nonplay(frequency_n)] = ...
            shifted_pac(phase_phasor, nonplay_cache);
    end
    if contrast_valid
        contrast_values(:, frequency_n) = block_label_contrast( ...
            theta_measure, phase_phasor, blocks, block_cache, ...
            block_permutations);
    end
    harmonic_values(:, frequency_n) = harmonic_shift_control( ...
        phase_phasor, theta_harmonic_frequency, ...
        theta_harmonic_fft_blocks, all_mask, ...
        surrogate_blocks, block_lags, frequency, ...
        options.theta_band_hz, options.maximum_harmonic);

    if frequency_n == 1 || mod(frequency_n, 5) == 0 || ...
            frequency_n == n_frequency
        fprintf('  scan %d/%d: %.2f Hz\n', ...
            frequency_n, n_frequency, frequency);
    end
end

unit = empty_unit();
unit.frequency_hz = frequencies;
unit.phase_theta_overlap = frequencies >= overlap_start_hz;
unit.phase_theta_overlap_start_hz = overlap_start_hz;
unit.analysis_duration_seconds = n_samples / target_fs;
unit.play_seconds = sum(play_mask) / target_fs;
unit.nonplay_seconds = sum(nonplay_mask) / target_fs;
unit.overlap_removed_seconds = sum(overlap) / target_fs;
unit.n_play_blocks = sum(block_labels);
unit.n_nonplay_blocks = sum(~block_labels);
unit.edge_trim_seconds = edge_seconds;
unit.surrogate_block_shifts_samples = block_lags(2:end, :);
unit.whole = score_against_null(whole_values, frequencies, false, ...
    options.primary_inference_max_hz);
unit.play = score_against_null(play_values, frequencies, false, ...
    options.primary_inference_max_hz);
unit.nonplay = score_against_null(nonplay_values, frequencies, false, ...
    options.primary_inference_max_hz);
unit.contrast = score_assignment_permutations( ...
    contrast_values, 1, frequencies, true, ...
    options.primary_inference_max_hz);
unit.harmonic = score_against_null(harmonic_values, frequencies, false, ...
    options.primary_inference_max_hz);
unit.preferred_phase_whole = preferred_phase_whole;
unit.preferred_phase_play = preferred_phase_play;
unit.preferred_phase_nonplay = preferred_phase_nonplay;
unit.slow_wavelet_power = slow_wavelet_power;
unit.notes = ['PAC null: independent intact circular shifts of theta power ' ...
    'inside contiguous surrogate blocks. ' ...
    'Condition contrast null: paired duration-matched play/non-play bouts ' ...
    'with within-pair label swaps. Headline p-values use max-statistic ' ...
    'FWER correction.'];
end

function cache = pac_normalization_cache( ...
        amplitude, mask, blocks, block_lags)
mask = logical(mask(:));
amplitude = double(amplitude(:));
n = sum(mask);
if n < 3
    error('PAC mask has fewer than three samples.');
end
n_draw = size(block_lags, 1);
sum_amplitude = zeros(n_draw, 1);
sum_square = zeros(n_draw, 1);
amplitude_fft = cell(numel(blocks), 1);
for block_n = 1:numel(blocks)
    index = blocks{block_n};
    block_amplitude = amplitude(index);
    block_mask = mask(index);
    amplitude_fft{block_n} = fft(block_amplitude);
    this_sum = circular_weighted_sum_fft( ...
        block_mask, amplitude_fft{block_n});
    this_square = circular_weighted_sum_fft( ...
        block_mask, fft(block_amplitude .^ 2));
    selected = block_lags(:, block_n) + 1;
    sum_amplitude = sum_amplitude + real(this_sum(selected));
    sum_square = sum_square + real(this_square(selected));
end
mean_amplitude = sum_amplitude / n;
variance = (sum_square - (sum_amplitude .^ 2) / n) / (n - 1);
variance = max(real(variance), 0);
cache = struct( ...
    'mask', mask, ...
    'blocks', {blocks}, ...
    'block_lags', block_lags, ...
    'n', n, ...
    'mean', mean_amplitude, ...
    'standard_deviation', sqrt(variance), ...
    'amplitude_fft', {amplitude_fft});
end

function [values, preferred_phase] = shifted_pac( ...
        phase_phasor, cache)
cross_sum = complex(zeros(size(cache.block_lags, 1), 1));
sum_phase = 0;
for block_n = 1:numel(cache.blocks)
    index = cache.blocks{block_n};
    weight = double(cache.mask(index)) .* phase_phasor(index);
    this_cross = circular_weighted_sum_fft( ...
        weight, cache.amplitude_fft{block_n});
    selected = cache.block_lags(:, block_n) + 1;
    cross_sum = cross_sum + this_cross(selected);
    sum_phase = sum_phase + sum(weight);
end
centered_sum = cross_sum - cache.mean .* sum_phase;
denominator = cache.n .* cache.standard_deviation;
values = abs(centered_sum) ./ denominator;
values(~isfinite(values) | denominator <= 0) = NaN;
preferred_phase = angle(centered_sum(1));
end

function sums = circular_weighted_sum_fft(weight, signal_fft)
% sums(k+1) = sum_t weight(t) * signal(mod(t+k-1,N)+1).
weight = weight(:);
if numel(weight) ~= numel(signal_fft)
    error('Circular-correlation inputs must have equal length.');
end
sums = ifft(conj(fft(conj(weight))) .* signal_fft);
end

function blocks = surrogate_block_partition( ...
        n_samples, block_samples, minimum_blocks)
if mod(n_samples, block_samples) ~= 0
    error('Surrogate analysis samples must divide into complete blocks.');
end
n_blocks = floor(n_samples / block_samples);
if n_blocks < minimum_blocks
    error('Too few complete surrogate blocks remain after edge trimming.');
end
blocks = cell(n_blocks, 1);
for block_n = 1:n_blocks
    first = (block_n - 1) * block_samples + 1;
    last = block_n * block_samples;
    blocks{block_n} = first:last;
end
end

function lags = independent_block_lags( ...
        blocks, n_permutations, minimum_shift)
lags = zeros(n_permutations + 1, numel(blocks));
for block_n = 1:numel(blocks)
    n = numel(blocks{block_n});
    eligible = minimum_shift:(n - minimum_shift);
    if isempty(eligible)
        error('A surrogate block is too short for the minimum shift.');
    end
    lags(2:end, block_n) = eligible(randi( ...
        numel(eligible), n_permutations, 1));
end
end

function transforms = block_signal_fft(signal, blocks)
transforms = cell(numel(blocks), 1);
for block_n = 1:numel(blocks)
    transforms{block_n} = fft(signal(blocks{block_n}));
end
end

function [harmonic_frequency, transforms] = harmonic_phase_transforms( ...
        signal_fft, n_signal, n_fft, sampling_rate, valid_index, blocks, ...
        slow_frequencies, theta_band, maximum_harmonic, n_cycles, ...
        half_width_sd)
harmonic_frequency = [];
for slow_frequency = slow_frequencies
    candidate = (2:maximum_harmonic) * slow_frequency;
    candidate = candidate(candidate >= theta_band(1) & ...
        candidate <= theta_band(2));
    harmonic_frequency = [harmonic_frequency candidate]; %#ok<AGROW>
end
harmonic_frequency = unique(round(harmonic_frequency * 1e9) / 1e9);
transforms = cell(numel(harmonic_frequency), 1);
for frequency_n = 1:numel(harmonic_frequency)
    phase_full = morlet_phase_from_fft( ...
        signal_fft, n_signal, n_fft, sampling_rate, ...
        harmonic_frequency(frequency_n), n_cycles, half_width_sd);
    phasor = exp(-1i * phase_full(valid_index));
    transforms{frequency_n} = block_signal_fft(phasor, blocks);
end
end

function [phase, power] = morlet_phase_from_fft(signal_fft, n_signal, n_fft, ...
        sampling_rate, frequency, n_cycles, half_width_sd)
sigma_seconds = n_cycles / (2 * pi * frequency);
half_samples = ceil(half_width_sd * sigma_seconds * sampling_rate);
wave_time = (-half_samples:half_samples)' / sampling_rate;
correction = exp(-0.5 * n_cycles ^ 2);
wavelet = (exp(2i * pi * frequency * wave_time) - correction) .* ...
    exp(-(wave_time .^ 2) / (2 * sigma_seconds ^ 2));
wavelet = wavelet / sqrt(sum(abs(wavelet) .^ 2));
convolution = ifft(signal_fft .* fft(wavelet, n_fft));
first = half_samples + 1;
analytic = convolution(first:first + n_signal - 1);
phase = angle(analytic);
power = abs(analytic) .^ 2;
end

function [blocks, labels] = matched_condition_bouts( ...
        play_mask, nonplay_mask)
play_blocks = contiguous_runs(play_mask);
nonplay_blocks = contiguous_runs(nonplay_mask);
matched_play = cell(0, 1);
matched_nonplay = cell(0, 1);
available_play = 1:numel(play_blocks);
available_nonplay = 1:numel(nonplay_blocks);
n_session_samples = numel(play_mask);
while ~isempty(available_play) && ~isempty(available_nonplay)
    cost = nan(numel(available_play), numel(available_nonplay));
    for play_n = 1:numel(available_play)
        play_index = play_blocks{available_play(play_n)};
        for nonplay_n = 1:numel(available_nonplay)
            nonplay_index = nonplay_blocks{available_nonplay(nonplay_n)};
            duration_cost = abs(log(numel(play_index) / ...
                numel(nonplay_index)));
            time_cost = abs(mean(play_index) - mean(nonplay_index)) / ...
                n_session_samples;
            cost(play_n, nonplay_n) = duration_cost + time_cost;
        end
    end
    [~, linear_index] = min(cost(:));
    [play_choice, nonplay_choice] = ind2sub(size(cost), linear_index);
    play_index = play_blocks{available_play(play_choice)};
    nonplay_index = nonplay_blocks{available_nonplay(nonplay_choice)};
    matched_length = min(numel(play_index), numel(nonplay_index));
    matched_play{end + 1, 1} = center_crop(play_index, matched_length); %#ok<AGROW>
    matched_nonplay{end + 1, 1} = ...
        center_crop(nonplay_index, matched_length); %#ok<AGROW>
    available_play(play_choice) = [];
    available_nonplay(nonplay_choice) = [];
end
blocks = [matched_play; matched_nonplay];
labels = [true(numel(matched_play), 1); ...
    false(numel(matched_nonplay), 1)];
end

function cropped = center_crop(index, target_length)
first = floor((numel(index) - target_length) / 2) + 1;
cropped = index(first:first + target_length - 1);
end

function blocks = contiguous_runs(mask)
mask = logical(mask(:));
change = diff([false; mask; false]);
run_start = find(change == 1);
run_end = find(change == -1) - 1;
blocks = cell(0, 1);
for run_n = 1:numel(run_start)
    blocks{end + 1, 1} = run_start(run_n):run_end(run_n); %#ok<AGROW>
end
end

function permutations = permuted_paired_labels(n_pairs, n_permutations)
if n_pairs <= 16
    n_assignment = 2 ^ n_pairs;
    codes = uint64((0:n_assignment - 2)');
    first_is_play = false(n_assignment - 1, n_pairs);
    for pair_n = 1:n_pairs
        first_is_play(:, pair_n) = logical(bitget(codes, pair_n));
    end
else
    first_is_play = rand(n_permutations, n_pairs) >= 0.5;
end
permutations = false(size(first_is_play, 1) + 1, 2 * n_pairs);
permutations(1, 1:n_pairs) = true;
permutations(2:end, 1:n_pairs) = first_is_play;
permutations(2:end, n_pairs + 1:end) = ~first_is_play;
end

function cache = block_normalization_cache(amplitude, blocks)
n_blocks = numel(blocks);
cache.n = zeros(n_blocks, 1);
cache.sum = zeros(n_blocks, 1);
cache.sum_square = zeros(n_blocks, 1);
for block_n = 1:n_blocks
    values = amplitude(blocks{block_n});
    cache.n(block_n) = numel(values);
    cache.sum(block_n) = sum(values);
    cache.sum_square(block_n) = sum(values .^ 2);
end
end

function difference = block_label_contrast(amplitude, phase_phasor, ...
        blocks, cache, permutations)
n_blocks = numel(blocks);
sum_phase = complex(zeros(n_blocks, 1));
sum_cross = complex(zeros(n_blocks, 1));
for block_n = 1:n_blocks
    index = blocks{block_n};
    sum_phase(block_n) = sum(phase_phasor(index));
    sum_cross(block_n) = sum(amplitude(index) .* phase_phasor(index));
end

label_double = double(permutations);
play = aggregate_block_pac(label_double, cache, sum_phase, sum_cross);
nonplay = aggregate_block_pac(1 - label_double, cache, ...
    sum_phase, sum_cross);
difference = play - nonplay;
end

function values = aggregate_block_pac(label_matrix, cache, ...
        sum_phase, sum_cross)
n = label_matrix * cache.n;
sum_amplitude = label_matrix * cache.sum;
sum_square = label_matrix * cache.sum_square;
phase_sum = label_matrix * sum_phase;
cross_sum = label_matrix * sum_cross;
mean_amplitude = sum_amplitude ./ n;
variance = (sum_square - (sum_amplitude .^ 2) ./ n) ./ (n - 1);
standard_deviation = sqrt(max(real(variance), 0));
centered_sum = cross_sum - mean_amplitude .* phase_sum;
values = abs(centered_sum) ./ (n .* standard_deviation);
values(~isfinite(values) | n < 3 | standard_deviation <= 0) = NaN;
end

function values = harmonic_shift_control( ...
        slow_phasor, theta_harmonic_frequency, ...
        theta_harmonic_fft_blocks, mask, blocks, block_lags, ...
        slow_frequency, theta_band, maximum_harmonic)
candidate = 2:maximum_harmonic;
candidate = candidate(candidate * slow_frequency >= theta_band(1) & ...
    candidate * slow_frequency <= theta_band(2));
values = nan(size(block_lags, 1), 1);
if isempty(candidate)
    return
end
for harmonic = candidate
    target_frequency = harmonic * slow_frequency;
    [distance, frequency_index] = min(abs( ...
        theta_harmonic_frequency - target_frequency));
    if distance > 1e-8
        error('Missing harmonic phase transform at %.9g Hz.', ...
            target_frequency);
    end
    frequency_transforms = ...
        theta_harmonic_fft_blocks{frequency_index};
    phase_sum = complex(zeros(size(block_lags, 1), 1));
    for block_n = 1:numel(blocks)
        index = blocks{block_n};
        weight = double(mask(index)) .* ...
            (slow_phasor(index) .^ harmonic);
        this_sum = circular_weighted_sum_fft( ...
            weight, frequency_transforms{block_n});
        selected = block_lags(:, block_n) + 1;
        phase_sum = phase_sum + this_sum(selected);
    end
    this_value = abs(phase_sum) / sum(mask);
    if all(isnan(values))
        values = this_value;
    else
        values = max(values, this_value);
    end
end
end

function scored = score_against_null( ...
        values, frequencies, two_sided, primary_max_hz)
if nargin < 4 || isempty(primary_max_hz)
    primary_max_hz = max(frequencies);
end
scored = score_against_null_core(values, frequencies, two_sided);
primary = frequencies <= primary_max_hz;
scored.primary = score_against_null_core( ...
    values(:, primary), frequencies(primary), two_sided);
end

function scored = score_against_null_core(values, frequencies, two_sided)
scored = empty_scored();
if isempty(values) || size(values, 1) < 2 || all(isnan(values(:)))
    return
end
observed = values(1, :);
null = values(2:end, :);
null_mean = mean(null, 1, 'omitnan');
null_sd = std(null, 0, 1, 'omitnan');
valid = isfinite(observed) & isfinite(null_mean) & ...
    isfinite(null_sd) & null_sd > 0;
if ~any(valid)
    return
end

observed_z = (observed - null_mean) ./ null_sd;
null_z = (null - null_mean) ./ null_sd;
p_point = nan(size(observed));
p_fwer = nan(size(observed));
if two_sided
    null_extreme = max(abs(null_z(:, valid)), [], 2, 'omitnan');
    for frequency_n = find(valid)
        p_point(frequency_n) = (1 + sum( ...
            abs(null_z(:, frequency_n)) >= abs(observed_z(frequency_n)))) / ...
            (size(null, 1) + 1);
        p_fwer(frequency_n) = (1 + sum( ...
            null_extreme >= abs(observed_z(frequency_n)))) / ...
            (size(null, 1) + 1);
    end
    [~, relative_peak] = max(abs(observed(valid)));
    [global_statistic, relative_max_stat] = max(abs(observed_z(valid)));
else
    null_extreme = max(null_z(:, valid), [], 2, 'omitnan');
    for frequency_n = find(valid)
        p_point(frequency_n) = (1 + sum( ...
            null_z(:, frequency_n) >= observed_z(frequency_n))) / ...
            (size(null, 1) + 1);
        p_fwer(frequency_n) = (1 + sum( ...
            null_extreme >= observed_z(frequency_n))) / ...
            (size(null, 1) + 1);
    end
    [~, relative_peak] = max(observed(valid));
    [global_statistic, relative_max_stat] = max(observed_z(valid));
end
valid_index = find(valid);
peak_index = valid_index(relative_peak);
max_stat_index = valid_index(relative_max_stat);
global_p = (1 + sum(null_extreme >= global_statistic)) / ...
    (size(null, 1) + 1);

scored.valid = true;
scored.two_sided = two_sided;
scored.observed_raw = observed;
scored.null_mean = null_mean;
scored.null_sd = null_sd;
scored.observed_z = observed_z;
scored.null_z = single(null_z);
scored.p_point = p_point;
scored.p_fwer = p_fwer;
scored.peak_frequency_hz = frequencies(peak_index);
scored.peak_z = observed_z(peak_index);
scored.peak_p_fwer = p_fwer(peak_index);
scored.max_stat_frequency_hz = frequencies(max_stat_index);
scored.global_p_fwer = global_p;
end

function hierarchy = build_hierarchy(units, options)
effect_names = {'whole', 'play', 'nonplay', 'contrast', 'harmonic'};
animal_ids = unique({units.animal_id}, 'stable');
hierarchy = struct();
hierarchy.frequency_hz = units(1).frequency_hz;
hierarchy.primary_inference_max_hz = options.primary_inference_max_hz;
hierarchy.primary_frequency = hierarchy.frequency_hz <= ...
    hierarchy.primary_inference_max_hz;
hierarchy.phase_theta_overlap_start_hz = ...
    phase_theta_overlap_start(options);
hierarchy.phase_theta_overlap = hierarchy.frequency_hz >= ...
    hierarchy.phase_theta_overlap_start_hz;
hierarchy.effect_names = effect_names;
hierarchy.intermediate = struct();
hierarchy.animals = repmat(struct('animal_id', ''), numel(animal_ids), 1);
hierarchy.group = struct();

rng(options.random_seed + 100000, 'twister');
for effect_n = 1:numel(effect_names)
    effect_name = effect_names{effect_n};
    hierarchy.intermediate.(effect_name) = combine_unit_hierarchy( ...
        units, effect_name, hierarchy.frequency_hz, ...
        options.n_group_permutations, options.primary_inference_max_hz, ...
        strcmp(effect_name, 'contrast'));
end

for animal_n = 1:numel(animal_ids)
    hierarchy.animals(animal_n).animal_id = animal_ids{animal_n};
    index = ismember({units.animal_id}, animal_ids{animal_n});
    for effect_n = 1:numel(effect_names)
        effect_name = effect_names{effect_n};
        hierarchy.animals(animal_n).(effect_name) = combine_unit_hierarchy( ...
            units(index), effect_name, hierarchy.frequency_hz, ...
            options.n_group_permutations, ...
            options.primary_inference_max_hz, ...
            strcmp(effect_name, 'contrast'));
    end
end

for effect_n = 1:numel(effect_names)
    effect_name = effect_names{effect_n};
    animal_effects = cell(0, 1);
    for animal_n = 1:numel(hierarchy.animals)
        effect = hierarchy.animals(animal_n).(effect_name);
        if effect.valid
            animal_effects{end + 1, 1} = effect; %#ok<AGROW>
        end
    end
    if numel(animal_effects) >= options.minimum_animals_group
        hierarchy.group.(effect_name) = combine_animal_sign_flips( ...
            animal_effects, hierarchy.frequency_hz, ...
            options.n_group_permutations, ...
            options.primary_inference_max_hz, ...
            strcmp(effect_name, 'contrast'));
    else
        hierarchy.group.(effect_name) = empty_scored();
        hierarchy.group.(effect_name).primary = empty_scored();
    end
    hierarchy.group.(effect_name).n_animals = numel(animal_effects);
end
hierarchy.interpretation = classify_hierarchy(hierarchy);
end

function interpretation = classify_hierarchy(hierarchy)
whole_full = hierarchy.group.whole;
whole = whole_full.primary;
harmonic = hierarchy.group.harmonic.primary;
contrast = hierarchy.group.contrast.primary;
if ~whole.valid || whole.global_p_fwer >= 0.05
    objective_one = 'not_confirmed';
elseif harmonic.valid && harmonic.global_p_fwer < 0.05
    objective_one = 'coupling_detected_but_waveform_harmonic_ambiguous';
else
    objective_one = 'coupling_detected_harmonic_screen_negative_not_exclusion';
end
if whole_full.valid && whole_full.global_p_fwer < 0.05 && ...
        (whole_full.max_stat_frequency_hz >= ...
        hierarchy.phase_theta_overlap_start_hz || ...
        whole_full.peak_frequency_hz >= ...
        hierarchy.phase_theta_overlap_start_hz)
    full_scan_note = 'full_scan_significant_and_theta_overlap_driven';
elseif whole_full.valid && whole_full.global_p_fwer < 0.05
    full_scan_note = 'full_scan_significant_below_overlap';
else
    full_scan_note = 'full_scan_not_significant';
end
if ~strcmp(objective_one, 'not_confirmed') && ...
        contrast.valid && contrast.global_p_fwer < 0.05
    objective_two = 'play_nonplay_difference_supported';
elseif ~strcmp(objective_one, 'not_confirmed') && contrast.valid
    objective_two = 'play_nonplay_difference_not_supported';
elseif strcmp(objective_one, 'not_confirmed')
    objective_two = 'not_interpreted_because_objective_one_failed';
else
    objective_two = 'insufficient_condition_data';
end
interpretation = struct( ...
    'objective_one', objective_one, ...
    'objective_two', objective_two, ...
    'full_scan_note', full_scan_note, ...
    'rule', ['Objective 1 requires equal-animal whole-session max-stat ' ...
    'p<0.05 in the primary 0.5-4.5 Hz family. Significant harmonic ' ...
    'locking marks the result ambiguous. The 0.5-6 Hz family is reported ' ...
    'separately and cannot replace the primary inference. ' ...
    'A negative harmonic screen does not prove independent oscillators. ' ...
    'Objective 2 requires the direct equal-animal condition interaction, ' ...
    'not significance in only one condition.']);
end

function combined = combine_unit_hierarchy( ...
        units, effect_name, frequencies, n_permutations, ...
        primary_max_hz, two_sided)
combined = empty_scored();
valid = arrayfun(@(x) x.(effect_name).valid, units);
units = units(valid);
combined.n_units = numel(units);
if isempty(units)
    return
end
n_frequency = numel(frequencies);
observed = nan(numel(units), n_frequency);
for unit_n = 1:numel(units)
    observed(unit_n, :) = units(unit_n).(effect_name).observed_z;
end

session_key = arrayfun(@(x) sprintf('%s|%d', ...
    x.recording, x.session_index), units, 'UniformOutput', false);
[unique_session, ~, session_group] = unique(session_key, 'stable');
null_draw = zeros(n_permutations, n_frequency);
total_weight = 0;
for session_n = 1:numel(unique_session)
    index = find(session_group == session_n);
    n_session_units = numel(index);
    n_null = min(arrayfun(@(x) size( ...
        x.(effect_name).null_z, 1), units(index)));
    session_null = zeros(n_null, n_frequency);
    for unit_n = index(:)'
        session_null = session_null + ...
            double(units(unit_n).(effect_name).null_z(1:n_null, :));
    end
    session_null = session_null / n_session_units;
    selected = randi(n_null, n_permutations, 1);
    null_draw = null_draw + ...
        n_session_units * session_null(selected, :);
    total_weight = total_weight + n_session_units;
end
observed_mean = mean(observed, 1, 'omitnan');
null_draw = null_draw / total_weight;
combined = score_against_null( ...
    [observed_mean; null_draw], frequencies, two_sided, primary_max_hz);
combined.n_units = numel(units);
combined.unit_observed_z = observed;
end

function combined = combine_animal_sign_flips( ...
        effects, frequencies, maximum_permutations, ...
        primary_max_hz, two_sided)
combined = empty_scored();
n_animals = numel(effects);
combined.n_animals = n_animals;
if n_animals < 2
    return
end
n_frequency = numel(frequencies);
observed = nan(n_animals, n_frequency);
for animal_n = 1:n_animals
    observed(animal_n, :) = effects{animal_n}.observed_z;
end

if n_animals <= 16
    n_assignment = 2 ^ n_animals;
    codes = uint64((0:n_assignment - 1)');
    signs = -ones(n_assignment, n_animals);
    for animal_n = 1:n_animals
        signs(:, animal_n) = ...
            2 * double(bitget(codes, animal_n)) - 1;
    end
    observed_row = n_assignment;
else
    signs = 2 * (rand(maximum_permutations, n_animals) >= 0.5) - 1;
    signs = [ones(1, n_animals); signs];
    observed_row = 1;
end

statistics = signed_animal_mean(signs, observed);
combined = score_assignment_permutations( ...
    statistics, observed_row, frequencies, two_sided, primary_max_hz);
combined.n_animals = n_animals;
combined.n_units = n_animals;
combined.unit_observed_z = observed;
end

function statistic = signed_animal_mean(signs, observed)
n_draw = size(signs, 1);
n_frequency = size(observed, 2);
statistic = nan(n_draw, n_frequency);
for frequency_n = 1:n_frequency
    valid = isfinite(observed(:, frequency_n));
    n = sum(valid);
    if n < 2
        continue
    end
    signed = signs(:, valid) .* observed(valid, frequency_n)';
    statistic(:, frequency_n) = mean(signed, 2);
end
end

function scored = score_assignment_permutations( ...
        statistics, observed_row, frequencies, two_sided, primary_max_hz)
if nargin < 5 || isempty(primary_max_hz)
    primary_max_hz = max(frequencies);
end
scored = score_assignment_permutations_core( ...
    statistics, observed_row, frequencies, two_sided);
primary = frequencies <= primary_max_hz;
scored.primary = score_assignment_permutations_core( ...
    statistics(:, primary), observed_row, frequencies(primary), two_sided);
end

function scored = score_assignment_permutations_core( ...
        statistics, observed_row, frequencies, two_sided)
scored = empty_scored();
assignment_mean = mean(statistics, 1, 'omitnan');
assignment_sd = std(statistics, 0, 1, 'omitnan');
z = (statistics - assignment_mean) ./ assignment_sd;
observed_raw = statistics(observed_row, :);
observed_z = z(observed_row, :);
valid = isfinite(observed_raw) & isfinite(observed_z) & ...
    assignment_sd > 0;
if ~any(valid)
    return
end

if two_sided
    extreme = max(abs(z(:, valid)), [], 2, 'omitnan');
    [~, relative_peak] = max(abs(observed_raw(valid)));
    [global_statistic, relative_max_stat] = max(abs(observed_z(valid)));
else
    extreme = max(z(:, valid), [], 2, 'omitnan');
    [~, relative_peak] = max(observed_raw(valid));
    [global_statistic, relative_max_stat] = max(observed_z(valid));
end

p_point = nan(size(observed_raw));
p_fwer = nan(size(observed_raw));
for frequency_n = find(valid)
    if two_sided
        p_point(frequency_n) = mean( ...
            abs(z(:, frequency_n)) >= abs(observed_z(frequency_n)));
        threshold = abs(observed_z(frequency_n));
    else
        p_point(frequency_n) = mean( ...
            z(:, frequency_n) >= observed_z(frequency_n));
        threshold = observed_z(frequency_n);
    end
    p_fwer(frequency_n) = mean(extreme >= threshold);
end
valid_index = find(valid);
peak_index = valid_index(relative_peak);
max_stat_index = valid_index(relative_max_stat);

scored.valid = true;
scored.two_sided = two_sided;
scored.observed_raw = observed_raw;
scored.null_mean = assignment_mean;
scored.null_sd = assignment_sd;
scored.observed_z = observed_z;
scored.null_z = single(z((1:size(z, 1)) ~= observed_row, :));
scored.p_point = p_point;
scored.p_fwer = p_fwer;
scored.peak_frequency_hz = frequencies(peak_index);
scored.peak_z = observed_z(peak_index);
scored.peak_p_fwer = p_fwer(peak_index);
scored.max_stat_frequency_hz = frequencies(max_stat_index);
scored.global_p_fwer = mean(extreme >= global_statistic);
end

function [unit_table, animal_table, group_table] = result_tables( ...
        units, hierarchy)
n_units = numel(units);
animal = strings(n_units, 1);
recording = strings(n_units, 1);
channel = nan(n_units, 1);
session = nan(n_units, 1);
whole_seconds = nan(n_units, 1);
play_seconds = nan(n_units, 1);
nonplay_seconds = nan(n_units, 1);
whole_peak_hz = nan(n_units, 1);
whole_peak_z = nan(n_units, 1);
whole_peak_p_fwer = nan(n_units, 1);
whole_maxstat_hz = nan(n_units, 1);
whole_p_fwer = nan(n_units, 1);
whole_primary_maxstat_hz = nan(n_units, 1);
whole_primary_p_fwer = nan(n_units, 1);
play_p_fwer = nan(n_units, 1);
play_primary_p_fwer = nan(n_units, 1);
nonplay_p_fwer = nan(n_units, 1);
nonplay_primary_p_fwer = nan(n_units, 1);
contrast_peak_hz = nan(n_units, 1);
contrast_peak_z = nan(n_units, 1);
contrast_peak_p_fwer = nan(n_units, 1);
contrast_maxstat_hz = nan(n_units, 1);
contrast_p_fwer = nan(n_units, 1);
contrast_primary_maxstat_hz = nan(n_units, 1);
contrast_primary_p_fwer = nan(n_units, 1);
harmonic_p_fwer = nan(n_units, 1);
harmonic_primary_p_fwer = nan(n_units, 1);
for unit_n = 1:n_units
    animal(unit_n) = string(units(unit_n).animal_id);
    recording(unit_n) = string(units(unit_n).recording);
    channel(unit_n) = units(unit_n).lpag_lfp_row;
    session(unit_n) = units(unit_n).session_index;
    whole_seconds(unit_n) = units(unit_n).analysis_duration_seconds;
    play_seconds(unit_n) = units(unit_n).play_seconds;
    nonplay_seconds(unit_n) = units(unit_n).nonplay_seconds;
    [whole_peak_hz(unit_n), whole_peak_z(unit_n), whole_p_fwer(unit_n)] = ...
        effect_peak(units(unit_n).whole);
    whole_peak_p_fwer(unit_n) = ...
        effect_value(units(unit_n).whole, 'peak_p_fwer');
    whole_maxstat_hz(unit_n) = ...
        effect_value(units(unit_n).whole, 'max_stat_frequency_hz');
    whole_primary_maxstat_hz(unit_n) = effect_primary_value( ...
        units(unit_n).whole, 'max_stat_frequency_hz');
    whole_primary_p_fwer(unit_n) = effect_primary_value( ...
        units(unit_n).whole, 'global_p_fwer');
    [~, ~, play_p_fwer(unit_n)] = effect_peak(units(unit_n).play);
    play_primary_p_fwer(unit_n) = effect_primary_value( ...
        units(unit_n).play, 'global_p_fwer');
    [~, ~, nonplay_p_fwer(unit_n)] = effect_peak(units(unit_n).nonplay);
    nonplay_primary_p_fwer(unit_n) = effect_primary_value( ...
        units(unit_n).nonplay, 'global_p_fwer');
    [contrast_peak_hz(unit_n), contrast_peak_z(unit_n), ...
        contrast_p_fwer(unit_n)] = effect_peak(units(unit_n).contrast);
    contrast_peak_p_fwer(unit_n) = ...
        effect_value(units(unit_n).contrast, 'peak_p_fwer');
    contrast_maxstat_hz(unit_n) = ...
        effect_value(units(unit_n).contrast, 'max_stat_frequency_hz');
    contrast_primary_maxstat_hz(unit_n) = effect_primary_value( ...
        units(unit_n).contrast, 'max_stat_frequency_hz');
    contrast_primary_p_fwer(unit_n) = effect_primary_value( ...
        units(unit_n).contrast, 'global_p_fwer');
    [~, ~, harmonic_p_fwer(unit_n)] = effect_peak(units(unit_n).harmonic);
    harmonic_primary_p_fwer(unit_n) = effect_primary_value( ...
        units(unit_n).harmonic, 'global_p_fwer');
end
unit_table = table(animal, recording, channel, session, whole_seconds, ...
    play_seconds, nonplay_seconds, whole_peak_hz, whole_peak_z, ...
    whole_peak_p_fwer, whole_maxstat_hz, whole_p_fwer, ...
    whole_primary_maxstat_hz, whole_primary_p_fwer, ...
    play_p_fwer, play_primary_p_fwer, ...
    nonplay_p_fwer, nonplay_primary_p_fwer, contrast_peak_hz, ...
    contrast_peak_z, contrast_peak_p_fwer, contrast_maxstat_hz, ...
    contrast_p_fwer, contrast_primary_maxstat_hz, ...
    contrast_primary_p_fwer, harmonic_p_fwer, ...
    harmonic_primary_p_fwer);

n_animals = numel(hierarchy.animals);
animal = strings(n_animals, 1);
n_units = nan(n_animals, 1);
whole_peak_hz = nan(n_animals, 1);
whole_peak_z = nan(n_animals, 1);
whole_peak_p_fwer = nan(n_animals, 1);
whole_maxstat_hz = nan(n_animals, 1);
whole_p_fwer = nan(n_animals, 1);
whole_primary_maxstat_hz = nan(n_animals, 1);
whole_primary_p_fwer = nan(n_animals, 1);
play_p_fwer = nan(n_animals, 1);
play_primary_p_fwer = nan(n_animals, 1);
nonplay_p_fwer = nan(n_animals, 1);
nonplay_primary_p_fwer = nan(n_animals, 1);
contrast_peak_hz = nan(n_animals, 1);
contrast_peak_z = nan(n_animals, 1);
contrast_peak_p_fwer = nan(n_animals, 1);
contrast_maxstat_hz = nan(n_animals, 1);
contrast_p_fwer = nan(n_animals, 1);
contrast_primary_maxstat_hz = nan(n_animals, 1);
contrast_primary_p_fwer = nan(n_animals, 1);
harmonic_p_fwer = nan(n_animals, 1);
harmonic_primary_p_fwer = nan(n_animals, 1);
for animal_n = 1:n_animals
    animal(animal_n) = string(hierarchy.animals(animal_n).animal_id);
    n_units(animal_n) = hierarchy.animals(animal_n).whole.n_units;
    [whole_peak_hz(animal_n), whole_peak_z(animal_n), ...
        whole_p_fwer(animal_n)] = effect_peak( ...
        hierarchy.animals(animal_n).whole);
    whole_peak_p_fwer(animal_n) = effect_value( ...
        hierarchy.animals(animal_n).whole, 'peak_p_fwer');
    whole_maxstat_hz(animal_n) = effect_value( ...
        hierarchy.animals(animal_n).whole, 'max_stat_frequency_hz');
    whole_primary_maxstat_hz(animal_n) = effect_primary_value( ...
        hierarchy.animals(animal_n).whole, 'max_stat_frequency_hz');
    whole_primary_p_fwer(animal_n) = effect_primary_value( ...
        hierarchy.animals(animal_n).whole, 'global_p_fwer');
    [~, ~, play_p_fwer(animal_n)] = effect_peak( ...
        hierarchy.animals(animal_n).play);
    play_primary_p_fwer(animal_n) = effect_primary_value( ...
        hierarchy.animals(animal_n).play, 'global_p_fwer');
    [~, ~, nonplay_p_fwer(animal_n)] = effect_peak( ...
        hierarchy.animals(animal_n).nonplay);
    nonplay_primary_p_fwer(animal_n) = effect_primary_value( ...
        hierarchy.animals(animal_n).nonplay, 'global_p_fwer');
    [contrast_peak_hz(animal_n), contrast_peak_z(animal_n), ...
        contrast_p_fwer(animal_n)] = effect_peak( ...
        hierarchy.animals(animal_n).contrast);
    contrast_peak_p_fwer(animal_n) = effect_value( ...
        hierarchy.animals(animal_n).contrast, 'peak_p_fwer');
    contrast_maxstat_hz(animal_n) = effect_value( ...
        hierarchy.animals(animal_n).contrast, 'max_stat_frequency_hz');
    contrast_primary_maxstat_hz(animal_n) = effect_primary_value( ...
        hierarchy.animals(animal_n).contrast, 'max_stat_frequency_hz');
    contrast_primary_p_fwer(animal_n) = effect_primary_value( ...
        hierarchy.animals(animal_n).contrast, 'global_p_fwer');
    [~, ~, harmonic_p_fwer(animal_n)] = effect_peak( ...
        hierarchy.animals(animal_n).harmonic);
    harmonic_primary_p_fwer(animal_n) = effect_primary_value( ...
        hierarchy.animals(animal_n).harmonic, 'global_p_fwer');
end
animal_table = table(animal, n_units, whole_peak_hz, whole_peak_z, ...
    whole_peak_p_fwer, whole_maxstat_hz, whole_p_fwer, ...
    whole_primary_maxstat_hz, whole_primary_p_fwer, ...
    play_p_fwer, play_primary_p_fwer, ...
    nonplay_p_fwer, nonplay_primary_p_fwer, contrast_peak_hz, ...
    contrast_peak_z, contrast_peak_p_fwer, contrast_maxstat_hz, ...
    contrast_p_fwer, contrast_primary_maxstat_hz, ...
    contrast_primary_p_fwer, harmonic_p_fwer, ...
    harmonic_primary_p_fwer);

levels = {'Intermediate pooled units'; 'Equal-animal group'};
effect_names = hierarchy.effect_names(:);
n_rows = numel(levels) * numel(effect_names);
level = strings(n_rows, 1);
effect = strings(n_rows, 1);
n_contributors = nan(n_rows, 1);
peak_hz = nan(n_rows, 1);
peak_z = nan(n_rows, 1);
peak_p_fwer = nan(n_rows, 1);
maxstat_hz = nan(n_rows, 1);
p_fwer = nan(n_rows, 1);
primary_peak_hz = nan(n_rows, 1);
primary_peak_z = nan(n_rows, 1);
primary_peak_p_fwer = nan(n_rows, 1);
primary_maxstat_hz = nan(n_rows, 1);
primary_p_fwer = nan(n_rows, 1);
row = 0;
for level_n = 1:numel(levels)
    for effect_n = 1:numel(effect_names)
        row = row + 1;
        level(row) = string(levels{level_n});
        effect(row) = string(effect_names{effect_n});
        if level_n == 1
            result = hierarchy.intermediate.(effect_names{effect_n});
            n_contributors(row) = result.n_units;
        else
            result = hierarchy.group.(effect_names{effect_n});
            n_contributors(row) = result.n_animals;
        end
        [peak_hz(row), peak_z(row), p_fwer(row)] = effect_peak(result);
        peak_p_fwer(row) = effect_value(result, 'peak_p_fwer');
        maxstat_hz(row) = effect_value(result, 'max_stat_frequency_hz');
        [primary_peak_hz(row), primary_peak_z(row), ...
            primary_p_fwer(row)] = effect_peak(result.primary);
        primary_peak_p_fwer(row) = effect_value( ...
            result.primary, 'peak_p_fwer');
        primary_maxstat_hz(row) = effect_value( ...
            result.primary, 'max_stat_frequency_hz');
    end
end
group_table = table(level, effect, n_contributors, peak_hz, peak_z, ...
    peak_p_fwer, maxstat_hz, p_fwer, primary_peak_hz, ...
    primary_peak_z, primary_peak_p_fwer, primary_maxstat_hz, ...
    primary_p_fwer);
end

function [peak, z, p] = effect_peak(effect)
peak = NaN;
z = NaN;
p = NaN;
if isfield(effect, 'valid') && effect.valid
    peak = effect.peak_frequency_hz;
    z = effect.peak_z;
    p = effect.global_p_fwer;
end
end

function value = effect_value(effect, field_name)
value = NaN;
if isfield(effect, 'valid') && effect.valid && isfield(effect, field_name)
    value = effect.(field_name);
end
end

function value = effect_primary_value(effect, field_name)
value = NaN;
if isfield(effect, 'primary') && ...
        isfield(effect.primary, 'valid') && effect.primary.valid && ...
        isfield(effect.primary, field_name)
    value = effect.primary.(field_name);
end
end

function print_summary(units, hierarchy)
fprintf('\n=== Maximum-statistic corrected summary ===\n');
effect_names = {'whole', 'play', 'nonplay', 'contrast', 'harmonic'};
for effect_n = 1:numel(effect_names)
    name = effect_names{effect_n};
    valid = arrayfun(@(x) x.(name).valid, units);
    significant = false(size(valid));
    significant_primary = false(size(valid));
    for unit_n = find(valid)
        significant(unit_n) = units(unit_n).(name).global_p_fwer < 0.05;
        if isfield(units(unit_n).(name).primary, 'valid') && ...
                units(unit_n).(name).primary.valid
            significant_primary(unit_n) = ...
                units(unit_n).(name).primary.global_p_fwer < 0.05;
        end
    end
    intermediate = hierarchy.intermediate.(name);
    group = hierarchy.group.(name);
    fprintf(['%s: unit primary %d/%d, full %d/%d significant | primary ' ...
        'intermediate p=%s | primary equal-animal p=%s | ' ...
        'full equal-animal p=%s (N=%d animals)\n'], ...
        name, sum(significant_primary), sum(valid), ...
        sum(significant), sum(valid), ...
        format_p(intermediate.primary), format_p(group.primary), ...
        format_p(group), group.n_animals);
end
fprintf(['Interpret PAC as independent delta-theta coupling only if the PAC ' ...
    'test survives correction and the harmonic control does not provide ' ...
    'a plausible waveform explanation. For play/non-play specificity, ' ...
    'require the contrast test; separate significance is insufficient.\n']);
fprintf(['Phase frequencies >= %.2f Hz overlap the theta band within two ' ...
    'wavelet spectral SD and are exploratory/ambiguous.\n'], ...
    hierarchy.phase_theta_overlap_start_hz);
fprintf('Objective 1 classification: %s\n', ...
    hierarchy.interpretation.objective_one);
fprintf('Objective 2 classification: %s\n', ...
    hierarchy.interpretation.objective_two);
fprintf('Full-scan note: %s\n', ...
    hierarchy.interpretation.full_scan_note);
end

function text = format_p(effect)
if ~isfield(effect, 'valid') || ~effect.valid
    text = 'NA';
else
    text = sprintf('%.4g', effect.global_p_fwer);
end
end

function make_summary_figures(units, hierarchy, output_dir)
frequency = hierarchy.frequency_hz;
figure_one = figure('Color', 'w', 'Position', [50 50 1200 850]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
unit_z = effect_matrix(units, 'whole');
imagesc(frequency, 1:size(unit_z, 1), unit_z);
axis xy
colorbar
xlabel('Candidate phase frequency (Hz)')
ylabel('Session x LPAG channel')
title('Whole partner sessions: PAC z-score versus circular-shift null')
nexttile
animal_z = animal_effect_matrix( ...
    hierarchy.animals, 'whole', numel(frequency));
imagesc(frequency, 1:size(animal_z, 1), animal_z);
axis xy
colorbar
xlabel('Candidate phase frequency (Hz)')
ylabel('Animal')
title('Animal-level PAC')
nexttile
plot_group_effect(frequency, animal_z, hierarchy.group.whole, ...
    'Equal-animal whole-session PAC', ...
    hierarchy.phase_theta_overlap_start_hz, ...
    hierarchy.primary_inference_max_hz);
save_figure_pair(figure_one, output_dir, 'theta_pac_whole_session');

figure_two = figure('Color', 'w', 'Position', [70 70 1200 850]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
effect_names = {'play', 'nonplay', 'contrast', 'harmonic'};
titles = {'Play PAC', 'Non-play PAC', ...
    'Play minus non-play interaction', 'Slow/theta harmonic control'};
for effect_n = 1:numel(effect_names)
    nexttile
    matrix = animal_effect_matrix( ...
        hierarchy.animals, effect_names{effect_n}, numel(frequency));
    plot_group_effect(frequency, matrix, ...
        hierarchy.group.(effect_names{effect_n}), titles{effect_n}, ...
        hierarchy.phase_theta_overlap_start_hz, ...
        hierarchy.primary_inference_max_hz);
end
save_figure_pair(figure_two, output_dir, 'theta_pac_conditions');
end

function matrix = effect_matrix(units, effect_name)
matrix = nan(numel(units), numel(units(1).frequency_hz));
for unit_n = 1:numel(units)
    effect = units(unit_n).(effect_name);
    if effect.valid
        matrix(unit_n, :) = effect.observed_z;
    end
end
end

function matrix = animal_effect_matrix(animals, effect_name, n_frequency)
if isempty(animals)
    matrix = nan(0, n_frequency);
    return
end
matrix = nan(numel(animals), n_frequency);
for animal_n = 1:numel(animals)
    effect = animals(animal_n).(effect_name);
    if effect.valid
        matrix(animal_n, :) = effect.observed_z;
    end
end
end

function plot_group_effect(frequency, animal_matrix, group, plot_title, ...
        overlap_start_hz, primary_max_hz)
hold on
plot(frequency, animal_matrix', 'Color', [0.82 0.82 0.82]);
mean_curve = mean(animal_matrix, 1, 'omitnan');
n = sum(isfinite(animal_matrix), 1);
sem = std(animal_matrix, 0, 1, 'omitnan') ./ sqrt(n);
ci = 1.96 * sem;
fill([frequency fliplr(frequency)], ...
    [mean_curve - ci fliplr(mean_curve + ci)], ...
    [0.25 0.45 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(frequency, mean_curve, 'k', 'LineWidth', 2);
yline(0, 'k:');
xline(primary_max_hz, 'b:');
xline(overlap_start_hz, 'r--');
if group.valid
    full_significant = group.p_fwer < 0.05;
    plot(frequency(full_significant), mean_curve(full_significant), ...
        'o', 'Color', [0.9 0.45 0], 'MarkerSize', 5);
    primary_frequency = frequency <= primary_max_hz;
    if isfield(group.primary, 'valid') && group.primary.valid
        primary_significant = group.primary.p_fwer < 0.05;
        primary_significant_full = false(size(frequency));
        primary_significant_full(primary_frequency) = primary_significant;
        plot(frequency(primary_significant_full), ...
            mean_curve(primary_significant_full), ...
            'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
        primary_p = group.primary.global_p_fwer;
    else
        primary_p = NaN;
    end
    title(sprintf('%s; primary p=%.4g, full p=%.4g', ...
        plot_title, primary_p, group.global_p_fwer));
else
    title([plot_title, '; insufficient animals']);
end
xlabel('Candidate phase frequency (Hz)')
ylabel('PAC z-score')
box off
end

function save_figure_pair(figure_handle, output_dir, base_name)
savefig(figure_handle, fullfile(output_dir, [base_name, '.fig']));
print(figure_handle, fullfile(output_dir, [base_name, '.png']), ...
    '-dpng', '-r200');
end

function intervals = behavior_intervals(behavior, types, merge_gap)
labels = as_labels(behavior.Type2);
keep = ismember(labels, types);
intervals = behavior{keep, {'Start', 'End'}};
intervals = intervals(all(isfinite(intervals), 2) & ...
    intervals(:, 2) > intervals(:, 1), :);
intervals = merge_intervals(intervals, merge_gap);
end

function merged = merge_intervals(intervals, gap)
if isempty(intervals)
    merged = zeros(0, 2);
    return
end
intervals = sortrows(intervals, 1);
merged = intervals(1, :);
for interval_n = 2:size(intervals, 1)
    if intervals(interval_n, 1) <= merged(end, 2) + gap
        merged(end, 2) = max(merged(end, 2), intervals(interval_n, 2));
    else
        merged(end + 1, :) = intervals(interval_n, :); %#ok<AGROW>
    end
end
end

function mask = interval_mask(time, intervals)
mask = false(size(time));
for interval_n = 1:size(intervals, 1)
    mask = mask | (time >= intervals(interval_n, 1) & ...
        time <= intervals(interval_n, 2));
end
end

function [behavior, partner_intervals] = load_behavior( ...
        recording_name, options)
behavior_file = fullfile(options.data_root, ...
    'Behavior backups', [recording_name, '.txt']);
sync_file = fullfile(options.data_root, 'Synch data', ...
    recording_name, 'synch_model_video2NPX.mat');
behavior = readtable(behavior_file);
behavior(:, 2) = [];
behavior.Properties.VariableNames = ...
    {'Animal', 'Start', 'End', 'Length', 'Type'};
behavior.Type2 = behavior.Type;
type = as_labels(behavior.Type2);
type(ismember(type, {'Pounce_A', 'Pounce_B'})) = {'Pounce'};
type(ismember(type, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
type(cellfun(@isempty, type)) = {'Other'};
behavior.Type2 = type;
behavior(ismember(as_labels(behavior.Animal), 'Reversal'), :) = [];

sync_data = load(sync_file);
if isfield(sync_data, 'synch_model_video2NPX')
    sync_model = sync_data.synch_model_video2NPX;
else
    fields = fieldnames(sync_data);
    if numel(fields) ~= 1
        error('Cannot identify the synchronization model in %s.', sync_file);
    end
    sync_model = sync_data.(fields{1});
end
behavior.Start = predict(sync_model, behavior.Start);
behavior.End = predict(sync_model, behavior.End);

animal = as_labels(behavior.Animal);
type = as_labels(behavior.Type);
session = ismember(animal, 'Session_structure') & ...
    ismember(type, 'Partners session');
if ~any(session)
    session = ismember(animal, 'Session_structure') & ...
        ~ismember(type, 'Tickling');
end
partner_intervals = behavior{session, {'Start', 'End'}};
partner_intervals = sortrows(partner_intervals, 1);
end

function labels = as_labels(values)
if iscell(values)
    labels = cellfun(@char, values, 'UniformOutput', false);
elseif isstring(values)
    labels = cellstr(values);
elseif iscategorical(values)
    labels = cellstr(values);
else
    labels = cellstr(string(values));
end
labels = labels(:);
end

function seed = deterministic_seed(base_seed, recording_name, session_index)
codes = double(char(recording_name));
weights = 1:numel(codes);
name_hash = mod(sum(codes .* weights), 1000000);
seed = mod(double(base_seed) + name_hash + ...
    1009 * double(session_index), 2 ^ 31 - 2) + 1;
end

function [recording_name, animal_id, repeated_animal] = ...
        recording_names(recording_dir)
[~, recording_name] = fileparts(recording_dir);
parts = strsplit(strtrim(recording_name));
parts = parts(~cellfun(@isempty, parts));
if numel(parts) < 3
    error('Unexpected Neuropixels folder name: %s', recording_name);
end
animal_batch = parts{1};
repeated_animal = parts{3};
if strcmp(repeated_animal, 'Single2')
    animal_id = ['Batch', animal_batch(2), repeated_animal];
else
    if numel(animal_batch) < 4
        error('Cannot derive animal identifier from %s.', recording_name);
    end
    animal_id = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
end

function [xcoords, ycoords, chan_map] = load_channel_map(recording_dir)
map_file = fullfile(recording_dir, 'chann_map_PAG.mat');
loaded = load(map_file, 'xcoords', 'ycoords', 'chanMap');
xcoords = loaded.xcoords;
ycoords = loaded.ycoords;
chan_map = loaded.chanMap;
end

function npx_type = neuropixels_type(xcoords)
hard_coded_x = [8 40; 258 290; 508 540; 758 790];
if any(ismember(xcoords(:), hard_coded_x(:)))
    npx_type = 2;
else
    npx_type = 1;
end
end

function [channel_lfp, n_samples] = read_selected_lfp( ...
        recording_dir, rows, n_lfp_channels)
rows = rows(:);
mat_file = fullfile(recording_dir, 'LFP_PAG.mat');
dat_file = fullfile(recording_dir, 'LFP_PAG.dat');
if isfile(mat_file)
    if is_mat_v73(mat_file)
        mapped = matfile(mat_file);
        info = whos(mapped, 'LFP');
        n_samples = info.size(2);
        channel_lfp = double(mapped.LFP(rows, 1:n_samples));
    else
        warning('FrequencyPAC:NonV73LFP', ...
            ['%s is not v7.3 and cannot be read by row. MATLAB must load ' ...
            'the complete LFP variable; convert this file to v7.3 if memory ' ...
            'is insufficient.'], mat_file);
        loaded = load(mat_file, 'LFP');
        channel_lfp = double(loaded.LFP(rows, :));
        n_samples = size(channel_lfp, 2);
    end
    return
end
if ~isfile(dat_file)
    error('Neither LFP_PAG.mat nor LFP_PAG.dat exists in %s.', recording_dir);
end
file_info = dir(dat_file);
n_samples = file_info.bytes / (n_lfp_channels * 2);
if n_samples ~= round(n_samples)
    error('LFP_PAG.dat size is incompatible with %d channels.', n_lfp_channels);
end
channel_lfp = zeros(numel(rows), n_samples);
file_id = fopen(dat_file, 'r', 'ieee-le');
if file_id < 0
    error('Could not open %s.', dat_file);
end
cleaner = onCleanup(@() fclose(file_id)); %#ok<NASGU>
for row_n = 1:numel(rows)
    fseek(file_id, 2 * (rows(row_n) - 1), 'bof');
    [samples, count] = fread(file_id, n_samples, ...
        'int16=>double', 2 * (n_lfp_channels - 1));
    if count ~= n_samples
        error('Read %d of %d samples from LFP row %d.', ...
            count, n_samples, rows(row_n));
    end
    channel_lfp(row_n, :) = samples(:)';
end
end

function answer = is_mat_v73(filename)
file_id = fopen(filename, 'r');
if file_id < 0
    error('Could not open %s to inspect its MAT version.', filename);
end
cleaner = onCleanup(@() fclose(file_id)); %#ok<NASGU>
header = fread(file_id, 128, '*char')';
answer = contains(header, 'MATLAB 7.3 MAT-file');
end

function save_checkpoint(output_file, checkpoint)
temporary_file = [output_file, '.tmp.mat'];
save(temporary_file, '-struct', 'checkpoint', '-v7.3');
[success, message] = movefile(temporary_file, output_file, 'f');
if ~success
    error('Could not replace checkpoint: %s', message);
end
end

function units = empty_pac_units()
units = repmat(empty_unit(), 0, 1);
end

function unit = empty_unit()
unit = struct( ...
    'animal_id', '', ...
    'recording', '', ...
    'repeated_animal', '', ...
    'lpag_channel', NaN, ...
    'lpag_lfp_row', NaN, ...
    'probe_number', NaN, ...
    'session_index', NaN, ...
    'session_bounds_seconds', [NaN NaN], ...
    'frequency_hz', [], ...
    'phase_theta_overlap', [], ...
    'phase_theta_overlap_start_hz', NaN, ...
    'analysis_duration_seconds', NaN, ...
    'play_seconds', NaN, ...
    'nonplay_seconds', NaN, ...
    'overlap_removed_seconds', NaN, ...
    'n_play_blocks', 0, ...
    'n_nonplay_blocks', 0, ...
    'edge_trim_seconds', NaN, ...
    'surrogate_block_shifts_samples', [], ...
    'whole', empty_scored(), ...
    'play', empty_scored(), ...
    'nonplay', empty_scored(), ...
    'contrast', empty_scored(), ...
    'harmonic', empty_scored(), ...
    'preferred_phase_whole', [], ...
    'preferred_phase_play', [], ...
    'preferred_phase_nonplay', [], ...
    'slow_wavelet_power', [], ...
    'notes', '');
end

function scored = empty_scored()
scored = struct( ...
    'valid', false, ...
    'two_sided', false, ...
    'observed_raw', [], ...
    'null_mean', [], ...
    'null_sd', [], ...
    'observed_z', [], ...
    'null_z', single([]), ...
    'p_point', [], ...
    'p_fwer', [], ...
    'peak_frequency_hz', NaN, ...
    'peak_z', NaN, ...
    'peak_p_fwer', NaN, ...
    'max_stat_frequency_hz', NaN, ...
    'global_p_fwer', NaN, ...
    'n_units', 0, ...
    'n_animals', 0, ...
    'unit_observed_z', [], ...
    'primary', struct());
end

function failures = empty_failures()
failures = repmat(struct( ...
    'recording', '', ...
    'channel', NaN, ...
    'session', NaN, ...
    'identifier', '', ...
    'message', '', ...
    'report', '', ...
    'time', ''), 0, 1);
end
