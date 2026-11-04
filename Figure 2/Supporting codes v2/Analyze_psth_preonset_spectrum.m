%% Analyze_psth_preonset_spectrum
% Spectrum of the pre-onset theta power trace.
%
% Run Estimate_psth_all_behaviors_theta_gamma.m first. This script loads the
% same psth_structure_*_all_behaviors_self_other_partner.mat files as
% Analyze_psth_all_behaviors_theta_gamma.m. Grand-average condition plots
% stay in that script.
%
% band_list is theta only. To run the other saved bands as well, set
% band_list = {'delta', 'theta', 'gamma'} in the Paths section.
%
% Run this file from the top, or run the sections in order inside the band
% loop. Each band reloads its own file.
%
% %% Pre-onset periodic spectrum
%     FOOOF periodic component of the electrode-1 trial average, from -5 s
%     to onset. Same figures previously drawn at the end of
%     Analyze_psth_all_behaviors_theta_gamma.m.
% %% Per-animal waveforms
%     One colored line per animal. The black line pools every trial, so an
%     animal with more trials pulls it harder. The gray line gives every
%     animal with at least one trial the same weight. The dashed line is
%     the pooled waveform after the animal with the most trials is removed.
%     If that animal is the only one, there is no dashed line.
%     Behaviors are split across figures: play 1-4, play 5-8, then the six
%     non-play behaviors.
% %% Per-animal spectra
%     FOOOF periodic component of each of those waveforms. The fit and the
%     panels use 0-20 Hz. The Command Window reports the mean periodic
%     component from 2 to 3 Hz. 20 Hz is below the Nyquist frequency of the
%     10 ms power step.
% %% Raw mean versus per-trial z-score
%     The preview mean, and the mean after each trial is z-scored on its own
%     -2 to 0 s window.
%
% Read the Command Window this way. "share" is that animal's fraction of
% electrode-1 trials. "r" is the correlation of its waveform with the
% pooled-trial waveform. "pooled", "own", "animals-equal", and "remaining"
% use the mean FOOOF periodic component from 2 to 3 Hz, the same curves as
% the spectrum panels. "remaining" divides that value after the animal is
% removed by the value with everyone included. Remaining near 1 means the
% peak survives without that animal. Remaining near 0 means that animal
% carries the peak. Remaining above 1 means that animal was out of phase
% with the pool. Remaining is 0 when that animal is the only one with
% trials. "animals-equal" is the gray line. A much smaller animals-equal
% value than pooled means animals with more trials inflated the peak.
% A whole-session z-score is not in these files. The last section is only
% the per-trial baseline z-score. Scratching is left out of the non-play
% spectrum stack, as in the grand-average script.

%% Paths
this_file = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Analyze_psth_preonset_spectrum.m';
repo_root = fileparts(fileparts(fileparts(this_file)));
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');
addpath(fileparts(this_file));

all_behaviors = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD', ...
    'Grooming', 'Rearing', 'Sniffing', 'Scratching', 'PounceI', 'Bite'};
animal_label = {'B1D1', 'B1S3', 'B2S2', 'B3D2', 'B4S2', 'B4D4'};
play_behaviors = 1:8;
other_behaviors = [9 10 11 13 14];
behavior_pages = {1:4, 5:8, 9:14};
behavior_page_names = {'play_1', 'play_2', 'nonplay'};
% Theta only. Uncomment the next line to also plot delta and gamma.
band_list = {'theta'};
% band_list = {'delta', 'theta', 'gamma'};
pre_window = [-5 0];
effect_band = [2 3];
% 20 Hz is the top of the fit and of the figures. The power trace is
% sampled every 10 ms, so the pre-onset FFT can reach 50 Hz. The loop
% prints that Nyquist frequency and will not draw past it.
fooof_fit_max_hz = 20;
plot_max_hz = fooof_fit_max_hz;
fooof_peak_widths = [0.5 12];
fooof_aperiodic_mode = 'fixed';
animal_colors = [0.12 0.47 0.71; 0.89 0.10 0.11; 0.20 0.63 0.17; ...
    0.60 0.31 0.64; 1.00 0.50 0.05; 0.09 0.75 0.81];

for band_n = 1:numel(band_list)
band_name = band_list{band_n};
structure_file = ['psth_structure_' band_name '_all_behaviors_self_other_partner.mat'];
names_file = ['animal_names_' band_name '_all_behaviors_self_other_partner.mat'];

%% Load the same PSTH file as the grand-average script
disp(['loading ' band_name])
load(fullfile(saving_folder, structure_file), 'psth_structure');
load(fullfile(saving_folder, names_file), 'animal_names');
disp('ready')

%% Merge sessions, self and partner
bin_size = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
psth_ranges = psth_structure(1).hist_range;

all_psth_onset_self = cell(size(psth_structure(1).all_behavior_psth_self));
all_session_index_self = cell(1, size(psth_structure(1).all_behavior_psth_self, 2));
all_electrode_index_self = cell(1, size(psth_structure(1).all_behavior_psth_self, 2));
all_psth_onset_other = cell(size(psth_structure(1).all_behavior_psth_other));
all_session_index_other = cell(1, size(psth_structure(1).all_behavior_psth_other, 2));
all_electrode_index_other = cell(1, size(psth_structure(1).all_behavior_psth_other, 2));

for j = 1:numel(psth_structure)
    if ~contains(animal_names{j, 1}, animal_label)
        continue
    end
    animal_num = find(cell2mat(cellfun(@(x) contains(animal_names{j, 1}, x), ...
        animal_label, 'UniformOutput', false)));
    electrode_num = animal_names{j, 2};

    for bn = 1:size(all_psth_onset_self, 2)
        n_obs = size(psth_structure(j).all_behavior_psth_self{1, bn}, 1);
        for tp = 1:3
            if tp < 3 && n_obs > 0
                block = psth_structure(j).all_behavior_psth_self{tp, bn}(:, 1:4000);
            else
                block = psth_structure(j).all_behavior_psth_self{tp, bn};
            end
            all_psth_onset_self{tp, bn} = [all_psth_onset_self{tp, bn}; block];
        end
        all_session_index_self{bn} = [all_session_index_self{bn}; animal_num * ones(n_obs, 1)];
        all_electrode_index_self{bn} = [all_electrode_index_self{bn}; electrode_num * ones(n_obs, 1)];
    end

    for bn = 1:size(all_psth_onset_other, 2)
        n_obs = size(psth_structure(j).all_behavior_psth_other{1, bn}, 1);
        for tp = 1:3
            if tp < 3 && n_obs > 0
                block = psth_structure(j).all_behavior_psth_other{tp, bn}(:, 1:4000);
            else
                block = psth_structure(j).all_behavior_psth_other{tp, bn};
            end
            all_psth_onset_other{tp, bn} = [all_psth_onset_other{tp, bn}; block];
        end
        all_session_index_other{bn} = [all_session_index_other{bn}; animal_num * ones(n_obs, 1)];
        all_electrode_index_other{bn} = [all_electrode_index_other{bn}; electrode_num * ones(n_obs, 1)];
    end
end

%% Pre-onset periodic spectrum
% Electrode 1. Trials are demeaned on -5 to 0 s, averaged, then Fourier
% transformed. FOOOF keeps the periodic part. The fit and the figures use
% 0-20 Hz, or the Nyquist frequency when that is lower. The gray band is 1-5 Hz.
n_col = 4000;
for b = 1:numel(all_behaviors)
    if ~isempty(all_psth_onset_self{1, b})
        n_col = size(all_psth_onset_self{1, b}, 2);
        break
    end
end
bin_offsets = floor(psth_ranges(1) / bin_size) + (0:n_col - 1);
time_psth = bin_offsets * bin_size;
pre_mask = time_psth >= pre_window(1) & time_psth < pre_window(2);
time_pre = time_psth(pre_mask);
n_pre = nnz(pre_mask);
nyquist_hz = 1 / (2 * bin_size);
freq_step_hz = (1 / bin_size) / n_pre;
plot_max_hz = min(fooof_fit_max_hz, nyquist_hz);
fprintf(['Pre-onset FFT: step %.0f ms, %d samples, %.2f Hz resolution, ' ...
    'Nyquist %.0f Hz. Figures show the fit to %.0f Hz.\n'], ...
    bin_size * 1000, n_pre, freq_step_hz, nyquist_hz, plot_max_hz);
freq_fit = (1:floor(n_pre / 2)) * ((1 / bin_size) / n_pre);
freq_fit = freq_fit(freq_fit <= fooof_fit_max_hz);
freq_hz = freq_fit(freq_fit <= plot_max_hz);

spec_self = nan(numel(all_behaviors), numel(freq_fit));
spec_other = nan(numel(all_behaviors), numel(freq_fit));
disp(['FOOOF periodic component, ' band_name])
for b = 1:numel(all_behaviors)
    linear_self = mean_trace_spectrum(all_psth_onset_self{1, b}, ...
        all_electrode_index_self{b}, pre_mask, freq_fit, 1 / bin_size);
    linear_other = mean_trace_spectrum(all_psth_onset_other{1, b}, ...
        all_electrode_index_other{b}, pre_mask, freq_fit, 1 / bin_size);
    spec_self(b, :) = periodic_component_fooof(freq_fit, linear_self, ...
        fooof_aperiodic_mode, fooof_peak_widths);
    spec_other(b, :) = periodic_component_fooof(freq_fit, linear_other, ...
        fooof_aperiodic_mode, fooof_peak_widths);
end
spec_self = spec_self(:, freq_fit <= plot_max_hz);
spec_other = spec_other(:, freq_fit <= plot_max_hz);

if isempty(freq_hz)
    disp(['no pre-onset samples for ' band_name])
else
    rows_play_nonplay = [spec_self(play_behaviors, :); spec_other(play_behaviors, :); ...
        spec_self(other_behaviors, :); spec_other(other_behaviors, :)];
    labels_play_nonplay = [strcat(all_behaviors(play_behaviors), ' self'), ...
        strcat(all_behaviors(play_behaviors), ' other'), ...
        strcat(all_behaviors(other_behaviors), ' self'), ...
        strcat(all_behaviors(other_behaviors), ' other')];
    figure
    name_psth_figure(band_name, 'onset_pre_spectrum_play_nonplay_periodic')
    plot_spectrum_stack(freq_hz, rows_play_nonplay, labels_play_nonplay, ...
        1:(2 * numel(play_behaviors)), ...
        (2 * numel(play_behaviors) + 1):size(rows_play_nonplay, 1), ...
        [0.85 0.15 0.15], [0.15 0.15 0.15], 'Play', 'Non-play')

    rows_self_other = [spec_self(play_behaviors, :); spec_other(play_behaviors, :)];
    labels_self_other = [strcat(all_behaviors(play_behaviors), ' self'), ...
        strcat(all_behaviors(play_behaviors), ' other')];
    figure
    name_psth_figure(band_name, 'onset_pre_spectrum_self_other_play_periodic')
    plot_spectrum_stack(freq_hz, rows_self_other, labels_self_other, ...
        1:numel(play_behaviors), numel(play_behaviors) + (1:numel(play_behaviors)), ...
        [0.10 0.25 0.80], [0.80 0.15 0.15], 'Self play', 'Other play')
end

%% Per-animal waveforms and spectra
% Each animal mean uses only that animal's electrode-1 trials. The pooled
% line uses every electrode-1 trial. Every spectrum below is the FOOOF
% periodic component of that waveform, the same peak_fit as the stack above.
n_animals = numel(animal_label);
n_freq = numel(freq_hz);
animal_wave_self = nan(numel(all_behaviors), n_animals, n_pre);
animal_wave_other = nan(numel(all_behaviors), n_animals, n_pre);
animal_spec_self = nan(numel(all_behaviors), n_animals, n_freq);
animal_spec_other = nan(numel(all_behaviors), n_animals, n_freq);
trial_wave_self = nan(numel(all_behaviors), n_pre);
trial_wave_other = nan(numel(all_behaviors), n_pre);
trial_spec_self = nan(numel(all_behaviors), n_freq);
trial_spec_other = nan(numel(all_behaviors), n_freq);
n_self = zeros(numel(all_behaviors), n_animals);
n_other = zeros(numel(all_behaviors), n_animals);
sample_rate = 1 / bin_size;

for b = 1:numel(all_behaviors)
    [trial_wave_self(b, :), ~] = preonset_mean_trace(all_psth_onset_self{1, b}, ...
        all_electrode_index_self{b}, [], [], pre_mask);
    [trial_wave_other(b, :), ~] = preonset_mean_trace(all_psth_onset_other{1, b}, ...
        all_electrode_index_other{b}, [], [], pre_mask);
    trial_spec_self(b, :) = periodic_spectrum(trial_wave_self(b, :), sample_rate, ...
        freq_fit, fooof_aperiodic_mode, fooof_peak_widths, plot_max_hz);
    trial_spec_other(b, :) = periodic_spectrum(trial_wave_other(b, :), sample_rate, ...
        freq_fit, fooof_aperiodic_mode, fooof_peak_widths, plot_max_hz);
    for animal_n = 1:n_animals
        [animal_wave_self(b, animal_n, :), n_self(b, animal_n)] = preonset_mean_trace( ...
            all_psth_onset_self{1, b}, all_electrode_index_self{b}, ...
            all_session_index_self{b}, animal_n, pre_mask);
        [animal_wave_other(b, animal_n, :), n_other(b, animal_n)] = preonset_mean_trace( ...
            all_psth_onset_other{1, b}, all_electrode_index_other{b}, ...
            all_session_index_other{b}, animal_n, pre_mask);
        animal_spec_self(b, animal_n, :) = periodic_spectrum( ...
            squeeze(animal_wave_self(b, animal_n, :)).', sample_rate, ...
            freq_fit, fooof_aperiodic_mode, fooof_peak_widths, plot_max_hz);
        animal_spec_other(b, animal_n, :) = periodic_spectrum( ...
            squeeze(animal_wave_other(b, animal_n, :)).', sample_rate, ...
            freq_fit, fooof_aperiodic_mode, fooof_peak_widths, plot_max_hz);
    end
end

[equal_wave_self, loo_wave_self] = reference_waves(animal_wave_self, trial_wave_self, n_self);
[equal_wave_other, loo_wave_other] = reference_waves(animal_wave_other, trial_wave_other, n_other);
equal_spec_self = periodic_spectra_of_rows(equal_wave_self, sample_rate, ...
    freq_fit, fooof_aperiodic_mode, fooof_peak_widths, plot_max_hz);
equal_spec_other = periodic_spectra_of_rows(equal_wave_other, sample_rate, ...
    freq_fit, fooof_aperiodic_mode, fooof_peak_widths, plot_max_hz);
loo_spec_self = periodic_spectra_of_rows(loo_wave_self, sample_rate, ...
    freq_fit, fooof_aperiodic_mode, fooof_peak_widths, plot_max_hz);
loo_spec_other = periodic_spectra_of_rows(loo_wave_other, sample_rate, ...
    freq_fit, fooof_aperiodic_mode, fooof_peak_widths, plot_max_hz);

fprintf('\n%s pre-onset, electrode 1. Scores are the mean FOOOF periodic component from %.0f to %.0f Hz.\n', ...
    band_name, effect_band(1), effect_band(2));
print_actor_table('implanted', all_behaviors, animal_label, n_self, ...
    animal_wave_self, animal_spec_self, trial_wave_self, trial_spec_self, equal_spec_self, ...
    freq_hz, effect_band, sample_rate, freq_fit, fooof_aperiodic_mode, fooof_peak_widths, plot_max_hz);
print_actor_table('partner', all_behaviors, animal_label, n_other, ...
    animal_wave_other, animal_spec_other, trial_wave_other, trial_spec_other, equal_spec_other, ...
    freq_hz, effect_band, sample_rate, freq_fit, fooof_aperiodic_mode, fooof_peak_widths, plot_max_hz);

%% Per-animal pre-onset waveforms, implanted animal
open_animal_pages(band_name, 'onset_pre_waveform_animals_self', time_pre, ...
    animal_wave_self, trial_wave_self, equal_wave_self, loo_wave_self, n_self, ...
    animal_label, animal_colors, all_behaviors, behavior_pages, behavior_page_names, ...
    'Time from onset (s)', 'Log power', 'Implanted animal', []);

%% Per-animal pre-onset waveforms, partner
open_animal_pages(band_name, 'onset_pre_waveform_animals_other', time_pre, ...
    animal_wave_other, trial_wave_other, equal_wave_other, loo_wave_other, n_other, ...
    animal_label, animal_colors, all_behaviors, behavior_pages, behavior_page_names, ...
    'Time from onset (s)', 'Log power', 'Partner', []);

%% Per-animal pre-onset spectra, implanted animal
open_animal_pages(band_name, 'onset_pre_spectrum_animals_self', freq_hz, ...
    animal_spec_self, trial_spec_self, equal_spec_self, loo_spec_self, n_self, ...
    animal_label, animal_colors, all_behaviors, behavior_pages, behavior_page_names, ...
    'Frequency (Hz)', 'Periodic power (log10)', 'Implanted animal', [1 5]);

%% Per-animal pre-onset spectra, partner
open_animal_pages(band_name, 'onset_pre_spectrum_animals_other', freq_hz, ...
    animal_spec_other, trial_spec_other, equal_spec_other, loo_spec_other, n_other, ...
    animal_label, animal_colors, all_behaviors, behavior_pages, behavior_page_names, ...
    'Frequency (Hz)', 'Periodic power (log10)', 'Partner', [1 5]);

%% Raw mean versus per-trial z-score
% Left axis: mean of the raw log power. Right axis: mean after each trial
% is z-scored on its own -2 to 0 s window. Both use electrode 1.
baseline_mask = time_psth >= -2 & time_psth <= 0;
raw_self = nan(numel(all_behaviors), n_pre);
raw_other = nan(numel(all_behaviors), n_pre);
z_self = nan(numel(all_behaviors), n_pre);
z_other = nan(numel(all_behaviors), n_pre);
for b = 1:numel(all_behaviors)
    raw_self(b, :) = raw_pre_mean(all_psth_onset_self{1, b}, ...
        all_electrode_index_self{b}, pre_mask);
    raw_other(b, :) = raw_pre_mean(all_psth_onset_other{1, b}, ...
        all_electrode_index_other{b}, pre_mask);
    z_self(b, :) = baseline_z_mean(all_psth_onset_self{1, b}, ...
        all_electrode_index_self{b}, baseline_mask, pre_mask);
    z_other(b, :) = baseline_z_mean(all_psth_onset_other{1, b}, ...
        all_electrode_index_other{b}, baseline_mask, pre_mask);
end

open_scale_pages(band_name, 'onset_pre_raw_vs_trial_zscore_self', time_pre, ...
    raw_self, z_self, all_behaviors, behavior_pages, behavior_page_names, 'Implanted animal');
open_scale_pages(band_name, 'onset_pre_raw_vs_trial_zscore_other', time_pre, ...
    raw_other, z_other, all_behaviors, behavior_pages, behavior_page_names, 'Partner');
end

function name_psth_figure(band_name, plot_id)
set(gcf, 'WindowState', 'maximized', 'NumberTitle', 'off', ...
    'Name', ['psth_' band_name '_' plot_id]);
end

function spec = mean_trace_spectrum(psth, elec, pre_mask, freq_out, fs)
[trace, ~] = preonset_mean_trace(psth, elec, [], [], pre_mask);
spec = spectrum_of_trace(trace, fs, freq_out);
end

function [trace, n_trials] = preonset_mean_trace(psth, elec, session_index, animal_id, pre_mask)
% Electrode-1 trials, each demeaned on the pre-onset samples, then averaged.
% animal_id empty keeps every animal. A scalar keeps that animal only.
trace = nan(1, nnz(pre_mask));
n_trials = 0;
if isempty(psth) || isempty(elec) || size(psth, 2) ~= numel(pre_mask)
    return
end
n_use = min(numel(elec), size(psth, 1));
psth = psth(1:n_use, :);
elec = elec(1:n_use);
keep = elec(:) == 1;
if ~isempty(animal_id)
    n_use = min(n_use, numel(session_index));
    psth = psth(1:n_use, :);
    keep = keep(1:n_use) & session_index(1:n_use) == animal_id;
end
psth = psth(keep, :);
n_trials = size(psth, 1);
if n_trials == 0 || size(psth, 2) < 8
    return
end
seg = psth(:, pre_mask);
seg = seg - mean(seg, 2, 'omitnan');
trace = mean(seg, 1, 'omitnan');
if all(~isfinite(trace))
    trace = nan(1, nnz(pre_mask));
    return
end
trace(~isfinite(trace)) = 0;
end

function spec = spectrum_of_trace(trace, fs, freq_out)
spec = nan(1, numel(freq_out));
trace = double(trace(:).');
if isempty(trace) || numel(trace) < 8 || all(~isfinite(trace))
    return
end
trace(~isfinite(trace)) = 0;
n = numel(trace);
transformed = fft(trace);
n_keep = floor(n / 2);
pow = abs(transformed(2:n_keep + 1)).^2 / n;
freq = (1:n_keep) * (fs / n);
spec = interp1(freq, pow, freq_out, 'linear', nan);
end

function value = band_mean(freq, spec, band)
keep = freq >= band(1) & freq <= band(2) & isfinite(spec);
if ~any(keep)
    value = NaN;
else
    value = mean(spec(keep));
end
end

function r = shape_correlation(a, b)
a = double(a(:));
b = double(b(:));
keep = isfinite(a) & isfinite(b);
r = NaN;
if nnz(keep) < 8
    return
end
a = a(keep) - mean(a(keep));
b = b(keep) - mean(b(keep));
denominator = sqrt(sum(a .^ 2) * sum(b .^ 2));
if denominator > 0
    r = sum(a .* b) / denominator;
end
end

function print_actor_table(actor_name, behavior_names, animal_names, n_trials, ...
        animal_wave, animal_spec, trial_wave, trial_spec, equal_spec, ...
        freq, band, fs, freq_fit, aperiodic_mode, peak_widths, plot_max_hz)
fprintf('\n%s\n', actor_name);
n_animals = numel(animal_names);
for behavior_n = 1:numel(behavior_names)
    total = sum(n_trials(behavior_n, :));
    pooled_power = band_mean(freq, trial_spec(behavior_n, :), band);
    block = reshape(animal_wave(behavior_n, :, :), n_animals, []);
    equal_power = band_mean(freq, equal_spec(behavior_n, :), band);
    fprintf('  %s   trials %d   pooled %.4g   animals-equal %.4g\n', ...
        behavior_names{behavior_n}, total, pooled_power, equal_power);
    if total == 0
        continue
    end
    for animal_n = 1:n_animals
        n_animal = n_trials(behavior_n, animal_n);
        if n_animal == 0
            continue
        end
        own_power = band_mean(freq, squeeze(animal_spec(behavior_n, animal_n, :)).', band);
        correlation = shape_correlation(squeeze(animal_wave(behavior_n, animal_n, :)), ...
            trial_wave(behavior_n, :));
        remaining = NaN;
        if n_animal == total
            remaining = 0;
        elseif total > n_animal
            left_out = (total * trial_wave(behavior_n, :) - ...
                n_animal * block(animal_n, :)) / (total - n_animal);
            left_periodic = periodic_spectrum(left_out, fs, freq_fit, ...
                aperiodic_mode, peak_widths, plot_max_hz);
            left_power = band_mean(freq, left_periodic, band);
            if pooled_power > 0 && isfinite(left_power)
                remaining = left_power / pooled_power;
            end
        end
        note = '';
        if isfinite(remaining) && remaining < 0.5
            note = '  carries the 2-3 Hz peak';
        end
        fprintf('    %-6s  n=%4d  share=%.2f  r=% .2f  own=%.4g  remaining=%.2f%s\n', ...
            animal_names{animal_n}, n_animal, n_animal / total, correlation, ...
            own_power, remaining, note);
    end
end
end

function [equal_wave, loo_wave] = reference_waves(animal_wave, trial_wave, n_trials)
n_behaviors = size(trial_wave, 1);
n_animals = size(n_trials, 2);
n_time = size(trial_wave, 2);
equal_wave = nan(n_behaviors, n_time);
loo_wave = nan(n_behaviors, n_time);
for behavior_n = 1:n_behaviors
    block = reshape(animal_wave(behavior_n, :, :), n_animals, n_time);
    present = n_trials(behavior_n, :) > 0;
    if any(present)
        equal_wave(behavior_n, :) = mean(block(present, :), 1, 'omitnan');
    end
    total = sum(n_trials(behavior_n, :));
    [n_animal, largest] = max(n_trials(behavior_n, :));
    if total > n_animal && n_animal > 0
        loo_wave(behavior_n, :) = (total * trial_wave(behavior_n, :) - ...
            n_animal * block(largest, :)) / (total - n_animal);
    end
end
end

function periodic = periodic_spectrum(trace, fs, freq_fit, aperiodic_mode, peak_widths, plot_max_hz)
keep = freq_fit <= plot_max_hz;
periodic = nan(1, nnz(keep));
if isempty(freq_fit)
    return
end
linear = spectrum_of_trace(trace, fs, freq_fit);
full = periodic_component_fooof(freq_fit, linear, aperiodic_mode, peak_widths);
if numel(full) == numel(freq_fit)
    periodic = full(keep);
end
end

function specs = periodic_spectra_of_rows(waves, fs, freq_fit, aperiodic_mode, peak_widths, plot_max_hz)
n_plot = nnz(freq_fit <= plot_max_hz);
specs = nan(size(waves, 1), n_plot);
for row = 1:size(waves, 1)
    specs(row, :) = periodic_spectrum(waves(row, :), fs, freq_fit, ...
        aperiodic_mode, peak_widths, plot_max_hz);
end
end

function open_animal_pages(band_name, plot_prefix, x, animal_curves, trial_curve, ...
        equal_curve, loo_curve, n_trials, animal_names, colors, behavior_names, ...
        pages, page_names, x_label, y_label, actor_title, shade_band)
for page_n = 1:numel(pages)
    idx = pages{page_n};
    figure
    name_psth_figure(band_name, [plot_prefix '_' page_names{page_n}])
    plot_animal_panels(x, animal_curves(idx, :, :), trial_curve(idx, :), ...
        equal_curve(idx, :), loo_curve(idx, :), n_trials(idx, :), ...
        animal_names, colors, behavior_names(idx), x_label, y_label, ...
        actor_title, shade_band);
end
end

function open_scale_pages(band_name, plot_prefix, time_pre, raw_mean, z_mean, ...
        behavior_names, pages, page_names, actor_title)
for page_n = 1:numel(pages)
    idx = pages{page_n};
    figure
    name_psth_figure(band_name, [plot_prefix '_' page_names{page_n}])
    plot_scale_panels(time_pre, raw_mean(idx, :), z_mean(idx, :), ...
        behavior_names(idx), actor_title);
end
end

function [n_row, n_col] = panel_grid(n_panels)
if n_panels <= 4
    n_row = 2;
    n_col = 3;
else
    n_row = 2;
    n_col = 4;
end
end

function plot_animal_panels(x, animal_curves, trial_curve, equal_curve, loo_curve, n_trials, ...
        animal_names, colors, behavior_names, x_label, y_label, actor_title, shade_band)
n_animals = numel(animal_names);
n_behaviors = numel(behavior_names);
[n_row, n_col] = panel_grid(n_behaviors);
handles = gobjects(n_animals + 3, 1);
for behavior_n = 1:n_behaviors
    subplot(n_row, n_col, behavior_n)
    hold on
    for animal_n = 1:n_animals
        if n_trials(behavior_n, animal_n) == 0
            continue
        end
        y = squeeze(animal_curves(behavior_n, animal_n, :));
        h = plot(x, y(:).', 'Color', colors(animal_n, :));
        if ~isgraphics(handles(animal_n))
            set(h, 'DisplayName', animal_names{animal_n})
            handles(animal_n) = h;
        end
    end
    if any(isfinite(trial_curve(behavior_n, :)))
        h = plot(x, trial_curve(behavior_n, :), 'k', 'LineWidth', 1.6);
        if ~isgraphics(handles(n_animals + 1))
            set(h, 'DisplayName', 'All trials')
            handles(n_animals + 1) = h;
        end
    end
    if any(isfinite(equal_curve(behavior_n, :)))
        h = plot(x, equal_curve(behavior_n, :), 'Color', [0.45 0.45 0.45], 'LineWidth', 1.6);
        if ~isgraphics(handles(n_animals + 2))
            set(h, 'DisplayName', 'Mean of animals')
            handles(n_animals + 2) = h;
        end
    end
    if any(isfinite(loo_curve(behavior_n, :)))
        h = plot(x, loo_curve(behavior_n, :), 'k--', 'LineWidth', 1.2);
        if ~isgraphics(handles(n_animals + 3))
            set(h, 'DisplayName', 'Without largest')
            handles(n_animals + 3) = h;
        end
    end
    if ~isempty(shade_band) && numel(x) > 1
        yl = ylim;
        shade = patch([shade_band(1) shade_band(2) shade_band(2) shade_band(1)], ...
            [yl(1) yl(1) yl(2) yl(2)], [0.75 0.75 0.75], ...
            'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        uistack(shade, 'bottom')
    end
    hold off
    title(behavior_names{behavior_n}, 'Interpreter', 'none')
    if behavior_n == 1
        xlabel(x_label)
        ylabel(y_label)
    end
    if numel(x) > 1
        xlim([x(1), x(end)])
    end
    set(gca, 'TickDir', 'out')
end
subplot(n_row, n_col, n_behaviors + 1)
axis off
present = handles(isgraphics(handles));
if ~isempty(present)
    legend(present, 'Location', 'northwest')
end
title(actor_title)
end

function trace = raw_pre_mean(psth, elec, pre_mask)
trace = nan(1, nnz(pre_mask));
if isempty(psth) || isempty(elec) || size(psth, 2) ~= numel(pre_mask)
    return
end
n_use = min(numel(elec), size(psth, 1));
psth = psth(1:n_use, :);
psth = psth(elec(1:n_use) == 1, :);
if isempty(psth)
    return
end
trace = mean(psth(:, pre_mask), 1, 'omitnan');
end

function trace = baseline_z_mean(psth, elec, baseline_mask, pre_mask)
trace = nan(1, nnz(pre_mask));
if isempty(psth) || isempty(elec) || size(psth, 2) ~= numel(pre_mask)
    return
end
if numel(baseline_mask) ~= size(psth, 2)
    return
end
n_use = min(numel(elec), size(psth, 1));
psth = psth(1:n_use, :);
psth = psth(elec(1:n_use) == 1, :);
if isempty(psth)
    return
end
center = mean(psth(:, baseline_mask), 2);
scale = std(psth(:, baseline_mask), 0, 2);
scale(scale == 0) = NaN;
valid = ~isnan(scale) & ~isnan(center);
if ~any(valid)
    return
end
z = (psth(valid, :) - center(valid)) ./ scale(valid);
trace = mean(z(:, pre_mask), 1, 'omitnan');
end

function plot_scale_panels(time_pre, raw_mean, z_mean, behavior_names, actor_title)
n_behaviors = numel(behavior_names);
[n_row, n_col] = panel_grid(n_behaviors);
for behavior_n = 1:n_behaviors
    subplot(n_row, n_col, behavior_n)
    yyaxis left
    plot(time_pre, raw_mean(behavior_n, :), 'k', 'LineWidth', 1.2)
    ylabel('Raw log power')
    yyaxis right
    plot(time_pre, z_mean(behavior_n, :), 'Color', [0.10 0.25 0.80], 'LineWidth', 1.2)
    ylabel('Trial z-score')
    title(behavior_names{behavior_n}, 'Interpreter', 'none')
    if behavior_n == 1
        xlabel('Time from onset (s)')
    end
    if numel(time_pre) > 1
        xlim([time_pre(1), time_pre(end)])
    end
    set(gca, 'TickDir', 'out')
end
subplot(n_row, n_col, n_behaviors + 1)
axis off
text(0, 0.6, actor_title, 'FontWeight', 'bold')
text(0, 0.35, {'Black: raw mean', 'Blue: each trial z-scored', 'on its own -2 to 0 s'})
end

function plot_spectrum_stack(freq, rows, labels, idx_a, idx_b, color_a, color_b, name_a, name_b)
ax2 = subplot(3, 1, 2:3);
hold on
draw_group_spectrum(freq, rows(idx_b, :), color_b)
draw_group_spectrum(freq, rows(idx_a, :), color_a)
plot(freq([1 end]), [0 0], 'k:')
axis tight
yl = ylim;
delta_patch = patch([1 5 5 1], [yl(1) yl(1) yl(2) yl(2)], [0.75 0.75 0.75], ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack(delta_patch, 'bottom')
xlim([freq(1) freq(end)])
xlabel('Frequency (Hz)')
ylabel('Periodic power (log10)')
title([name_a ' vs ' name_b])
set(gca, 'TickDir', 'out')

ax1 = subplot(3, 1, 1);
imagesc(freq, 1:size(rows, 1), rows)
axis xy
finite_vals = rows(isfinite(rows));
if ~isempty(finite_vals)
    lo = prctile(finite_vals, 5);
    hi = prctile(finite_vals, 95);
    if hi > lo
        clim([lo hi])
    end
end
hold on
plot(freq([1 end]), [1 1] * numel(idx_a), 'w', 'LineWidth', 2)
yticks(1:size(rows, 1))
yticklabels(labels)
set(gca, 'TickDir', 'out')
colormap(ax1, jet(256))
colorbar
ax2.Position(1) = ax1.Position(1);
ax2.Position(3) = ax1.Position(3);
end

function draw_group_spectrum(freq, rows, color)
rows = rows(all(isfinite(rows), 2), :);
if isempty(rows)
    return
end
if size(rows, 1) >= 2
    [~, ~, ci] = ttest(rows);
    fill([freq fliplr(freq)], [ci(1, :) fliplr(ci(2, :))], color, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
end
plot(freq, mean(rows, 1), 'Color', color, 'LineWidth', 1.5)
end

function periodic = periodic_component_fooof(freq, linear_power, aperiodic_mode, peak_widths)
periodic = nan(1, numel(freq));
freq = freq(:).';
linear_power = linear_power(:).';
if numel(freq) < 8 || numel(linear_power) ~= numel(freq)
    return
end
if any(~isfinite(freq)) || any(~isfinite(linear_power))
    return
end
if all(linear_power <= 0)
    periodic = zeros(1, numel(freq));
    return
end
linear_power = max(linear_power, realmin('double'));
out = fooof_matlab(freq(:), linear_power(:), ...
    'AperiodicMode', aperiodic_mode, ...
    'PeakWidthLimits', peak_widths);
periodic = out.peak_fit(:).';
end
