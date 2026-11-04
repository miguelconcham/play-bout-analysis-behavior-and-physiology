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
%     One column per behavior, four or six columns per figure (play 1-4,
%     play 5-8, then the six non-play behaviors). The upper image spans two
%     rows and has one row per animal. The bottom panel draws each animal
%     as a thin line, the trial-weighted mean as a thick black line, and
%     the equal-animal mean as a thick gray line.
% %% Per-animal spectra
%     Same columns. The image and the lines are the FOOOF periodic component.
%     The fit and the panels use 0-20 Hz. The Command Window reports the
%     mean periodic component from 2 to 3 Hz. 20 Hz is below the Nyquist
%     frequency of the 10 ms power step.
% %% Pooled play and non-play
%     Self and partner are pooled. Each row is one behavior from one animal.
%     Within each animal those behaviors are averaged, then the thick line
%     is the mean across animals. The band is the 95% t interval across
%     animals, so a second behavior from the same animal is a repeated
%     measure and does not narrow the interval by itself. The same figure
%     is drawn for the log-power traces and for the periodic spectra.
%     Scratching stays out of non-play.
% %% Trial-level circular-shift null
%     Each -5-to-0 s trial is independently circularly rotated. This keeps
%     that trial's spectrum exactly and removes common phase relative to
%     onset. For every animal x behavior x actor unit, the trial mean is
%     transformed to a spectrum and FOOOF. The log10 spectrum minus the
%     fitted aperiodic component is ranked from 0 to 1 at every frequency
%     against 1,000 circular null fits. Unit percentiles are then averaged
%     within animal and animals are averaged equally. Target-band,
%     corrected-peak, and leave-one-animal-out summaries use the same hierarchy.
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
% A whole-session z-score is not in these files. The scale comparison is
% only the per-trial baseline z-score. Scratching is left out of the
% non-play spectrum stack and pooled controls, as in the grand-average script.

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
% Circular-shift inference. The target band is predefined above. Peak
% searches are corrected by the permutation maximum inside search_band.
circular_null_n_perm = 1000;
circular_null_seed = 29;
circular_null_min_trials_per_unit = 5;
circular_null_min_units_per_animal = 2;
circular_null_search_band = [0.4 8];
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

[equal_wave_self, ~] = reference_waves(animal_wave_self, trial_wave_self, n_self);
[equal_wave_other, ~] = reference_waves(animal_wave_other, trial_wave_other, n_other);
equal_spec_self = periodic_spectra_of_rows(equal_wave_self, sample_rate, ...
    freq_fit, fooof_aperiodic_mode, fooof_peak_widths, plot_max_hz);
equal_spec_other = periodic_spectra_of_rows(equal_wave_other, sample_rate, ...
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
    animal_wave_self, trial_wave_self, equal_wave_self, n_self, ...
    animal_label, animal_colors, all_behaviors, behavior_pages, behavior_page_names, ...
    'Time from onset (s)', 'Log power', 'Implanted animal', []);

%% Per-animal pre-onset waveforms, partner
open_animal_pages(band_name, 'onset_pre_waveform_animals_other', time_pre, ...
    animal_wave_other, trial_wave_other, equal_wave_other, n_other, ...
    animal_label, animal_colors, all_behaviors, behavior_pages, behavior_page_names, ...
    'Time from onset (s)', 'Log power', 'Partner', []);

%% Per-animal pre-onset spectra, implanted animal
open_animal_pages(band_name, 'onset_pre_spectrum_animals_self', freq_hz, ...
    animal_spec_self, trial_spec_self, equal_spec_self, n_self, ...
    animal_label, animal_colors, all_behaviors, behavior_pages, behavior_page_names, ...
    'Frequency (Hz)', 'Periodic power (log10)', 'Implanted animal', [1 5]);

%% Per-animal pre-onset spectra, partner
open_animal_pages(band_name, 'onset_pre_spectrum_animals_other', freq_hz, ...
    animal_spec_other, trial_spec_other, equal_spec_other, n_other, ...
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

%% Pooled play and non-play, repeated measures
% One row per animal and behavior, self and partner together. Behaviors are
% averaged within the animal. The line and the band are the mean and the
% 95% t interval across animals.
[play_wave, play_ids] = stack_actor_behaviors(animal_wave_self, animal_wave_other, ...
    n_self, n_other, play_behaviors);
[nonplay_wave, nonplay_ids] = stack_actor_behaviors(animal_wave_self, animal_wave_other, ...
    n_self, n_other, other_behaviors);
[play_spec, play_spec_ids] = stack_actor_behaviors(animal_spec_self, animal_spec_other, ...
    n_self, n_other, play_behaviors);
[nonplay_spec, nonplay_spec_ids] = stack_actor_behaviors(animal_spec_self, animal_spec_other, ...
    n_self, n_other, other_behaviors);

[play_mu, play_lo, play_hi] = repeated_animal_mean(play_wave, play_ids);
[nonplay_mu, nonplay_lo, nonplay_hi] = repeated_animal_mean(nonplay_wave, nonplay_ids);
[play_spec_mu, play_spec_lo, play_spec_hi] = repeated_animal_mean(play_spec, play_spec_ids);
[nonplay_spec_mu, nonplay_spec_lo, nonplay_spec_hi] = repeated_animal_mean(nonplay_spec, nonplay_spec_ids);

figure
name_psth_figure(band_name, 'onset_pre_waveform_play_pooled')
plot_pooled_repeated(time_pre, play_wave, play_ids, animal_label, ...
    play_mu, play_lo, play_hi, 'Time from onset (s)', 'Log power', ...
    'Play, self and partner', []);
figure
name_psth_figure(band_name, 'onset_pre_waveform_nonplay_pooled')
plot_pooled_repeated(time_pre, nonplay_wave, nonplay_ids, animal_label, ...
    nonplay_mu, nonplay_lo, nonplay_hi, 'Time from onset (s)', 'Log power', ...
    'Non-play, self and partner', []);
figure
name_psth_figure(band_name, 'onset_pre_spectrum_play_pooled')
plot_pooled_repeated(freq_hz, play_spec, play_spec_ids, animal_label, ...
    play_spec_mu, play_spec_lo, play_spec_hi, 'Frequency (Hz)', ...
    'Periodic power (log10)', 'Play, self and partner', [1 5]);
figure
name_psth_figure(band_name, 'onset_pre_spectrum_nonplay_pooled')
plot_pooled_repeated(freq_hz, nonplay_spec, nonplay_spec_ids, animal_label, ...
    nonplay_spec_mu, nonplay_spec_lo, nonplay_spec_hi, 'Frequency (Hz)', ...
    'Periodic power (log10)', 'Non-play, self and partner', [1 5]);

%% Trial-level circular-shift FOOOF null for onset locking
% Each complete -5-to-0 s trial is demeaned, then independently circularly
% shifted. Thus every trial keeps exactly the same Fourier magnitude while
% its phase relative to onset is randomized. FOOOF is fitted separately to
% every observed and permuted animal x behavior x actor mean. The statistic
% is log10(real spectrum) minus the final fitted aperiodic component, ranked
% against that unit's null before pooling. There is no inferential test on
% the mean theta-power waveform.
previous_rng = rng;
rng(circular_null_seed, 'twister');
disp(['trial-level circular-shift FOOOF null, ' ...
    num2str(circular_null_n_perm) ' permutations'])
play_null = circular_trial_null(all_psth_onset_self, all_psth_onset_other, ...
    all_electrode_index_self, all_electrode_index_other, ...
    all_session_index_self, all_session_index_other, ...
    play_behaviors, numel(animal_label), pre_mask, sample_rate, freq_hz, ...
    effect_band, circular_null_search_band, circular_null_n_perm, ...
    circular_null_min_trials_per_unit, circular_null_min_units_per_animal, ...
    fooof_aperiodic_mode, fooof_peak_widths);
nonplay_null = circular_trial_null(all_psth_onset_self, all_psth_onset_other, ...
    all_electrode_index_self, all_electrode_index_other, ...
    all_session_index_self, all_session_index_other, ...
    other_behaviors, numel(animal_label), pre_mask, sample_rate, freq_hz, ...
    effect_band, circular_null_search_band, circular_null_n_perm, ...
    circular_null_min_trials_per_unit, circular_null_min_units_per_animal, ...
    fooof_aperiodic_mode, fooof_peak_widths);
rng(previous_rng);

print_circular_results([band_name ' play'], play_null, animal_label, ...
    all_behaviors, effect_band);
print_circular_results([band_name ' non-play'], nonplay_null, animal_label, ...
    all_behaviors, effect_band);

figure
name_psth_figure(band_name, 'onset_pre_spectrum_circular_null')
plot_circular_spectrum_hierarchy(freq_hz, play_null, nonplay_null, ...
    animal_label, effect_band);
end

function result = circular_trial_null(psth_self, psth_other, elec_self, elec_other, ...
        animal_self, animal_other, behavior_index, n_animals, pre_mask, fs, ...
        freq, effect_band, search_band, n_perm, min_trials_per_unit, ...
        min_units_per_animal, aperiodic_mode, peak_widths)
% One unit is one animal x behavior x actor. FOOOF is fitted before pooling.
units = struct('segments', {}, 'animal', {}, 'behavior', {}, ...
    'actor', {}, 'n_trials', {});
n_units = zeros(n_animals, 1);
n_trials = zeros(n_animals, 1);
for animal_n = 1:n_animals
    for behavior_n = behavior_index
        self_seg = trial_pre_segments(psth_self{1, behavior_n}, ...
            elec_self{behavior_n}, animal_self{behavior_n}, animal_n, pre_mask);
        if size(self_seg, 1) >= min_trials_per_unit
            unit_n = numel(units) + 1;
            units(unit_n).segments = self_seg;
            units(unit_n).animal = animal_n;
            units(unit_n).behavior = behavior_n;
            units(unit_n).actor = 'self';
            units(unit_n).n_trials = size(self_seg, 1);
            n_trials(animal_n) = n_trials(animal_n) + size(self_seg, 1);
        end
        other_seg = trial_pre_segments(psth_other{1, behavior_n}, ...
            elec_other{behavior_n}, animal_other{behavior_n}, animal_n, pre_mask);
        if size(other_seg, 1) >= min_trials_per_unit
            unit_n = numel(units) + 1;
            units(unit_n).segments = other_seg;
            units(unit_n).animal = animal_n;
            units(unit_n).behavior = behavior_n;
            units(unit_n).actor = 'other';
            units(unit_n).n_trials = size(other_seg, 1);
            n_trials(animal_n) = n_trials(animal_n) + size(other_seg, 1);
        end
    end
    n_units(animal_n) = sum([units.animal] == animal_n);
end

n_total_units = numel(units);
n_freq = numel(freq);
unit_periodic = nan(n_total_units, n_freq);
unit_percentile = nan(n_total_units, n_freq);
unit_null_percentile = nan(n_total_units, n_perm, n_freq);
unit_target_percentile = nan(n_total_units, 1);
unit_null_target_percentile = nan(n_total_units, n_perm);
unit_peak_hz = nan(n_total_units, 1);
unit_peak_p = nan(n_total_units, 1);

for unit_n = 1:n_total_units
    observed_wave = mean(units(unit_n).segments, 1);
    observed_raw = spectrum_of_trace(observed_wave, fs, freq);
    observed_periodic = periodic_residual_fooof( ...
        freq, observed_raw, aperiodic_mode, peak_widths);
    null_periodic = nan(n_perm, n_freq);
    for permutation = 1:n_perm
        null_wave = circular_shift_mean(units(unit_n).segments);
        null_raw = spectrum_of_trace(null_wave, fs, freq);
        null_periodic(permutation, :) = periodic_residual_fooof( ...
            freq, null_raw, aperiodic_mode, peak_widths);
        if mod(permutation, 100) == 0 || permutation == n_perm
            fprintf('  unit %d/%d FOOOF null %d/%d\n', ...
                unit_n, n_total_units, permutation, n_perm);
        end
    end
    unit_periodic(unit_n, :) = observed_periodic;
    unit_percentile(unit_n, :) = percentile_spectrum( ...
        observed_periodic, null_periodic);
    null_percentiles = null_rank_percentiles(null_periodic);
    unit_null_percentile(unit_n, :, :) = ...
        reshape(null_percentiles, 1, n_perm, n_freq);
    observed_target = band_mean(freq, observed_periodic, effect_band);
    null_target = band_rows(freq, null_periodic, effect_band);
    unit_target_percentile(unit_n) = empirical_percentile( ...
        observed_target, null_target);
    unit_null_target_percentile(unit_n, :) = ...
        null_rank_percentiles(null_target).';
    [unit_peak_hz(unit_n), unit_peak_p(unit_n)] = score_peak_inference( ...
        freq, unit_percentile(unit_n, :), null_percentiles, search_band);
end

unit_animal = reshape([units.animal], [], 1);
unit_behavior = reshape([units.behavior], [], 1);
unit_actor = {units.actor}.';
unit_trials = reshape([units.n_trials], [], 1);
included = n_units >= min_units_per_animal;

animal_score = nan(n_animals, n_freq);
animal_percentile = nan(n_animals, n_freq);
animal_null_score = nan(n_animals, n_perm, n_freq);
animal_target_score = nan(n_animals, 1);
animal_target_percentile = nan(n_animals, 1);
animal_null_target_score = nan(n_animals, n_perm);
animal_p_effect = nan(n_animals, 1);
animal_peak_hz = nan(n_animals, 1);
animal_peak_p = nan(n_animals, 1);
for animal_n = 1:n_animals
    selected = unit_animal == animal_n;
    if ~any(selected)
        continue
    end
    animal_score(animal_n, :) = mean(unit_percentile(selected, :), 1, 'omitnan');
    null_scores = squeeze(mean(unit_null_percentile(selected, :, :), 1, 'omitnan'));
    null_scores = reshape(null_scores, n_perm, n_freq);
    animal_null_score(animal_n, :, :) = ...
        reshape(null_scores, 1, n_perm, n_freq);
    animal_percentile(animal_n, :) = percentile_spectrum( ...
        animal_score(animal_n, :), null_scores);
    animal_target_score(animal_n) = ...
        mean(unit_target_percentile(selected), 'omitnan');
    null_target_scores = mean( ...
        unit_null_target_percentile(selected, :), 1, 'omitnan').';
    animal_null_target_score(animal_n, :) = null_target_scores.';
    animal_target_percentile(animal_n) = empirical_percentile( ...
        animal_target_score(animal_n), null_target_scores);
    animal_p_effect(animal_n) = empirical_upper_p( ...
        animal_target_score(animal_n), null_target_scores);
    [animal_peak_hz(animal_n), animal_peak_p(animal_n)] = ...
        score_peak_inference(freq, animal_score(animal_n, :), ...
        null_scores, search_band);
end

group_score = nan(1, n_freq);
group_percentile = nan(1, n_freq);
group_null_score = nan(n_perm, n_freq);
group_target_score = NaN;
group_target_percentile = NaN;
group_null_target_score = nan(n_perm, 1);
p_effect = NaN;
group_peak_hz = NaN;
group_peak_p = NaN;
group_max_threshold = NaN;
significant_frequency = false(size(freq));
if nnz(included) >= 2
    group_score = mean(animal_score(included, :), 1, 'omitnan');
    group_null_score = squeeze(mean( ...
        animal_null_score(included, :, :), 1, 'omitnan'));
    group_null_score = reshape(group_null_score, n_perm, n_freq);
    group_percentile = percentile_spectrum(group_score, group_null_score);
    group_target_score = mean(animal_target_score(included), 'omitnan');
    group_null_target_score = mean( ...
        animal_null_target_score(included, :), 1, 'omitnan').';
    group_target_percentile = empirical_percentile( ...
        group_target_score, group_null_target_score);
    p_effect = empirical_upper_p(group_target_score, group_null_target_score);
    [group_peak_hz, group_peak_p, group_max_threshold, significant_frequency] = ...
        score_peak_inference(freq, group_score, group_null_score, search_band);
end

loo_target_score = nan(n_animals, 1);
loo_p_effect = nan(n_animals, 1);
for animal_n = 1:n_animals
    keep = included;
    keep(animal_n) = false;
    if nnz(keep) < 2
        continue
    end
    loo_target_score(animal_n) = mean(animal_target_score(keep), 'omitnan');
    loo_null = mean(animal_null_target_score(keep, :), 1, 'omitnan').';
    loo_p_effect(animal_n) = empirical_upper_p(loo_target_score(animal_n), loo_null);
end

result.unit_periodic = unit_periodic;
result.unit_percentile = unit_percentile;
result.unit_null_percentile = unit_null_percentile;
result.unit_target_percentile = unit_target_percentile;
result.unit_peak_hz = unit_peak_hz;
result.unit_peak_p = unit_peak_p;
result.unit_animal = unit_animal;
result.unit_behavior = unit_behavior;
result.unit_actor = unit_actor;
result.unit_trials = unit_trials;
result.animal_score = animal_score;
result.animal_percentile = animal_percentile;
result.animal_null_score = animal_null_score;
result.animal_target_score = animal_target_score;
result.animal_target_percentile = animal_target_percentile;
result.animal_null_target_score = animal_null_target_score;
result.animal_p_effect = animal_p_effect;
result.animal_peak_hz = animal_peak_hz;
result.animal_peak_p = animal_peak_p;
result.group_score = group_score;
result.group_percentile = group_percentile;
result.group_null_score = group_null_score;
result.group_target_score = group_target_score;
result.group_target_percentile = group_target_percentile;
result.group_null_target_score = group_null_target_score;
result.p_effect = p_effect;
result.group_peak_hz = group_peak_hz;
result.group_peak_p = group_peak_p;
result.group_max_threshold = group_max_threshold;
result.significant_frequency = significant_frequency;
result.loo_target_score = loo_target_score;
result.loo_p_effect = loo_p_effect;
result.included = included;
result.n_units = n_units;
result.n_trials = n_trials;
result.effect_band = effect_band;
result.search_band = search_band;
end

function [peak_hz, peak_p, max_threshold, significant_frequency] = ...
        score_peak_inference(freq, observed_score, null_score, search_band)
search = freq >= search_band(1) & freq <= search_band(2);
significant_frequency = false(size(freq));
peak_hz = NaN;
peak_p = NaN;
max_threshold = NaN;
if ~any(search)
    return
end
freq_search = freq(search);
[observed_max, index] = finite_max(observed_score(search));
if isfinite(observed_max)
    peak_hz = freq_search(index);
end
null_max = finite_row_max(null_score(:, search));
peak_p = empirical_upper_p(observed_max, null_max);
valid = null_max(isfinite(null_max));
if ~isempty(valid)
    max_threshold = prctile(valid, 95);
    significant_frequency(search) = observed_score(search) > max_threshold;
end
end

function percentiles = percentile_spectrum(observed, null_curves)
percentiles = nan(size(observed));
for frequency = 1:numel(observed)
    percentiles(frequency) = empirical_percentile( ...
        observed(frequency), null_curves(:, frequency));
end
end

function percentile = empirical_percentile(observed, null_values)
valid = null_values(isfinite(null_values));
if ~isfinite(observed) || isempty(valid)
    percentile = NaN;
    return
end
less = sum(valid < observed);
equal = sum(valid == observed);
percentile = (less + 0.5 * equal + 0.5) / (numel(valid) + 1);
end

function percentiles = null_rank_percentiles(values)
was_vector = isvector(values);
if was_vector
    values = values(:);
end
percentiles = nan(size(values));
for column = 1:size(values, 2)
    x = values(:, column);
    valid_index = find(isfinite(x));
    if isempty(valid_index)
        continue
    end
    [sorted, order] = sort(x(valid_index));
    rank_value = nan(size(sorted));
    first = 1;
    while first <= numel(sorted)
        last = first;
        while last < numel(sorted) && sorted(last + 1) == sorted(first)
            last = last + 1;
        end
        midrank = (first + last) / 2;
        rank_value(first:last) = (midrank - 0.5) / numel(sorted);
        first = last + 1;
    end
    restored = nan(size(sorted));
    restored(order) = rank_value;
    percentiles(valid_index, column) = restored;
end
if was_vector
    percentiles = percentiles(:);
end
end

function values = band_rows(freq, spectra, band)
keep = freq >= band(1) & freq <= band(2);
if ~any(keep)
    values = nan(size(spectra, 1), 1);
else
    values = mean(spectra(:, keep), 2, 'omitnan');
end
end

function p = empirical_upper_p(observed, null_values)
valid = null_values(isfinite(null_values));
if ~isfinite(observed) || isempty(valid)
    p = NaN;
    return
end
tolerance = max(1e-12, abs(observed) * 1e-12);
p = (1 + sum(valid >= observed - tolerance)) / (1 + numel(valid));
end

function [value, index] = finite_max(values)
values = values(:).';
valid = isfinite(values);
value = NaN;
index = NaN;
if ~any(valid)
    return
end
valid_index = find(valid);
[value, local_index] = max(values(valid));
index = valid_index(local_index);
end

function values = finite_row_max(matrix)
values = nan(size(matrix, 1), 1);
for row = 1:size(matrix, 1)
    values(row) = finite_max(matrix(row, :));
end
end

function segments = trial_pre_segments(psth, electrode, animal_id, selected_animal, pre_mask)
segments = zeros(0, nnz(pre_mask));
if isempty(psth) || isempty(electrode) || isempty(animal_id) || ...
        size(psth, 2) ~= numel(pre_mask)
    return
end
n_use = min([size(psth, 1), numel(electrode), numel(animal_id)]);
electrode = electrode(1:n_use);
animal_id = animal_id(1:n_use);
keep = electrode(:) == 1 & animal_id(:) == selected_animal;
segments = double(psth(keep, pre_mask));
segments = segments(all(isfinite(segments), 2), :);
if isempty(segments)
    return
end
segments = segments - mean(segments, 2);
end

function trace = circular_shift_mean(segments)
% Independent circular rotation of every trial. The periodogram of each
% row is unchanged; only phase alignment between trials is removed.
[n_trials, n_time] = size(segments);
shifts = randi(n_time, n_trials, 1) - 1;
indices = mod(bsxfun(@minus, 0:n_time - 1, shifts), n_time) + 1;
rows = repmat((1:n_trials)', 1, n_time);
shifted = segments(sub2ind([n_trials, n_time], rows, indices));
trace = mean(shifted, 1);
end

function print_circular_results(name, result, animal_names, behavior_names, effect_band)
fprintf('\n%s hierarchical circular-shift FOOOF result\n', name);
fprintf('  group %.1f-%.1f Hz percentile = %.3f, p = %.4g\n', ...
    effect_band(1), effect_band(2), result.group_target_percentile, result.p_effect);
fprintf('  group peak %.2f Hz, max-frequency corrected p = %.4g\n', ...
    result.group_peak_hz, result.group_peak_p);
fprintf('  unit results (periodic power ranked against that unit''s null)\n');
fprintf('  animal   behavior     actor   trials  target_pct  peak_hz  peak_p\n');
for unit_n = 1:numel(result.unit_animal)
    fprintf('  %-7s  %-11s  %-6s  %6d     %.3f    %7.2f  %.4g\n', ...
        animal_names{result.unit_animal(unit_n)}, ...
        behavior_names{result.unit_behavior(unit_n)}, ...
        result.unit_actor{unit_n}, result.unit_trials(unit_n), ...
        result.unit_target_percentile(unit_n), result.unit_peak_hz(unit_n), ...
        result.unit_peak_p(unit_n));
end
fprintf('  animal summaries (equal behavior/actor weight)\n');
fprintf('  animal   units  trials  target_pct  target_p  peak_hz  peak_p  LOO_score  LOO_p\n');
for animal_n = 1:numel(animal_names)
    if ~result.included(animal_n)
        fprintf('  %-7s %5d %7d   excluded\n', animal_names{animal_n}, ...
            result.n_units(animal_n), result.n_trials(animal_n));
        continue
    end
    fprintf('  %-7s %5d %7d      %.3f   %8.4g  %7.2f  %7.4g      %.3f  %7.4g\n', ...
        animal_names{animal_n}, result.n_units(animal_n), ...
        result.n_trials(animal_n), result.animal_target_percentile(animal_n), ...
        result.animal_p_effect(animal_n), result.animal_peak_hz(animal_n), ...
        result.animal_peak_p(animal_n), result.loo_target_score(animal_n), ...
        result.loo_p_effect(animal_n));
end
end

function plot_circular_spectrum_hierarchy(freq, play_result, nonplay_result, ...
        animal_names, effect_band)
results = {play_result, nonplay_result};
names = {'Play', 'Non-play'};
cmap = jet(256);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for column = 1:2
    result = results{column};
    n_rows = size(result.unit_percentile, 1);
    image_axis = nexttile(column, [2 1]);
    shown = imagesc(freq, 1:n_rows, result.unit_percentile);
    set(shown, 'AlphaData', isfinite(result.unit_percentile));
    set(image_axis, 'YDir', 'reverse', 'TickDir', 'out', ...
        'YLim', [0.5 n_rows + 0.5]);
    clim(image_axis, [0 1]);
    colormap(image_axis, cmap);
    animals = unique(result.unit_animal, 'stable');
    tick_at = nan(size(animals));
    hold(image_axis, 'on');
    for animal_index = 1:numel(animals)
        rows = find(result.unit_animal == animals(animal_index));
        tick_at(animal_index) = mean(rows);
        if animal_index > 1
            boundary = rows(1) - 0.5;
            plot(image_axis, [freq(1) freq(end)], boundary * [1 1], ...
                'w-', 'LineWidth', 1.2);
        end
    end
    plot(image_axis, effect_band(1) * [1 1], [0.5 n_rows + 0.5], ...
        'k:', 'LineWidth', 0.8);
    plot(image_axis, effect_band(2) * [1 1], [0.5 n_rows + 0.5], ...
        'k:', 'LineWidth', 0.8);
    hold(image_axis, 'off');
    yticks(image_axis, tick_at);
    yticklabels(image_axis, animal_names(animals));
    if column == 1
        ylabel(image_axis, 'Animal (rows are behavior x actor units)');
    end
    title(image_axis, [names{column} ' unit FOOOF percentiles']);
    xlim(image_axis, [freq(1) freq(end)]);
    colorbar(image_axis);

    line_axis = nexttile(4 + column);
    hold(line_axis, 'on');
    for animal_n = 1:numel(animal_names)
        if result.included(animal_n)
            plot(line_axis, freq, result.animal_score(animal_n, :), ...
                'Color', [0.82 0.82 0.82], 'LineWidth', 0.8, ...
                'HandleVisibility', 'off');
        end
    end
    null_lo = prctile(result.group_null_score, 2.5, 1);
    null_hi = prctile(result.group_null_score, 97.5, 1);
    null_center = median(result.group_null_score, 1, 'omitnan');
    valid = isfinite(null_lo) & isfinite(null_hi);
    if any(valid)
        x_band = freq(valid);
        fill(line_axis, [x_band fliplr(x_band)], ...
            [null_lo(valid) fliplr(null_hi(valid))], [0.25 0.45 0.85], ...
            'FaceAlpha', 0.38, 'EdgeColor', 'none', ...
            'DisplayName', 'Group null 95%');
    end
    plot(line_axis, freq, null_center, '--', 'Color', [0.15 0.35 0.75], ...
        'LineWidth', 1.1, 'DisplayName', 'Null median');
    plot(line_axis, freq, result.group_score, 'k', 'LineWidth', 2.4, ...
        'DisplayName', 'Equal-animal mean percentile');
    if isfinite(result.group_max_threshold)
        yline(line_axis, result.group_max_threshold, 'r--', ...
            'Max-frequency 5% threshold', 'HandleVisibility', 'off');
    end
    significant = result.significant_frequency & isfinite(result.group_score);
    if any(significant)
        plot(line_axis, freq(significant), result.group_score(significant), ...
            'r.', 'MarkerSize', 12, 'DisplayName', 'Max-T significant');
    end
    yl = ylim(line_axis);
    target = patch(line_axis, ...
        [effect_band(1) effect_band(2) effect_band(2) effect_band(1)], ...
        [yl(1) yl(1) yl(2) yl(2)], [0.65 0.65 0.65], ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(target, 'bottom');
    hold(line_axis, 'off');
    xlabel(line_axis, 'Frequency (Hz)');
    if column == 1
        ylabel(line_axis, 'Mean FOOOF percentile');
    end
    xlim(line_axis, [freq(1) freq(end)]);
    ylim(line_axis, [0 1]);
    set(line_axis, 'TickDir', 'out', 'Box', 'off');
    title(line_axis, sprintf(['target pct=%.3f, p=%.4g; peak %.2f Hz, ' ...
        'corrected p=%.4g'], result.group_target_percentile, ...
        result.p_effect, result.group_peak_hz, result.group_peak_p));
    legend(line_axis, 'Location', 'best');
end
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
        equal_curve, n_trials, animal_names, colors, behavior_names, ...
        pages, page_names, x_label, y_label, actor_title, shade_band)
for page_n = 1:numel(pages)
    idx = pages{page_n};
    figure
    name_psth_figure(band_name, [plot_prefix '_' page_names{page_n}])
    plot_animal_panels(x, animal_curves(idx, :, :), trial_curve(idx, :), ...
        equal_curve(idx, :), n_trials(idx, :), ...
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

function plot_animal_panels(x, animal_curves, trial_curve, equal_curve, n_trials, ...
        animal_names, colors, behavior_names, x_label, y_label, actor_title, shade_band)
% One column per behavior. The image occupies the top two rows, one row per
% animal. The bottom row is thin animal lines and thick means.
n_animals = numel(animal_names);
n_behaviors = numel(behavior_names);
if n_behaviors == 0
    return
end
is_spectrum = ~isempty(shade_band);
[c_lim, cmap] = page_color_scale(animal_curves, is_spectrum);
handles = gobjects(n_animals + 2, 1);
image_axes = gobjects(n_behaviors, 1);
line_axes = gobjects(n_behaviors, 1);
tiledlayout(3, n_behaviors, 'TileSpacing', 'compact', 'Padding', 'compact');
for behavior_n = 1:n_behaviors
    block = reshape(animal_curves(behavior_n, :, :), n_animals, []);
    image_axes(behavior_n) = nexttile(behavior_n, [2 1]);
    shown = imagesc(x, 1:n_animals, block);
    set(shown, 'AlphaData', isfinite(block));
    set(gca, 'YDir', 'reverse', 'TickDir', 'out', 'YLim', [0.5, n_animals + 0.5]);
    if ~isempty(c_lim)
        clim(c_lim);
    end
    colormap(gca, cmap);
    yticks(1:n_animals);
    if behavior_n == 1
        yticklabels(animal_names);
        ylabel('Animal');
    else
        yticklabels({});
    end
    title(behavior_names{behavior_n}, 'Interpreter', 'none');
    if numel(x) > 1
        xlim([x(1), x(end)]);
    end
    if is_spectrum && numel(shade_band) == 2
        hold on
        plot(shade_band(1) * [1 1], [0.5, n_animals + 0.5], 'w-', 'LineWidth', 0.8);
        plot(shade_band(2) * [1 1], [0.5, n_animals + 0.5], 'w-', 'LineWidth', 0.8);
        hold off
        ylim([0.5, n_animals + 0.5]);
    end

    line_axes(behavior_n) = nexttile(2 * n_behaviors + behavior_n);
    hold on
    for animal_n = 1:n_animals
        if n_trials(behavior_n, animal_n) == 0 || ~any(isfinite(block(animal_n, :)))
            continue
        end
        h = plot(x, block(animal_n, :), 'Color', colors(animal_n, :), 'LineWidth', 0.8);
        if ~isgraphics(handles(animal_n))
            set(h, 'DisplayName', animal_names{animal_n});
            handles(animal_n) = h;
        end
    end
    if any(isfinite(trial_curve(behavior_n, :)))
        h = plot(x, trial_curve(behavior_n, :), 'k', 'LineWidth', 2.2);
        if ~isgraphics(handles(n_animals + 1))
            set(h, 'DisplayName', 'All trials');
            handles(n_animals + 1) = h;
        end
    end
    if any(isfinite(equal_curve(behavior_n, :)))
        h = plot(x, equal_curve(behavior_n, :), 'Color', [0.35 0.35 0.35], 'LineWidth', 2.2);
        if ~isgraphics(handles(n_animals + 2))
            set(h, 'DisplayName', 'Mean of animals');
            handles(n_animals + 2) = h;
        end
    end
    if is_spectrum && numel(x) > 1 && numel(shade_band) == 2
        yl = ylim;
        shade = patch([shade_band(1) shade_band(2) shade_band(2) shade_band(1)], ...
            [yl(1) yl(1) yl(2) yl(2)], [0.75 0.75 0.75], ...
            'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        uistack(shade, 'bottom');
    end
    hold off
    xlabel(x_label);
    if behavior_n == 1
        ylabel(y_label);
    end
    if numel(x) > 1
        xlim([x(1), x(end)]);
    end
    set(gca, 'TickDir', 'out', 'Box', 'off');
end
if isgraphics(image_axes(end))
    scale_bar = colorbar(image_axes(end));
    scale_bar.Layout.Tile = 'east';
end
present = handles(isgraphics(handles));
if ~isempty(present) && isgraphics(line_axes(1))
    line_legend = legend(line_axes(1), present, 'Orientation', 'horizontal');
    line_legend.Layout.Tile = 'south';
end
sgtitle(actor_title);
end

function [c_lim, cmap] = page_color_scale(animal_curves, is_spectrum)
vals = animal_curves(isfinite(animal_curves));
c_lim = [];
if is_spectrum
    cmap = jet(256);
    if numel(vals) >= 8
        lo = prctile(vals, 5);
        hi = prctile(vals, 95);
        if hi > lo
            c_lim = [lo hi];
        end
    end
else
    cmap = diverging_colormap(256);
    if ~isempty(vals)
        lim = prctile(abs(vals), 98);
        if lim > 0
            c_lim = [-lim lim];
        end
    end
end
end

function cmap = diverging_colormap(n)
n_half = floor(n / 2);
blue = [linspace(0.12, 1, n_half)', linspace(0.25, 1, n_half)', ones(n_half, 1)];
red = [ones(n_half, 1), linspace(1, 0.15, n_half)', linspace(1, 0.12, n_half)'];
cmap = [blue; red];
end

function [curves, animal_ids] = stack_actor_behaviors(curves_self, curves_other, ...
        n_self, n_other, behavior_index)
n_animals = size(n_self, 2);
n_x = size(curves_self, 3);
curves = zeros(0, n_x);
animal_ids = zeros(0, 1);
for animal_n = 1:n_animals
    for behavior_n = behavior_index
        if n_self(behavior_n, animal_n) > 0
            row = reshape(curves_self(behavior_n, animal_n, :), 1, []);
            if any(isfinite(row))
                curves(end + 1, :) = row; %#ok<AGROW>
                animal_ids(end + 1, 1) = animal_n; %#ok<AGROW>
            end
        end
        if n_other(behavior_n, animal_n) > 0
            row = reshape(curves_other(behavior_n, animal_n, :), 1, []);
            if any(isfinite(row))
                curves(end + 1, :) = row; %#ok<AGROW>
                animal_ids(end + 1, 1) = animal_n; %#ok<AGROW>
            end
        end
    end
end
end

function [mu, lo, hi] = repeated_animal_mean(curves, animal_ids)
% Average behaviors within each animal, then take the mean and the 95% t
% interval across animals. Extra behaviors from the same animal do not add
% degrees of freedom.
n_x = size(curves, 2);
mu = nan(1, n_x);
lo = nan(1, n_x);
hi = nan(1, n_x);
if isempty(curves)
    return
end
for t = 1:n_x
    y = curves(:, t);
    valid = isfinite(y);
    if ~any(valid)
        continue
    end
    [estimate, se, df] = animal_level_mean(y(valid), animal_ids(valid));
    mu(t) = estimate;
    if isfinite(se) && df > 0
        crit = tinv(0.975, df);
        lo(t) = estimate - crit * se;
        hi(t) = estimate + crit * se;
    end
end
end

function [estimate, se, df] = animal_level_mean(y, ids)
y = y(:);
ids = ids(:);
animals = unique(ids, 'stable');
a = numel(animals);
m = zeros(a, 1);
for animal_n = 1:a
    m(animal_n) = mean(y(ids == animals(animal_n)));
end
estimate = mean(m);
df = a - 1;
if a < 2
    se = NaN;
    return
end
se = std(m, 0) / sqrt(a);
end

function plot_pooled_repeated(x, curves, animal_ids, animal_names, mu, lo, hi, ...
        x_label, y_label, figure_title, shade_band)
if isempty(curves)
    axis off
    title(figure_title, 'Interpreter', 'none')
    return
end
is_spectrum = ~isempty(shade_band);
[c_lim, cmap] = page_color_scale(curves, is_spectrum);
n_rows = size(curves, 1);
animals = unique(animal_ids, 'stable');
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(1, [2 1])
shown = imagesc(x, 1:n_rows, curves);
set(shown, 'AlphaData', isfinite(curves));
set(gca, 'YDir', 'reverse', 'TickDir', 'out', 'YLim', [0.5, n_rows + 0.5]);
if ~isempty(c_lim)
    clim(c_lim);
end
colormap(gca, cmap);
hold on
tick_at = zeros(1, numel(animals));
for animal_n = 1:numel(animals)
    rows = find(animal_ids == animals(animal_n));
    tick_at(animal_n) = mean(rows);
    if animal_n > 1
        boundary = rows(1) - 0.5;
        plot([x(1) x(end)], boundary * [1 1], 'w-', 'LineWidth', 1.2);
    end
end
hold off
yticks(tick_at);
yticklabels(animal_names(animals));
ylabel('Animal');
title(figure_title, 'Interpreter', 'none');
if numel(x) > 1
    xlim([x(1), x(end)]);
end
if is_spectrum && numel(shade_band) == 2
    hold on
    plot(shade_band(1) * [1 1], [0.5, n_rows + 0.5], 'w-', 'LineWidth', 0.8);
    plot(shade_band(2) * [1 1], [0.5, n_rows + 0.5], 'w-', 'LineWidth', 0.8);
    hold off
    ylim([0.5, n_rows + 0.5]);
end
colorbar

nexttile(3)
hold on
for row = 1:n_rows
    plot(x, curves(row, :), 'Color', [0.82 0.82 0.82], 'LineWidth', 0.6, ...
        'HandleVisibility', 'off');
end
finite_band = isfinite(lo) & isfinite(hi);
if any(finite_band)
    x_band = x(finite_band);
    fill([x_band, fliplr(x_band)], [lo(finite_band), fliplr(hi(finite_band))], ...
        [0.20 0.35 0.70], 'FaceAlpha', 0.45, 'EdgeColor', 'none', ...
        'DisplayName', '95% CI');
end
plot(x, mu, 'k', 'LineWidth', 2.4, 'DisplayName', 'Grand mean');
if is_spectrum && numel(x) > 1 && numel(shade_band) == 2
    yl = ylim;
    shade = patch([shade_band(1) shade_band(2) shade_band(2) shade_band(1)], ...
        [yl(1) yl(1) yl(2) yl(2)], [0.75 0.75 0.75], ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(shade, 'bottom');
end
hold off
xlabel(x_label);
ylabel(y_label);
if numel(x) > 1
    xlim([x(1), x(end)]);
end
set(gca, 'TickDir', 'out', 'Box', 'off');
legend('Location', 'best');
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

function periodic = periodic_residual_fooof(freq, linear_power, aperiodic_mode, peak_widths)
% Observed log10 spectrum minus FOOOF's final aperiodic fit. This is the
% inferential periodic residual requested for the circular null.
periodic = nan(1, numel(freq));
freq = freq(:).';
linear_power = linear_power(:).';
if numel(freq) < 8 || numel(linear_power) ~= numel(freq)
    return
end
if any(~isfinite(freq)) || any(~isfinite(linear_power)) || all(linear_power <= 0)
    return
end
linear_power = max(linear_power, realmin('double'));
out = fooof_matlab(freq(:), linear_power(:), ...
    'AperiodicMode', aperiodic_mode, ...
    'PeakWidthLimits', peak_widths);
periodic = (out.log_power(:) - out.aperiodic_fit(:)).';
end
