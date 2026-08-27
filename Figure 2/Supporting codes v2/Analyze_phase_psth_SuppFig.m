%% 1 Analyze_phase_psth
% Phase reset analysis at bout onset (circular stats + surrogate test).
%
% Precomputed file combinations (frequency band, behavior, calls, CV) are listed in:
%   Figure 2/Figure 2 Psth animal names and result combinations.txt
%
% This script uses the PHASE row. For band-limited power PSTHs (delta/theta/gamma
% play bout, calls, etc.), swap the load filenames per that list — see Figure 2.m.

this_file = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Analyze_phase_psth_SuppFig.m';
repo_root = fileparts(fileparts(fileparts(this_file)));
figure2_dir = fullfile(repo_root, 'Figure 2');
combo_file = fullfile(figure2_dir, 'Figure 2 Psth animal names and result combinations.txt');
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');
figure_2_folder = fullfile(figure2_dir, 'outputs');

%% 2 Load precomputed phase-onset data
% Default: PHASE section in combo file. Other examples (commented):
%   Delta play bout  → psth_structure_delta_updated.mat / animal_names_delta_updated.mat
%   Delta exploratory → psth_structure_delta_exploratory_bout.mat / animal_names_delta_exploratory_bout.mat
%   Theta play bout  → psth_structure_theta_updated.mat / animal_names_theta_updated.mat
%   Gamma calls      → psth_structure_call_gamma.mat / animal_names_call_gamma.mat

load(fullfile(saving_folder, 'phase_onset_playbout.mat'), 'psth_structure');
load(fullfile(saving_folder, 'animal_names_dphase_onset_playbout.mat'), 'animal_names');

hist_range = psth_structure(1).hist_range;
time = linspace(hist_range(1), hist_range(2), size(psth_structure(1).Phase_data_onset, 2));

%% 3 Concatenate phase across sessions
all_phase_onset = [];
all_phase_offset = [];

for j = [1:11, 13]
    all_phase_onset = [all_phase_onset; psth_structure(j).Phase_data_onset];
    all_phase_offset = [all_phase_offset; psth_structure(j).Phase_data_offset];
end

%% 4 Mean phase and resultant length over time
figure
subplot(2, 1, 1)
plot(time, circ_mean(all_phase_onset))
xlim([-1 2])

subplot(2, 1, 2)
plot(time, circ_r(all_phase_onset))
xlim([-1 2])

%% 5 Permutation test against circular surrogate
time_range = [-5 5];
time4perm = time(time >= time_range(1) & time <= time_range(2));
all_phase_onset4perm = all_phase_onset(:, time >= time_range(1) & time <= time_range(2));

load(fullfile(saving_folder, 'phase_surrogate.mat'), 'surrogate_r');
n_perm = size(surrogate_r, 1);
pctl = sum(surrogate_r > circ_r(all_phase_onset4perm)) / (n_perm + 1);

figure
plot(time4perm, 1 - pctl)
hold on
yyaxis right
plot(time4perm, circ_r(all_phase_onset4perm))
xlim([-2 4])

%% 6 Peak/trough phase with significance-colored resultant length
mean_angle = circ_mean(all_phase_onset4perm);
trough_angles = mean_angle;
trough_angles(trough_angles < pi/2 & trough_angles > -pi/2) = NaN;
peak_angles = mean_angle;
peak_angles(~(peak_angles < pi/2 & peak_angles > -pi/2)) = NaN;

figure
plot(time4perm, 1 - pctl, 'k-')
hold on
yyaxis right
plot(time4perm, trough_angles, '-b')
plot(time4perm, peak_angles, 'r-')
xlim([-2 4])

%% 7 Detailed phase-reset figure with zoom panel
x_lim = [-.5 1];
mean_angle = circ_mean(all_phase_onset4perm);
jumps = abs(diff(mean_angle)) > pi;
time_broken = time4perm;
angle_broken = mean_angle;
time_broken([false, jumps]) = NaN;
angle_broken([false, jumps]) = NaN;

trough_angles = angle_broken;
trough_angles(trough_angles < pi/2 & trough_angles > -pi/2) = NaN;
peak_angles = angle_broken;
peak_angles(~(peak_angles < pi/2 & peak_angles > -pi/2)) = NaN;

blue_colored_pctl = 1 - pctl;
blue_colored_pctl(mean_angle < pi/2 & mean_angle > -pi/2) = NaN;
red_colored_pctl = 1 - pctl;
red_colored_pctl(~(mean_angle < pi/2 & mean_angle > -pi/2)) = NaN;

figure('units', 'normalized', 'outerposition', [0 0 .5 1]);
subplot(3, 1, 1:2)
plot(time4perm, blue_colored_pctl, 'b-', 'LineWidth', 2)
hold on
plot(time4perm, red_colored_pctl, 'r-', 'LineWidth', 2)
ylim([0 1])
yyaxis right
plot(time_broken, trough_angles, ':b')
plot(time_broken, trough_angles + 2*pi, ':b')
plot(time_broken, trough_angles - 2*pi, ':b')
plot(time_broken, peak_angles, ':r')
plot(time_broken, peak_angles + 2*pi, ':r')
plot(time_broken, peak_angles - 2*pi, ':r')
ylabel('Phase (rad)')
ylim([-pi pi])
xlim(x_lim)
set(gca, 'TickDir', 'out')

subplot(3, 1, 3)
plot(time4perm, blue_colored_pctl, 'b-', 'LineWidth', 2)
hold on
plot(time4perm, red_colored_pctl, 'r-', 'LineWidth', 2)
ylim([.95 1])
yyaxis right
plot(time_broken, trough_angles, ':b')
plot(time_broken, trough_angles + 2*pi, ':b')
plot(time_broken, trough_angles - 2*pi, ':b')
plot(time_broken, peak_angles, ':r')
plot(time_broken, peak_angles + 2*pi, ':r')
plot(time_broken, peak_angles - 2*pi, ':r')
ylabel('Phase (rad)')
ylim([-pi pi])
xticks([-.5 0 .312 .5 .66 1])
xlim(x_lim)
set(gca, 'TickDir', 'out')

if ~exist(figure_2_folder, 'dir'), mkdir(figure_2_folder); end
print(gcf, '-vector', '-dsvg', fullfile(figure_2_folder, 'phasse reset during play.svg'))

%% 8 (AS FORMER SUPP FIGURE) Multi-panel phase reset with polar histograms at key times
pctl = sum(surrogate_r > circ_r(all_phase_onset4perm)) / (n_perm + 1);
r_lim = [0 .28];
n_bins = 36;
y_lim = [0 .12];
angle_values = [-.5 0 .42 .5 1];

mean_phases = circ_mean(all_phase_onset4perm);
trough_angles = ~(mean_phases < pi/2 & mean_phases > -pi/2);
peak_angles = ~trough_angles;

figure
subplot(5, 1, 1)
plot(time4perm, circ_mean(all_phase_onset4perm), 'k')

subplot(5, 1, 2)
hold on
ylim(y_lim)
start_end = find_beg_end(trough_angles);
for j = 1:size(start_end, 1)
    fill(time4perm(start_end(j, [1 2 2 1])), y_lim([1 1 2 2]), 'b', 'FaceAlpha', .2, 'EdgeColor', 'none')
end
start_end = find_beg_end(peak_angles);
for j = 1:size(start_end, 1)
    fill(time4perm(start_end(j, [1 2 2 1])), y_lim([1 1 2 2]), 'r', 'FaceAlpha', .2, 'EdgeColor', 'none')
end
plot(time4perm, circ_r(all_phase_onset4perm), 'k')
significnat_r = circ_r(all_phase_onset4perm);
significnat_r(pctl > 0.05) = NaN;
plot(time4perm, significnat_r, 'k', 'LineWidth', 3)
plot([0 0], y_lim, ':k')
xlim([-1 2])

for j = 1:numel(angle_values)
    subplot(5, numel(angle_values), j + [2 3 4] * numel(angle_values))
    these_angles = all_phase_onset4perm(:, round(time4perm, 3) == angle_values(j));
    these_angles = these_angles(:);
    median_angle = circ_mean(these_angles);
    polarhistogram(these_angles, linspace(-pi, pi, n_bins), 'FaceColor', 'k', 'FaceAlpha', .2, 'EdgeColor', 'none', 'Normalization', 'pdf')
    [pdf_estimate] = circ_ksdensity(these_angles, linspace(-pi, pi, n_bins));
    pdf_estimate = movmean([pdf_estimate; pdf_estimate; pdf_estimate], 5);
    pdf_estimate = pdf_estimate((n_bins + 1):(2 * n_bins));
    hold on
    polarplot(linspace(-pi, pi, n_bins), 5 * pdf_estimate, 'k', 'LineWidth', 2)
    rlim(r_lim);
    r = circ_r(these_angles);
    polarplot([median_angle median_angle], r_lim * r / 0.1, 'r', 'LineWidth', 2)
    title([num2str(r), ' at ', num2str(angle_values(j))])
end

print(gcf, '-vector', '-dsvg', fullfile(figure_2_folder, 'phasse reset during play witn alge description.svg'))
