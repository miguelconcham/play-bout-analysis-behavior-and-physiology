%% Analyze_phase_exploratory_onset
% Delta phase reset at exploratory-behavior bout onset.

saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Theta psth';
figure_2_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Codes\Figure codes\Figure 2 Inputs';

%% Load precomputed phase-onset data
load([saving_folder,'\phase_onset_playbout.mat'],'psth_structure');
load([saving_folder,'\animal_names_dphase_onset_playbout.mat'],'animal_names');

hist_range = psth_structure(1).hist_range;
lfp_sr = 2500;
time = linspace(hist_range(1), hist_range(2), size(psth_structure(1).Phase_data_onset, 2));

%% Concatenate phase across sessions
all_phase_onset = [];
all_phase_offset = [];

for j = [1:11, 13]
    all_phase_onset = [all_phase_onset; psth_structure(j).Phase_data_onset];
    all_phase_offset = [all_phase_offset; psth_structure(j).Phase_data_offset];
end

%% Mean phase and resultant length over time
figure
subplot(2, 1, 1)
plot(time, circ_mean(all_phase_onset))
xlim([-1 2])

subplot(2, 1, 2)
plot(time, circ_r(all_phase_onset))
xlim([-1 2])

%% Permutation test against circular surrogate
time_range = [-5 5];
time4perm = time(time >= time_range(1) & time <= time_range(2));
all_phase_onset4perm = all_phase_onset(:, time >= time_range(1) & time <= time_range(2));

load([saving_folder,'\phase_surrogate.mat'],'surrogate_r');
n_perm = size(surrogate_r, 1);
pctl = sum(surrogate_r > circ_r(all_phase_onset4perm)) / (n_perm + 1);

figure
plot(time4perm, 1 - pctl)
hold on
yyaxis right
plot(time4perm, circ_r(all_phase_onset4perm))
xlim([-2 4])

%% Peak/trough phase with significance-colored resultant length
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

%% Detailed phase-reset figure with zoom panel
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

print(gcf, '-vector', '-dsvg', [figure_2_folder, '/phasse reset during play.svg'])

%% Multi-panel phase reset with polar histograms at key times
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

print(gcf, '-vector', '-dsvg', [figure_2_folder, '/phasse reset during play witn alge description.svg'])
