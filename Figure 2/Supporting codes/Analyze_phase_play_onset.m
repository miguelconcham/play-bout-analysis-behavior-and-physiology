

npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
figure_2_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure 2 Inputs';

animal_file_names =  cellfun(@(x) ['B', x],strsplit([animal_list.name], 'B'), 'UniformOutput',false)';
animal_file_names(1) = [];
% animal2exclude = {'B4D4 0826 Dual'};

freq_range_1    = [1 5];
sr              = 2500;
filter_order    = 2000;



% Parameters for delta
Hd_freq = designfilt('bandpassfir', ...
'FilterOrder', filter_order, ...
'CutoffFrequency1', freq_range_1(1), ...
'CutoffFrequency2', freq_range_1(2), ...
'SampleRate', sr, ...
'DesignMethod', 'window', ...
'Window', 'hamming');
bin_size_freq = 0.01;

% Parameters for theta
% Hd_freq = designfilt('bandpassfir', ...
% 'FilterOrder', filter_order, ...
% 'CutoffFrequency1', freq_range_2(1), ...
% 'CutoffFrequency2', freq_range_2(2), ...
% 'SampleRate', sr, ...
% 'DesignMethod', 'window', ...
% 'Window', 'hamming');
% bin_size_freq = 0.001;



%%
% 
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';


load([saving_folder,'\phase_onset_playbout.mat'],'psth_structure');
load([saving_folder,'\animal_names_dphase_onset_playbout.mat'],'animal_names');
%%

all_phase_onset = [];
all_phase_offset = [];


for j=[1:11,13]
    all_phase_onset = [all_phase_onset;psth_structure(j).Phase_data_onset];
    all_phase_offset = [all_phase_offset;psth_structure(j).Phase_data_offset];
end



%%


figure


subplot(2,1,1)
plot(time,circ_mean(all_phase_onset))
xlim([-1 2])

subplot(2,1,2)
plot(time,circ_r(all_phase_onset))
xlim([-1 2])
%%
time_range = [-5 5];
time4perm = time(time>= time_range(1) & time<= time_range(2));
all_phase_onset4perm = all_phase_onset(:,time>= time_range(1) & time<= time_range(2));
%% generate=surroage r
% n_perm = 10000;
% [num_rows, num_cols] = size(all_phase_onset4perm);
% surrogate_r = nan(n_perm, num_cols);
% 
% % Pre-calculate column indices for the whole matrix to avoid re-creating them
% col_idx = 1:num_cols;
%%
for p = 1:n_perm
    % 1. Generate one random shift for every row at once
    shifts = randi([0, num_cols-1], num_rows, 1);
    
    % 2. Create a matrix of shifted column indices
    % Use bsxfun or broadcasting to shift the col_idx by the random amounts
    shifted_indices = mod(bsxfun(@plus, col_idx-1, shifts), num_cols) + 1;
    
    % 3. Use linear indexing to "shuffle" the rows in one shot
    % Convert row and shifted column indices into linear indices
    row_idx = repmat((1:num_rows)', 1, num_cols);
    lin_idx = sub2ind([num_rows, num_cols], row_idx, shifted_indices);
    
    shuffled_data = all_phase_onset4perm(lin_idx);
    
    % 4. Calculate R for all columns at once
    surrogate_r(p, :) = circ_r(shuffled_data);
    
    if mod(p, 500) == 0
        fprintf('Permutation %d/%d completed...\n', p, n_perm);
    end
end
save([saving_folder,'\phase_surrogate.mat'],'surrogate_r', '-v7.3');

%%
load([saving_folder,'\phase_surrogate.mat'],'surrogate_r');

sub_surroage_r = surrogate_r;

pctl = sum(sub_surroage_r>circ_r(all_phase_onset4perm))/(n_perm+1);


figure_2_folder
plot(time4perm,1-pctl)
hold on
yyaxis right
plot(time4perm,circ_r(all_phase_onset4perm))
xlim([-2 4])

%%

mean_angle = circ_mean(all_phase_onset4perm);

trough_angles = mean_angle;
trough_angles(trough_angles<pi/2 & trough_angles>-pi/2) = NaN;
 peak_angles = mean_angle;
peak_angles(~(peak_angles<pi/2 & peak_angles>-pi/2)) = NaN;
 


plot(time4perm,1-pctl, 'k-')
hold on
yyaxis right
plot(time4perm, trough_angles, '-b')
plot(time4perm, peak_angles, 'r-')

xlim([-2 4])




%% 
x_lim = [-.5 1];
% 1. Calculate the mean angle
mean_angle = circ_mean(all_phase_onset4perm);

% 2. Detect jumps (phase wrapping) and insert NaNs
% We find where the absolute difference between points is > pi
jumps = abs(diff(mean_angle)) > pi;

% Create a version of the time and angle vectors with NaNs at jump points
% This 'breaks' the line in the plot
time_broken = time4perm;
angle_broken = mean_angle;

% Insert NaNs where the jump occurs
time_broken([false, jumps]) = NaN;
angle_broken([false, jumps]) = NaN;

% 3. Split into peaks and troughs as before
trough_angles = angle_broken;
trough_angles(trough_angles < pi/2 & trough_angles > -pi/2) = NaN;

peak_angles = angle_broken;
peak_angles(~(peak_angles < pi/2 & peak_angles > -pi/2)) = NaN;

blue_colored_pctl = 1-pctl;
blue_colored_pctl(mean_angle < pi/2 & mean_angle > -pi/2) = NaN;

red_colored_pctl = 1-pctl;
red_colored_pctl(~(mean_angle < pi/2 & mean_angle > -pi/2)) = NaN;

% 4. Plot
figure_2_folder('units','normalized','outerposition',[0 0 .5 1]);
subplot(3,1,1:2)
% plot(time4perm, 1-pctl, 'k-')
plot(time4perm, blue_colored_pctl, 'b-', 'LineWidth',2)
hold on

plot(time4perm,red_colored_pctl, 'r-', 'LineWidth',2)
ylim([0 1])
yyaxis right
plot(time_broken, trough_angles, ':b') % Use time_broken to respect the NaNs
plot(time_broken, trough_angles+2*pi, ':b') % Use time_broken to respect the NaNs
plot(time_broken, trough_angles-2*pi, ':b') 
plot(time_broken, peak_angles, ':r')
plot(time_broken, peak_angles+2*pi, ':r')
plot(time_broken, peak_angles -2*pi, ':r')
ylabel('Phase (rad)')
ylim([-pi pi])

xlim(x_lim)
set(gca, 'TickDir', 'out'); % Set tick direction to out


subplot(3,1,3)
% plot(time4perm, 1-pctl, 'k-')
plot(time4perm, blue_colored_pctl, 'b-', 'LineWidth',2)
hold on

plot(time4perm,red_colored_pctl, 'r-', 'LineWidth',2)
ylim([.95 1])
yyaxis right
plot(time_broken, trough_angles, ':b') % Use time_broken to respect the NaNs
plot(time_broken, trough_angles+2*pi, ':b') % Use time_broken to respect the NaNs
plot(time_broken, trough_angles-2*pi, ':b') 
plot(time_broken, peak_angles, ':r')
plot(time_broken, peak_angles+2*pi, ':r')
plot(time_broken, peak_angles -2*pi, ':r')
ylabel('Phase (rad)')
ylim([-pi pi])
xticks([-.5 0 .312 .5 .66 1 ])

xlim(x_lim)
set(gca, 'TickDir', 'out'); % Set tick direction to out

%% save if needed

print(gcf,'-vector','-dsvg',[figure_2_folder,'/phasse reset during play.svg'])

%%

%% ploting all together

n_perm = 10000;
pctl = sum(surrogate_r>circ_r(all_phase_onset4perm))/(n_perm+1);


r_lim = [0 .28];
figure
sp_n = 1;
n_bins = 36;
y_lim = [0 .12];




mean_phases     = circ_mean(all_phase_onset4perm);
trough_angles   = ~(mean_phases < pi/2 & mean_phases > -pi/2);
peak_angles     = ~trough_angles;


subplot(5,1,1)
plot(time4perm,circ_mean(all_phase_onset4perm), 'k')


subplot(5,1,2)
hold on
ylim([0 .12])
start_end = find_beg_end(trough_angles);
for j=1:size(start_end,1)
    fill(time4perm(start_end(j,[1 2 2 1])), y_lim([1 1 2 2]), 'b', 'FaceAlpha',.2, 'EdgeColor','none')
end

start_end = find_beg_end(peak_angles);
for j=1:size(start_end,1)
    fill(time4perm(start_end(j,[1 2 2 1])), y_lim([1 1 2 2]), 'r', 'FaceAlpha',.2, 'EdgeColor','none')
end  
plot(time4perm,circ_r(all_phase_onset4perm), 'k')

significnat_r = circ_r(all_phase_onset4perm);
significnat_r(pctl>0.05) = NaN;

plot(time4perm,significnat_r, 'k', 'LineWidth',3)
plot([0 0], y_lim, ':k')



xlim([-1 2])
angle_values = [ -1 -.5 0 .42 .5 1 2 ]
angle_values = [-.5 0 .42 .5 1 ]
for j=1:numel(angle_values)
    subplot(5,numel(angle_values),j +[2 3 4]*numel(angle_values)  )
    these_angles = all_phase_onset4perm(:,round(time4perm,3)==angle_values(j) );
    these_angles = these_angles(:);
    median_angle = circ_mean(these_angles);
    polarhistogram(these_angles, linspace(-pi, pi, n_bins), 'FaceColor','k', 'FaceAlpha',.2, 'EdgeColor', 'none', 'Normalization','pdf')
    [counts, edges] = histcounts(these_angles, 'BinLimits', [-pi, pi], 'NumBins', 36);
    [pdf_estimate] = circ_ksdensity(these_angles,  linspace(-pi, pi, n_bins));
    pdf_estimate = [pdf_estimate; pdf_estimate; pdf_estimate];
    pdf_estimate = movmean(pdf_estimate,5);
    pdf_estimate = pdf_estimate((n_bins+1):(2*n_bins));

    hold on    
    polarplot(linspace(-pi, pi, n_bins),5*pdf_estimate, 'k', 'LineWidth',2)
    rlim(r_lim);
    r = circ_r(these_angles);
    hold on
    polarplot([median_angle median_angle], r_lim*r/0.1, 'r', 'LineWidth',2)
    title([num2str(r), ' at ', num2str(angle_values(j))])
    sp_n = sp_n+1;
end

%%

print(gcf,'-vector','-dsvg',[figure_2_folder,'/phasse reset during play witn alge description.svg'])
