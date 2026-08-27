clear all

saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\Theta psth';

% % %PSTH TO THE TROUGH
% load([saving_folder,'\delta_phase_couplig_structure_updated_with_non_playbouts_trough_peaks.mat']);
% load([saving_folder,'\delta_phase_couplig_animal_names_updated_with_non_playbouts_trough_peaks.mat']);
% %PSTH TO THE PEAK
load([saving_folder,'\delta_phase_couplig_structure_updated_with_non_playbouts.mat'],'phase_struct');

load([saving_folder,'\delta_phase_couplig_animal_names_updated_with_non_playbouts.mat'],'animal_names');

phase_prop_names = {'PreferedAngle','MVL','MVLPval','PPC','PPCPval','MeanRate','Id'};

all_neurons      = [];
all_play_psth    = [];
all_preplay_psth = [];

alpha = 0.05;

for ns = 1:numel(phase_struct)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PSTH (SESSION-AVERAGED)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    play_psth_this = squeeze(phase_struct(ns).play_psth);
    preplay_psth_this = squeeze(phase_struct(ns).pre_play_psth);

    play_psth_mean    = squeeze(nanmean(play_psth_this, 1));
    preplay_psth_mean = squeeze(nanmean(preplay_psth_this, 1));

    all_play_psth    = cat(1, all_play_psth, play_psth_mean);
    all_preplay_psth = cat(1, all_preplay_psth, preplay_psth_mean);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CLUSTER INFO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    this_session_cluster_info = phase_struct(ns).cluster_info;

    % Partner 1 / 2
    sub_Table = array2table(squeeze(phase_struct(ns).session_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.Partner1 = sub_Table;

    sub_Table = array2table(squeeze(phase_struct(ns).session_phase_stats(2,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.Partner2 = sub_Table;

    % Play / PrePlay / NonPlay / PreNonPlay
    sub_Table = array2table(squeeze(phase_struct(ns).play_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.Play = sub_Table;

    sub_Table = array2table(squeeze(phase_struct(ns).pre_play_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.PrePlay = sub_Table;

    sub_Table = array2table(squeeze(phase_struct(ns).non_play_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.NonPlay = sub_Table;

    sub_Table = array2table(squeeze(phase_struct(ns).pre_non_play_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.PreNonPlay = sub_Table;

    % Entire session
    sub_Table = array2table(squeeze(phase_struct(ns).entire_recording_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.EntireSession = sub_Table;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TRAINED CELL FLAG (based on ENTIRE SESSION PPC p-value)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    entire_stats = squeeze(phase_struct(ns).entire_recording_phase_stats(1,:,:));
    PPCp_entire  = entire_stats(:,5);

    this_session_cluster_info.IsTrainedCell = PPCp_entire < alpha;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SESSION LABEL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    this_session_cluster_info.session = repmat(animal_names(ns,1), ...
        size(this_session_cluster_info,1), 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STACK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_neurons = [all_neurons; this_session_cluster_info];

end

%%
  psth_time = phase_struct(1).edges_freq;
  psth_time = (psth_time(1:end-1)+psth_time(2:end))/2;

trough_neurons = all_neurons.EntireSession.PPCPval<0.01 & (all_neurons.EntireSession.PreferedAngle>pi/2 | all_neurons.EntireSession.PreferedAngle<-pi/2);


figure

plot(psth_time,mean(all_play_psth(trough_neurons,:) - mean(mean(all_play_psth(trough_neurons,:)),2), 'omitmissing'), 'b')
hold on
plot(psth_time,mean(all_preplay_psth(trough_neurons,:) - mean(mean(all_preplay_psth(trough_neurons,:)),2), 'omitmissing'), 'k')
%% trough nd peak oscillation component
y_lim = [-0.015 0.015];
data_play    = all_play_psth(trough_neurons,:) ...
    - mean(all_play_psth(trough_neurons,:),2,'omitmissing');

data_preplay = all_preplay_psth(trough_neurons,:) ...
    - mean(all_preplay_psth(trough_neurons,:),2,'omitmissing');

% mean across neurons
m_play    = mean(data_play,1,'omitmissing');
m_preplay = mean(data_preplay,1,'omitmissing');

% SEM across neurons
sem_play    = std(data_play,0,1,'omitmissing') ./ sqrt(sum(~isnan(data_play(:,1))));
sem_preplay = std(data_preplay,0,1,'omitmissing') ./ sqrt(sum(~isnan(data_preplay(:,1))));

% CI (95% or 99%)
ci_play    = 2.58 * sem_play; %1.96 = 95% 2.58 = 99%
ci_preplay = 2.58 * sem_preplay;

figure
subplot(1,2,1)

% PLAY
plot(psth_time, m_play, 'b', 'LineWidth', 2)
hold on
fill([psth_time fliplr(psth_time)], ...
     [m_play-ci_play fliplr(m_play+ci_play)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

% PRE-PLAY
plot(psth_time, m_preplay, 'k', 'LineWidth', 2)
fill([psth_time fliplr(psth_time)], ...
     [m_preplay-ci_preplay fliplr(m_preplay+ci_preplay)], ...
     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot([0 0],y_lim, 'r')
ylim(y_lim)
% peak oscillation component

peak_neurons = all_neurons.EntireSession.PPCPval<0.01 & ~(all_neurons.EntireSession.PreferedAngle>pi/2 | all_neurons.EntireSession.PreferedAngle<-pi/2);



data_play    = all_play_psth(peak_neurons,:) ...
    - mean(all_play_psth(peak_neurons,:),2,'omitmissing');

data_preplay = all_preplay_psth(peak_neurons,:) ...
    - mean(all_preplay_psth(peak_neurons,:),2,'omitmissing');

% mean across neurons
m_play    = mean(data_play,1,'omitmissing');
m_preplay = mean(data_preplay,1,'omitmissing');

% SEM across neurons
sem_play    = std(data_play,0,1,'omitmissing') ./ sqrt(sum(~isnan(data_play(:,1))));
sem_preplay = std(data_preplay,0,1,'omitmissing') ./ sqrt(sum(~isnan(data_preplay(:,1))));

% CI (95% or 99%)
ci_play    = 2.58 * sem_play;
ci_preplay = 2.58 * sem_preplay;


subplot(1,2,2)

% PLAY
plot(psth_time, m_play, 'b', 'LineWidth', 2)
hold on
fill([psth_time fliplr(psth_time)], ...
     [m_play-ci_play fliplr(m_play+ci_play)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

% PRE-PLAY
plot(psth_time, m_preplay, 'k', 'LineWidth', 2)
fill([psth_time fliplr(psth_time)], ...
     [m_preplay-ci_preplay fliplr(m_preplay+ci_preplay)], ...
     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

plot([0 0],y_lim, 'r')
ylim(y_lim)
%% trough and peak NON oscillation component

data_play    = all_play_psth(trough_neurons,:);

data_preplay = all_preplay_psth(trough_neurons,:);

% mean across neurons
m_play    = mean(data_play,1,'omitmissing');
m_preplay = mean(data_preplay,1,'omitmissing');

% SEM across neurons
sem_play    = std(data_play,0,1,'omitmissing') ./ sqrt(sum(~isnan(data_play(:,1))));
sem_preplay = std(data_preplay,0,1,'omitmissing') ./ sqrt(sum(~isnan(data_preplay(:,1))));

% CI (95%)
ci_play    = 1.96 * sem_play;
ci_preplay = 1.96 * sem_preplay;

figure
subplot(1,2,1)

% PLAY
plot(psth_time, m_play, 'b', 'LineWidth', 2)
hold on
fill([psth_time fliplr(psth_time)], ...
     [m_play-ci_play fliplr(m_play+ci_play)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

% PRE-PLAY
plot(psth_time, m_preplay, 'k', 'LineWidth', 2)
fill([psth_time fliplr(psth_time)], ...
     [m_preplay-ci_preplay fliplr(m_preplay+ci_preplay)], ...
     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
% ylim([-.02 .02])
% peak oscillation component

peak_neurons = all_neurons.EntireSession.PPCPval<0.01 & ~(all_neurons.EntireSession.PreferedAngle>pi/2 | all_neurons.EntireSession.PreferedAngle<-pi/2);



data_play    = all_play_psth(peak_neurons,:) ;

data_preplay = all_preplay_psth(peak_neurons,:);

% mean across neurons
m_play    = mean(data_play,1,'omitmissing');
m_preplay = mean(data_preplay,1,'omitmissing');

% SEM across neurons
sem_play    = std(data_play,0,1,'omitmissing') ./ sqrt(sum(~isnan(data_play(:,1))));
sem_preplay = std(data_preplay,0,1,'omitmissing') ./ sqrt(sum(~isnan(data_preplay(:,1))));

% CI (95%)
ci_play    = 1.96 * sem_play;
ci_preplay = 1.96 * sem_preplay;


subplot(1,2,2)

% PLAY
plot(psth_time, m_play, 'b', 'LineWidth', 2)
hold on
fill([psth_time fliplr(psth_time)], ...
     [m_play-ci_play fliplr(m_play+ci_play)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

% PRE-PLAY
plot(psth_time, m_preplay, 'k', 'LineWidth', 2)
fill([psth_time fliplr(psth_time)], ...
     [m_preplay-ci_preplay fliplr(m_preplay+ci_preplay)], ...
     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

% ylim([-.02 .02])

