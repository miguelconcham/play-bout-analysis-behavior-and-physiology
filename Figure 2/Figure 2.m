%% Figure 2 — main plotting script
% PSTH / results file combinations (frequency band, behavior, calls, CV):
%   Figure 2 Psth animal names and result combinations.txt

repo_root = fileparts(mfilename('fullpath'));
combo_file = fullfile(repo_root, 'Figure 2 Psth animal names and result combinations.txt');
data_root = fullfile(repo_root, '..', 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');

figure_dir = fullfile(repo_root, 'outputs');
synch_directory = fullfile(data_root, 'Synch data');
hmm_raw_data = fullfile(data_root, 'HMM data', 'HMM raw data');
call_folder = fullfile(data_root, 'CallDetectionBackup');
behavior_folder = fullfile(data_root, 'Behavior backups');
npx_folder = fullfile(data_root, 'NPX data', 'NPX raw data');
area_limit_table = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\Area_limits_GoodLooking.xlsx';

%% Load PSTH — example: DELTA play bout (combo file row DELTA #1)
% Other combinations: swap filenames per combo file (see combo_file).

disp('loading')
load(fullfile(saving_folder, 'psth_structure_delta_updated.mat'), 'psth_structure'); % Struct of struct with one day each. Struct(1)  has band power and play bout.
load(fullfile(saving_folder, 'animal_names_delta_updated.mat'), 'animal_names');     % Animals names correspond to structures
% load(fullfile(saving_folder, 'psth_structure_delta_olny_aggr.mat'), 'psth_structure');
% load(fullfile(saving_folder, 'animal_names_delta_olny_aggr.mat'), 'animal_names');



%% mergind psth data 
smooth_wind = 20;
baseline_range = [-5 -2];
animal_label = {'B1D1','B1S3','B2S2','B3D2', 'B4S2', 'B4D4'};
electorde_numner = [1 2];
bin_size = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
psth_ranges = psth_structure(1).hist_range;
wrap_range = psth_structure(1).range_time_wrap;
time = psth_ranges(1):bin_size:psth_ranges(2)+bin_size;
baseline_index = time<baseline_range(2) & time>baseline_range(1);
baseline_index_time_wrap = 1:round((abs(wrap_range(1))/bin_size));

all_psth_onset          = [];
all_psth_onset_behavior = [];
all_psth_onset_only_playobut    = [];

all_psth_offset         = [];
all_psth_tw             = [];
all_psth_tw_3points     = [];
all_play_bouts          = [];
time_wrap_time          = [(baseline_index_time_wrap*bin_size) + wrap_range(1),linspace(0,1,psth_structure(1).n_bins_time_wrap),1 + (1:round((abs(wrap_range(2))/bin_size)))*bin_size];
time_wrap_3_points      = [(baseline_index_time_wrap*bin_size) + wrap_range(1),linspace(0,1-1/psth_structure(1).n_bins_time_wrap,psth_structure(1).n_bins_time_wrap), ...
    linspace(1,2-1/psth_structure(1).n_bins_time_wrap,psth_structure(1).n_bins_time_wrap),2 + (1:round((abs(wrap_range(2))/bin_size)))*bin_size];


animal_index = [];
electrode_index = []
session_index = [];
play_bout_numbers = [];
for j=1:numel(psth_structure)

    if contains(animal_names{j},animal_label)

        animal_num      = find(cell2mat(cellfun(@(x) contains(animal_names{j},x), animal_label, 'UniformOutput',false)));
        electrode_num   = animal_names{j,2}  ;
        this_animal_playbouts = psth_structure(j).play_bouts_table;
        this_animal_lengths = diff(this_animal_playbouts');

        this_psth_onset         = psth_structure(j).play_bout_onset;
        this_psth_onset_onlypb  = this_psth_onset;

        animal_index = [animal_index;repmat(animal_num,size(this_psth_onset,1),1)];
        electrode_index = [electrode_index;repmat(electrode_num,size(this_psth_onset,1),1)];
         session_index = [session_index;repmat(j,size(this_psth_onset,1),1)];
         play_bout_numbers = [play_bout_numbers;(1:size(this_psth_onset,1))'];

        for trial=1:size(this_psth_onset,1)
            this_psth_onset(trial,:) = ( this_psth_onset(trial,:) - mean( this_psth_onset(trial,baseline_index)))/std( this_psth_onset(trial,baseline_index));
            this_psth_onset(trial,:) = movmean(this_psth_onset(trial,:), smooth_wind);
            this_psth_onset_onlypb(trial,:) = this_psth_onset(trial,:);
            % this_psth_onset(trial,time> this_animal_lengths(trial)) = NaN;
        end
        all_psth_onset      = [all_psth_onset; this_psth_onset];
        all_psth_onset_only_playobut = [all_psth_onset_only_playobut; this_psth_onset_onlypb];

        this_psth_onset     = psth_structure(j).play_bout_onset;
        this_psth_offset    = psth_structure(j).play_bout_offset;
        for trial=1:size(this_psth_offset,1)
            this_psth_offset(trial,:) = ( this_psth_offset(trial,:) - mean( this_psth_onset(trial,baseline_index)))/std( this_psth_onset(trial,baseline_index));
            this_psth_offset(trial,:) = movmean(this_psth_offset(trial,:), smooth_wind);
        end
        all_psth_offset = [all_psth_offset; this_psth_offset];

        this_psth_tw = psth_structure(j).play_bout_tw_this;
        for trial=1:size(this_psth_tw,1)
            this_psth_tw(trial,:) = ( this_psth_tw(trial,:) - mean( this_psth_tw(trial,baseline_index_time_wrap)))/std( this_psth_tw(trial,baseline_index_time_wrap));
        end
        all_psth_tw = [all_psth_tw; this_psth_tw];


        this_psth_tw = psth_structure(j).three_point_tw;
        for trial=1:size(this_psth_tw,1)
            this_psth_tw(trial,:) = ( this_psth_tw(trial,:) - mean( this_psth_tw(trial,baseline_index_time_wrap)))/std( this_psth_tw(trial,baseline_index_time_wrap));
        end
        all_psth_tw_3points = [all_psth_tw_3points; this_psth_tw];



        this_psth_ab = psth_structure(j).animal_behavior_onset;
        for trial=1:size(this_psth_ab,1)
            this_psth_ab(trial,:) = ( this_psth_ab(trial,:) - mean( this_psth_ab(trial,baseline_index_time_wrap)))/std( this_psth_ab(trial,baseline_index_time_wrap));
            this_psth_ab(trial,:) = movmean(this_psth_ab(trial,:), smooth_wind);
        end
       
        all_psth_onset_behavior = [all_psth_onset_behavior; this_psth_ab(:, 1:4000)];

        all_play_bouts = [all_play_bouts;this_animal_playbouts];
    end
end

play_bout_length = diff(all_play_bouts')';


[sorted_play_bout_length, order] = sort(play_bout_length);
%% plot single animals  (and obtain mean response per animal)
X_lim = [-3 3]
figure
min_length = 0;
stacked_mean_onset = [];
% what2plot = all_psth_onset_behavior; %Select what to plot
what2plot = all_psth_onset;
% what2plot = all_psth_tw_3points;
zscore_limit = Inf;
time2use = time;
% time2use = time_wrap_3_points;

artifcat_removal = max(abs(what2plot(:,time2use>X_lim(1))),[],2,'omitmissing')<4;
for an= 1:numel(animal_label)
    animal_bool = animal_index==an;
    length_bool = play_bout_length>min_length;
    electrode_bool = electrode_index==1 & artifcat_removal;
    [sorted_play_bout_length, order] = sort(play_bout_length(animal_bool & length_bool & electrode_bool,:));
    subplot(5,numel(animal_label),(1:numel(animal_label):2*numel(animal_label)) + an-1)

    array = what2plot(animal_bool & length_bool & electrode_bool,:);
    imagesc(time2use,1:numel(sorted_play_bout_length),array(order,:) )
    xlim(X_lim)
    clim([-2 2])
    axis xy
    hold on
    plot([0 0],[1 numel(sorted_play_bout_length)], 'w')
    plot(sorted_play_bout_length,1:numel(sorted_play_bout_length), 'w')
    title(animal_label{an})

    subplot(5,numel(animal_label),((2*numel(animal_label) + 1):numel(animal_label):5*numel(animal_label)) + an-1)

    [~, ~, ci]  = ttest(array);
    no_nan = ~any(isnan(ci));
  
    fill([time2use(no_nan) fliplr(time2use(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
    hold on
    plot(time2use,mean(array, 'omitmissing'), 'k')
    xlim(X_lim)
    stacked_mean_onset = [stacked_mean_onset;mean(array, 'omitmissing')];
end
%% plot lfp together with power
min_length = .0;
animal =2;
session_number = 6;
y_lim = [-3 3];

animal_bool = animal_index==animal;
length_bool = play_bout_length>min_length;
electrode_bool = electrode_index==1 & artifcat_removal;
this_sessions = session_index(animal_bool & length_bool & electrode_bool);
[sorted_play_bout_length, order] = sort(play_bout_length(animal_bool & length_bool & electrode_bool,:));

all_session_6 = find(this_sessions(order)==session_number)';
plot_beh = true;
plot_playbout = true;
trial_n         = 21;
for trial = all_session_6(trial_n)


    % x_lim = [time2use(1) time2use(end)];
    x_lim = [-2 12];
    what2plot = all_psth_onset;
    figure('units','normalized','outerposition',[0 0.5 1 .5]);
    subplot(2,1,1)
    animal_bool = animal_index==animal;
    length_bool = play_bout_length>min_length;
    electrode_bool = electrode_index==1 & artifcat_removal;
    [sorted_play_bout_length, order] = sort(play_bout_length(animal_bool & length_bool & electrode_bool,:));


    array = what2plot(animal_bool & length_bool & electrode_bool,:);
    plot(time2use, array(order(trial),:))
    hold on





    this_sessions = session_index(animal_bool & length_bool & electrode_bool);
    title(animal_names{this_sessions(order(trial)),1})
    this_anima_play_bout = all_play_bouts(animal_bool & length_bool & electrode_bool,:);
    this_playbout_numbers = play_bout_numbers(animal_bool & length_bool & electrode_bool);
    this_trial_number = this_playbout_numbers(order(trial));

    implanted_animal =animal_names{this_sessions(order(trial)),1};
    implanted_animal = strsplit(implanted_animal,' ');
    implanted_animal = implanted_animal{end};
    Behavior = psth_structure(this_sessions(order(trial))).Behavior;

    animal_list = unique(Behavior.Animal);
    animal_list(ismember(animal_list,{'Session_structure'}))=[];
    partners = animal_list;
    partners(ismember(partners, implanted_animal)) = [];
    animal_cell = {{implanted_animal},partners};
    play_bouts = psth_structure(this_sessions(order(trial))).play_bouts_table;



    if plot_beh
        beh_list = find(Behavior.End>this_anima_play_bout(order(trial),1)+time2use(1) & Behavior.End<this_anima_play_bout(order(trial),1)+time2use(end));
        beh_list(ismember(Behavior.Animal(beh_list),'Session_structure')) = [];
        start_end = Behavior(beh_list,{'Start','End'})-this_anima_play_bout(1)
        for bn= 1:numel(beh_list)
            animal_n = find(cell2mat(cellfun(@(x) ismember(Behavior.Animal(beh_list(bn)),x),animal_cell,'UniformOutput',false)));
            fill([Behavior.Start(beh_list(bn)) Behavior.End(beh_list(bn))  Behavior.End(beh_list(bn))  Behavior.Start(beh_list(bn)) ]-this_anima_play_bout(order(trial),1), ...
                [y_lim([1 1]), y_lim([1 1])+range(y_lim)/2] + (range(y_lim)/2)*(animal_n==1), 'r', 'FaceAlpha',.5, 'EdgeColor', 'k')
            text(Behavior.Start(beh_list(bn))-this_anima_play_bout(order(trial),1),y_lim(1)+range(y_lim)/2 + (range(y_lim)/2)*(animal_n==1),Behavior.Type(beh_list(bn)) )
        end

    end
    if plot_playbout

        play_bout_list = find(play_bouts(:,2)>this_anima_play_bout(order(trial),1)+time2use(1) & play_bouts(:,2)<this_anima_play_bout(order(trial),1)+time2use(end));

        start_end = play_bouts(play_bout_list,:)-this_anima_play_bout(order(trial));
        for bn= 1:numel(play_bout_list)
            fill(start_end(bn,[1 2 2 1]), ...
                y_lim([1 1 2 2]), 'g', 'FaceAlpha',.5, 'EdgeColor', 'k')
            text(start_end(bn,1),y_lim(2),'PlayBOut' )
        end
    end





    plot([0 0],y_lim, 'g')
    plot(sorted_play_bout_length([trial trial]),y_lim, 'g')
    ylim(y_lim)
    xlim(x_lim)
    pause(.1)
end
%% load lfp 

this_animal_dir = [npx_folder,'\',animal_names{this_sessions(order(trial)),1}];

disp('LOADING LFP')
if exist([npx_folder,'\',animal_names{this_sessions(order(trial)),1},'\LFP_PAG.mat'], 'file')==2
    % Preprocessed LFP file exists
    NPX_Type = 2;
    load([npx_folder,'\',animal_names{this_sessions(order(trial)),1},'\LFP_PAG.mat'], 'LFP')
elseif exist([npx_folder,'\',animal_names{this_sessions(order(trial)),1},'\LFP_PAG.dat'], 'file')==2
    % Load raw binary LFP file
    NPX_Type = 1;
    file_pointer = fopen([npx_folder,'\',animal_names{this_sessions(order(trial)),1},'\LFP_PAG.dat'], 'r');
    LFP = fread(file_pointer,'int16');
    LFP = reshape(LFP, 384, numel(LFP)/384);
end
disp('LFP LOADED')

%% -------------------- SELECT PAG CHANNEL(S) --------------------
disp('Loading Channel Map')
hard_coded_x_coords = [8 40;258 290; 508 540; 758 790];
area_limit = readtable(area_limit_table);

animal_code_params = strsplit(animal_names{this_sessions(order(trial)),1}, ' ');
animal_batch       = animal_code_params{1};
date               = animal_code_params{2};
repeated_animal    = animal_code_params{3};
% Build animal identifier for area selection
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal,animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName,this_animal),:);

if NPX_Type == 1
    % Raw LFP: select channel range for LPAG region
    PAG_channels = area_limit{ismember(area_limit.area, {'LPAG'}), {'ch_start', 'ch_end'}};
    PAG_channels = str2double(PAG_channels);
    channel_Range = [min(PAG_channels(:)) max(PAG_channels(:))];
    mid_PAG_channel = round(mean(channel_Range));
else
    % Preprocessed: use ChannelMap.mat to locate mid-PAG channel
    load([npx_folder,'\',animal_code,'\ChannelMap.mat'], 'xcoords', 'ycoords','chanMap')
    Y_Range = area_limit{ismember(area_limit.area, {'LPAG'}), {'ProbeNum','depth_start', 'depth_end'}};
    mid_PAG_channel = nan(size(Y_Range,1),1);
    figure
    plot(xcoords,ycoords, 'k.'); hold on
    for j=1:size(Y_Range,1)
        this_indexes = ycoords>=Y_Range(j,2) & ycoords<=Y_Range(j,3) & ismember(xcoords,hard_coded_x_coords(Y_Range(j,1),:));
        all_locs = [xcoords(this_indexes) ycoords(this_indexes)];
        plot(all_locs(:,1),all_locs(:,2), 'r.')
        mean_loc = mean(all_locs);
        [~, closest_channel]= min(sum(abs([xcoords ycoords]-repmat(mean_loc,numel(ycoords),1)),2));
        plot(xcoords(closest_channel), ycoords(closest_channel), 'xb')
        mid_PAG_channel(j) = chanMap(closest_channel);
    end
end

%% load synch function

load([synch_directory,'\', animal_names{this_sessions(order(trial)),1}, '\synch_model_video2NPX.mat'])
load([synch_directory,'\', animal_names{this_sessions(order(trial)),1}, '\synch_model_audio2NPX.mat'])
%% load raw data
lfp_ylim_left = [-2000 2000];
lfp_ylim_right = [-1000 1000];

low_pass_freq = 50;


start_event = this_anima_play_bout(order(trial),1)+time2use(1); %in video time
end_event   =this_anima_play_bout(order(trial),1)+time2use(end); %in video time


sr = 2500;
start_event_index = round(start_event*sr);
end_event_index     = round(end_event*sr);
slected_lfp = LFP(mid_PAG_channel,start_event_index:end_event_index);
slected_lfp = movmean(slected_lfp,sr/low_pass_freq );
time = (start_event_index:end_event_index)/sr;
theta_band = [6 12];
delta_band = [1 5];

subplot(2,1,2)
hold off


filter_order = 2000;
% Hd_theta = designfilt('bandpassfir', ...
% 'FilterOrder', filter_order, ...
% 'CutoffFrequency1', theta_band(1), ...
% 'CutoffFrequency2', theta_band(2), ...
% 'SampleRate', sr, ...
% 'DesignMethod', 'window', ...
% 'Window', 'hamming');

Hd_delta = designfilt('bandpassfir', ...
'FilterOrder', filter_order, ...
'CutoffFrequency1', delta_band(1), ...
'CutoffFrequency2', delta_band(2), ...
'SampleRate', sr, ...
'DesignMethod', 'window', ...
'Window', 'hamming');

% filtered_signal_theta = filtfilt(Hd_theta.Coefficients, 1, LFP(mid_PAG_channel,start_event_index:end_event_index));
filtered_signal_delta = filtfilt(Hd_delta.Coefficients, 1, LFP(mid_PAG_channel,start_event_index:end_event_index));

delta_ampl = abs(hilbert(filtered_signal_delta));
hold off
yyaxis left
plot(time-this_anima_play_bout(order(trial),1),slected_lfp,':k')
ylim(lfp_ylim_left)
yyaxis right
hold on
% plot(time-precited_playbout_onset,filtered_signal_theta,'b')
plot(time-this_anima_play_bout(order(trial),1),filtered_signal_delta,'r')

    plot(time-this_anima_play_bout(order(trial),1), delta_ampl, 'b')
    y_lim = ylim;
    for bn= 1:numel(beh_list)

    fill([Behavior.Start(beh_list(bn)) Behavior.End(beh_list(bn))  Behavior.End(beh_list(bn))  Behavior.Start(beh_list(bn)) ]-this_anima_play_bout(order(trial),1), ...
        y_lim([1 1 2 2]), 'r', 'FaceAlpha',.5, 'EdgeColor', 'k')
    text(Behavior.Start(beh_list(bn))-this_anima_play_bout(order(trial),1),y_lim(round(.5*(-1)^bn +1.5)),Behavior.Type(beh_list(bn)) )
end
xlim(x_lim)
ylim(lfp_ylim_right)

%% Load explained variance — example: DELTA + THETA CV (combo file DELTA #7, THETA #3)
% GAMMA #3 → cvResults_mean_calls_gamma_AlllVar_play_baut.mat
% load(fullfile(saving_folder, 'cvResults_single_calls.mat'), 'cvResults')
% load(fullfile(saving_folder, 'cvResults_mean_calls_delta_v2.mat'), 'cvResults')
load(fullfile(saving_folder, 'cvResults_mean_calls_delta_v2_AlllVar.mat'), 'cvResults')
load(fullfile(saving_folder, 'cvResults_mean_calls_theta_v2_AlllVar.mat'), 'cvResults')





k = cvResults.k;
predictors = cvResults.predictors;


R2_full_allfolds =cvResults.Call.r2_full;
avgR2_full = mean(R2_full_allfolds);

y_for_sig = 20;
figure('units','normalized','outerposition',[.5 0 .5 1]);
hold on
randx = .25*(rand(k,2)-.5);


re_order = 1:numel(predictors);
p_vals = nan(numel(re_order),1);
mean_diff = nan(numel(predictors),1);


for j=1:numel(predictors)
    % Extract R² values for this predictor
    r2_full = cvResults.(predictors{re_order(j)}).r2_full;
    r2_reduced = cvResults.(predictors{re_order(j)}).r2_reduced;
    mean_diff(j) =mean( r2_full - r2_reduced)/mean(r2_full); % difference per fold
end
% [~, re_order] =  sort(mean_diff, 'descend');

re_order = [14 7 1 15  2 13 6 4 3 8 9 10 5 11 12]; % manual order
% re_order =1:numel(re_order) % manual order

for j=1:numel(re_order)
    % Extract R² values for this predictor
    r2_full = cvResults.(predictors{re_order(j)}).r2_full;
    r2_reduced = cvResults.(predictors{re_order(j)}).r2_reduced;
    r2_diff = r2_full - r2_reduced; % difference per fold

    % Plot swarm of differences
    swarmchart((2*j-1)*ones(k,1), 100*r2_diff/avgR2_full, 'o', 'MarkerFaceColor','flat', 'MarkerFaceAlpha',.25);
    plot((2*j-1), 100*mean(r2_diff/avgR2_full), '_r','MarkerSize',10,'LineWidth',2);

    % Paired t-test (full vs reduced)
    % [~,p_val] = ttest(r2_diff);
    p_val = signrank(r2_diff);
    p_val = p_val * numel(predictors); % Bonferroni correction
    p_vals(re_order(j)) = p_val;

    % Annotate significance stars
    if p_val<0.001
        text(2*j -1, y_for_sig, '* * *','HorizontalAlignment','center')
    elseif p_val<0.01
        text(2*j -1, y_for_sig, ' * * ','HorizontalAlignment','center')
    elseif p_val<0.05
        text(2*j -1, y_for_sig, '  *  ','HorizontalAlignment','center')
    else
        text(2*j -1, y_for_sig, ' n.s ','HorizontalAlignment','center')
    end

    % Annotate actual p-value
    text(2*j -1, y_for_sig+ 5, num2str(p_val,'%.5f'), 'Rotation',45,'HorizontalAlignment','center')
end

% Reference line at 0
plot([0 2*numel(predictors)], [0 0], 'k')

% Formatting
xticks(1:2:2*numel(predictors))
xticks(1:2:2*numel(predictors))
xticklabels(predictors(re_order))
ylabel('\Delta R^2 (Full - Reduced)')
ylim([-1 y_for_sig+10])
title('Cross-validated \DeltaR^2 when removing each predictor')
set(gca, 'FontSize', 24)

%% Load PSTH for LFP power panels — example: DELTA play bout (combo file DELTA #1)
disp('Loading')
load(fullfile(saving_folder, 'psth_structure_delta_updated.mat'), 'psth_structure');
load(fullfile(saving_folder, 'animal_names_delta_updated.mat'), 'animal_names');
disp('Loading Ready')

%% merging_psth
smooth_wind = 20;
baseline_range = [-2 0]
animal_label = {'B1D1','B1S3','B2S2','B3D2', 'B4S2', 'B4D4'};
electorde_numner = [1 2];
bin_size = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
psth_ranges = psth_structure(1).hist_range;
wrap_range = psth_structure(1).range_time_wrap;
time = psth_ranges(1):bin_size:psth_ranges(2)+bin_size;
baseline_index = time<baseline_range(2) & time>baseline_range(1);
baseline_index_time_wrap = 1:round((abs(wrap_range(1))/bin_size));
all_psth_onset          = [];
all_psth_onset_behavior = [];
all_psth_onset_only_playobut    = [];

all_psth_offset         = [];
all_psth_tw             = [];
all_psth_tw_3points     = [];
all_play_bouts          = [];
time_wrap_time          = [(baseline_index_time_wrap*bin_size) + wrap_range(1),linspace(0,1,psth_structure(1).n_bins_time_wrap),1 + (1:round((abs(wrap_range(2))/bin_size)))*bin_size];
time_wrap_3_points      = [(baseline_index_time_wrap*bin_size) + wrap_range(1),linspace(0,1-1/psth_structure(1).n_bins_time_wrap,psth_structure(1).n_bins_time_wrap), ...
    linspace(1,2-1/psth_structure(1).n_bins_time_wrap,psth_structure(1).n_bins_time_wrap),2 + (1:round((abs(wrap_range(2))/bin_size)))*bin_size];


animal_index = [];
electrode_index = []
for j=1:numel(psth_structure)

    if contains(animal_names{j,1},animal_label)

        animal_num      = find(cell2mat(cellfun(@(x) contains(animal_names{j},x), animal_label, 'UniformOutput',false)));
        electrode_num   = animal_names{j,2}  ;
        this_animal_playbouts = psth_structure(j).play_bouts_table;
        this_animal_lengths = diff(this_animal_playbouts');

        this_psth_onset         = psth_structure(j).play_bout_onset;
        this_psth_onset_onlypb  = this_psth_onset;

        animal_index = [animal_index;repmat(animal_num,size(this_psth_onset,1),1)];
        electrode_index = [electrode_index;repmat(electrode_num,size(this_psth_onset,1),1)];

        for trial=1:size(this_psth_onset,1)
            this_psth_onset(trial,:) = ( this_psth_onset(trial,:) - mean( this_psth_onset(trial,baseline_index)))/std( this_psth_onset(trial,baseline_index));
            this_psth_onset(trial,:) = movmean(this_psth_onset(trial,:), smooth_wind);
            this_psth_onset_onlypb(trial,:) = this_psth_onset(trial,:);
            this_psth_onset_onlypb(trial,time> this_animal_lengths(trial)) = NaN;
        end
        all_psth_onset      = [all_psth_onset; this_psth_onset];
        all_psth_onset_only_playobut = [all_psth_onset_only_playobut; this_psth_onset_onlypb];

        this_psth_onset     = psth_structure(j).play_bout_onset;
        this_psth_offset    = psth_structure(j).play_bout_offset;
        for trial=1:size(this_psth_offset,1)
            this_psth_offset(trial,:) = ( this_psth_offset(trial,:) - mean( this_psth_onset(trial,baseline_index)))/std( this_psth_onset(trial,baseline_index));
            this_psth_offset(trial,:) = movmean(this_psth_offset(trial,:), smooth_wind);
        end
        all_psth_offset = [all_psth_offset; this_psth_offset];





        all_play_bouts = [all_play_bouts;this_animal_playbouts];
    end
end

play_bout_length = diff(all_play_bouts')';


[sorted_play_bout_length, order] = sort(play_bout_length);

%% plot single animals  (and obtain mean response per animal)
X_lim = [-2 2]
figure
min_length = .0;
stacked_mean_onset = [];
% what2plot = all_psth_onset_behavior; %Select what to plot
what2plot = all_psth_onset_only_playobut;
for an= 1:numel(animal_label)
    animal_bool = animal_index==an;
    length_bool = play_bout_length>min_length;
    electrode_bool = electrode_index==1;
    [sorted_play_bout_length, order] = sort(play_bout_length(animal_bool & length_bool & electrode_bool,:));
    subplot(5,numel(animal_label),(1:numel(animal_label):2*numel(animal_label)) + an-1)

    array = what2plot(animal_bool & length_bool & electrode_bool,:);
    imagesc(time,1:numel(sorted_play_bout_length),array(order,:) )
    xlim(X_lim)
    clim([-2 2])
    axis xy
    hold on
    plot([0 0],[1 numel(sorted_play_bout_length)], 'w')
    plot(sorted_play_bout_length,1:numel(sorted_play_bout_length), 'w')
    title(animal_label{an})

    subplot(5,numel(animal_label),((2*numel(animal_label) + 1):numel(animal_label):5*numel(animal_label)) + an-1)

    [~, ~, ci]  = ttest(array);
    no_nan = ~any(isnan(ci));
    fill([time(no_nan) fliplr(time(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
    hold on
    plot(time,mean(array, 'omitmissing'), 'k')
    xlim(X_lim)
    stacked_mean_onset = [stacked_mean_onset;mean(array, 'omitmissing')];
end

%% Load mixed-model results — example: DELTA play bout (combo file DELTA #1)
% load(fullfile(saving_folder, 'results_play_bout.mat'), 'results'); % all time points
load(fullfile(saving_folder, 'results_play_bout_PBonly_zscore4_updated.mat'), 'results'); % play-bout-only bins
est = results.est;
ci = results.ci;
pvals_fdr = results.pvals;
%% plot all together

X_lim = [-1 3];
alpha = 0.05;
figure
min_length = .0;
fill_lim = [.35 .4];

length_bool = play_bout_length>min_length;
electrode_bool = electrode_index==1;
[sorted_play_bout_length, order] = sort(play_bout_length( length_bool & electrode_bool,:));

subplot(2,1,1)
array = what2plot( length_bool & electrode_bool,:);
imagesc(time,1:numel(sorted_play_bout_length),array(order,:) )

xlim(X_lim)
clim([-1 2])
axis xy
hold on
plot([0 0],[1 numel(sorted_play_bout_length)], 'w')
plot(sorted_play_bout_length,1:numel(sorted_play_bout_length), 'w')
title('All data')

subplot(2,1,2)

[~, ~, ci]  = ttest(array);
fill([time fliplr(time)], [ci(1,:) fliplr(ci(2,:))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
hold on
plot(time,mean(array, 'omitmissing'), 'k')
xlim(X_lim)

plot(time,stacked_mean_onset, 'k:')
hold on
plot(time,mean(stacked_mean_onset), 'k:')
xlim(X_lim)

sig_idx = pvals_fdr < alpha;
limited_time = results.time;

borders = [ 0 sig_idx' 0];
beg_index = find(diff(borders)==1);

end_index = find(diff(borders)==-1)-1;

for bn = 1:numel(beg_index)
    time2plot = limited_time(beg_index(bn):end_index(bn))';
    fill([time2plot fliplr(time2plot)], [time2plot*0+fill_lim(1) time2plot*0+fill_lim(2)], 'r', 'EdgeColor','none')
end
plot(limited_time,est)

hold on
y_lim =ylim;
plot([0 0],y_lim, 'k' )
ylim tight
%% Call responses

%% Load call PSTH — example: GAMMA calls (combo file GAMMA #2)
load(fullfile(saving_folder, 'psth_structure_call_gamma.mat'), 'psth_structure');
load(fullfile(saving_folder, 'animal_names_call_gamma.mat'), 'animal_names');

%% mergin data

smooth_wind             = 20;
baseline_range          = [-2 0]
animal_label            = {'B1D1','B1S3','B2S2','B3D2', 'B4S2'};
electorde_number        = [1 2];
bin_size                = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
psth_ranges             = psth_structure(1).hist_range;
time                    = psth_ranges(1):bin_size:psth_ranges(2)+bin_size;
baseline_index          = time<baseline_range(2) & time>baseline_range(1);
all_psth_onset          = [];
all_psth_offset         = [];
all_onset_regressors    = [];
all_offset_regressors   = [];
all_Calls               = [];

electrode_index  = [];
animal_index = [];
for j=1:numel(psth_structure)

    if contains(animal_names{j},animal_label)

        animal_num      = find(cell2mat(cellfun(@(x) contains(animal_names{j},x), animal_label, 'UniformOutput',false)));
        electrode_num   = animal_names{j,2};

        this_psth_onset         = psth_structure(j).call_onset;
        animal_index = [animal_index;repmat(animal_num,size(this_psth_onset,1),1)];
        electrode_index = [electrode_index;ones(size(this_psth_onset,1),1)*electrode_num];
        for trial=1:size(this_psth_onset,1)
            this_psth_onset(trial,:) = ( this_psth_onset(trial,:) - mean( this_psth_onset(trial,baseline_index)))/std( this_psth_onset(trial,baseline_index));
            this_psth_onset(trial,:) = movmean(this_psth_onset(trial,:), smooth_wind);
        end
        all_psth_onset      = [all_psth_onset; this_psth_onset];

        this_psth_onset     = psth_structure(j).call_onset;
        this_psth_offset    = psth_structure(j).call_offset;
        for trial=1:size(this_psth_offset,1)
            this_psth_offset(trial,:) = ( this_psth_offset(trial,:) - mean( this_psth_onset(trial,baseline_index)))/std( this_psth_onset(trial,baseline_index));
            this_psth_offset(trial,:) = movmean(this_psth_offset(trial,:), smooth_wind);
        end
        all_psth_offset = [all_psth_offset; this_psth_offset];


        all_onset_regressors = [all_onset_regressors; psth_structure(j).call_onset_regressor];
        all_offset_regressors = [all_offset_regressors; psth_structure(j).call_onset_regressor];


        all_Calls = [all_Calls;psth_structure(j).CallStats];
    end
end

call_lengths = all_Calls.CallLengths;


[sorted_call_lengths, order] = sort(call_lengths);

disp('merge done')
%% ploting each animal (and obtain mean response per animal)

stacked_mean = [];
X_lim = [-2 .5];
figure
min_length = .0;
for an= 1:numel(animal_label)
    animal_bool = animal_index==an;
    length_bool = call_lengths>min_length;
    electrode_bool = electrode_index==1;
    [sorted_call_lengths, order] = sort(call_lengths(animal_bool & length_bool & electrode_bool,:));
    subplot(5,numel(animal_label),(1:numel(animal_label):2*numel(animal_label)) + an-1)

    array = all_psth_onset(animal_bool & length_bool & electrode_bool,:);
    imagesc(time,1:numel(sorted_call_lengths),array(order,:) )
    xlim(X_lim)
    clim([-2 2])
    axis xy
    hold on
    plot([0 0],[1 numel(sorted_call_lengths)], 'w')
    plot(sorted_call_lengths,1:numel(sorted_call_lengths), 'w')
    title(animal_label{an})

    subplot(5,numel(animal_label),((2*numel(animal_label) + 1):numel(animal_label):5*numel(animal_label)) + an-1)

    [~, ~, ci]  = ttest(array);
    no_nan = ~any(isnan(ci));
    fill([time(no_nan) fliplr(time(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
    hold on
    plot(time,mean(array,'omitmissing'), 'k')
    yyaxis right
    plot(time, mean(all_onset_regressors(animal_bool & length_bool,:)), 'r')
    xlim(X_lim)

    stacked_mean = [stacked_mean;mean(array,'omitmissing')];
end

%% Load call mixed-model results — example: GAMMA calls (combo file GAMMA #2)
load(fullfile(saving_folder, 'results_call_updated_gamma.mat'), 'results');

limited_time = results.time;
est = results.est;
ci = results.ci;
pvals_fdr = results.pvals;
%% now plot all together
figure
an=2;
alpha = 0.01;
fill_lim = [.15 .18]
X_lim = [-1 2];
min_length = 0;
animal_bool = animal_index==an;
length_bool = call_lengths>min_length;
electrode_bool = electrode_index==1;
[sorted_call_lengths, order] = sort(call_lengths( length_bool & electrode_bool,:));
subplot(2,1,1)

array = all_psth_onset( length_bool & electrode_bool,:);
imagesc(time,1:numel(sorted_call_lengths),array(order,:) )
xlim(X_lim)
clim([-2 2])
axis xy
hold on
plot([0 0],[1 numel(sorted_call_lengths)], 'w')
plot(sorted_call_lengths,1:numel(sorted_call_lengths), 'w')
title(animal_label{an})

subplot(2,1,2)

[~, ~, ci]  = ttest(array);
no_nan = ~any(isnan(ci));
fill([time(no_nan) fliplr(time(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
hold on
plot(time,mean(array, 'omitmissing'), 'k')



plot(time,stacked_mean, 'k:')
hold on
plot(time,mean(stacked_mean), 'k:')

sig_idx = pvals_fdr < alpha;
limited_time = results.time;

borders = [ 0 sig_idx' 0];
beg_index = find(diff(borders)==1);

end_index = find(diff(borders)==-1)-1;

for bn = 1:numel(beg_index)
    time2plot = limited_time(beg_index(bn):end_index(bn))';
    fill([time2plot fliplr(time2plot)], [time2plot*0+fill_lim(1) time2plot*0+fill_lim(2)], 'r', 'EdgeColor','none', 'FaceAlpha',.2)
end

hold on
y_lim =ylim;
plot([0 0],y_lim, 'k' )
ylim tight
xlim(X_lim)
yyaxis right
plot(time, mean(all_onset_regressors( length_bool,:)), 'r')
xlim(X_lim)

%% plot theta and speed


load([saving_folder,'\psth_structure_speed_theta.mat'],'psth_structure')
load([saving_folder,'\animal_names_speed_theta.mat'],'animal_names')

%% obtain theta increas eduirng play (mixed effect linear model)

stasts_table = cell(numel(psth_structure),6);

for j=1:numel(psth_structure)

    stasts_table(j,:) = { psth_structure(j).lm.Coefficients.Estimate(4), ...
        psth_structure(j).lm.Coefficients.pValue(4),...
        psth_structure(j).lm.Coefficients.tStat(4),...
        animal_names{j,:}};

end

stasts_table =cell2table(stasts_table);
stasts_table.Properties.VariableNames = {'Estimate','pValue','tStat','Animal','Partner','Electrode'};

stasts_table = stasts_table(stasts_table.Electrode==1,:);

lme = fitlme(stasts_table, 'Estimate ~ 1 + (1|Animal)');
coef = fixedEffects(lme);       % estimated mean
ci_estimate = coefCI(lme);               % confidence interval
pval = lme.Coefficients.pValue; % p-value for intercept

fprintf('Mean = %.3f, 95%% CI [%.3f, %.3f], p = %.4f\n', coef, ci(1), ci(2), pval);



%% Compute speed and theta  
n_grid = 101;
distance_grid = linspace(-10, 10, n_grid);

data_together = [];
mean_powers_play    = nan(numel(psth_structure) ,n_grid);
mean_powers_noplay  = nan(numel(psth_structure) ,n_grid);
for j=1:numel(psth_structure)


    tbl = psth_structure(j).model_data;
    tbl = tbl(~isnan(tbl.Speed),:);

    Y = tbl.Power;
    Y = (Y - mean(Y, 'omitmissing'))/std(Y, 'omitmissing');
    X = tbl.Speed;
    Xa = X(tbl.Play=='true');
    Ya = Y(tbl.Play=='true');


    Xb = X(tbl.Play=='false');
    Yb =  Y(tbl.Play=='false');


    data_together = [data_together;[X Y]];

    edges = [distance_grid Inf];

    % Assign each sample to a bin
    [~,~,binA] = histcounts(Xa, edges);
    [~,~,binB] = histcounts(Xb, edges);

    % Pre-allocate as NaN
    meanA = nan(size(distance_grid));
    meanB = nan(size(distance_grid));

    % Compute mean per bin (only for bins with data)
    for i = 1:numel(distance_grid)
        meanA(i) = mean(Ya(binA == i), 'omitnan');
        meanB(i) = mean(Yb(binB == i), 'omitnan');
    end


    mean_powers_play(j,:)    = meanA;
    mean_powers_noplay(j,:)  = meanB;
end
%% now plot
[bin_count, ~, samples2plot] = histcounts(data_together(:,1), -5:.1:10);

figure
bins_with_values = unique(samples2plot)'
max_count = 100;
min_count = 20;

indexes2plot = [];

for j=bins_with_values
    indexes= find(samples2plot==j);
    if numel(indexes)>=max_count

        indexes  = datasample(indexes,max_count );
    end
    if numel(indexes)>=min_count
        indexes2plot = [indexes2plot;indexes];
    end
end




plot(data_together(indexes2plot,1), data_together(indexes2plot,2), '.k', 'MarkerSize',    .1)
axis xy
