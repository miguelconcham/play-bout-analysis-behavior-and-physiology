function cc_structure = GENERATE_STRUCUTRE_animal_synch_cross_correlogram(directory,Hd_freq)
%% set parameters

%time range for the delta psth
time_range = [-5 10];
cross_correlogram_lenth = [-3 3];


%lfp sampling rate
sr = 2500;
%params mutual information

%params playbout
% non_play_behaviors   = {'PounceI',  'Bite'};
play_behaviors      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};
% play_behaviors      = non_play_behaviors;
bin_size            = 0.01;
conv_length         = 1;

%params spectrogram
wind_length         = 1;
wind_overlap        = .990;

freq_pow_range      = [1 4;
                       5 6;
                       7 12;
                       12 20;
                       20 40];

f                   = .1:.1:max(freq_pow_range(:));
% wind_length         = .250;
% wind_overlap        = .240;
% f                   = 5:.1:14;
% freq_pow_range      = [6 12];
%% load data: animal info 

animal_info =directory;
animal_info = strsplit(animal_info, '\');
animal_info = animal_info{end};
animal_info = strsplit(animal_info, ' ');



%% load  data: Synch model
file_loc    = fullfile(directory,'synch_model_video2NPX.mat');
load(file_loc, 'synch_model_video2NPX')

%% load  data: Behavior
Behavior_file = dir([directory, '/ELAN*']);
Behavior_file =Behavior_file.name;
file_loc    = fullfile(directory,Behavior_file);

animal_1                            = animal_info{1};
Behavior                            = readtable(file_loc);
Behavior(:,2)                       = [];
Behavior.Properties.VariableNames   = {'Animal', 'Start', 'End', 'Length', 'Type'};



Behavior.Type2          = Behavior.Type;
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_A', 'Pounce_B'}))      = {'Pounce'}; %% Merging behaviors to Type2
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'}))    = {'PounceI'};
Behavior.Type2(ismember( Behavior.Type2,''))                            = {'Other'};
Behavior(ismember(Behavior.Animal, 'Reversal'),:)                       = [];

animal_types            = unique(Behavior.Animal);

animal_types(ismember(animal_types,'Session_structure'))                =[];

Behavior.Start          = predict(synch_model_video2NPX, Behavior.Start);
Behavior.End            = predict(synch_model_video2NPX, Behavior.End);

%% load  data: LFP
file_loc    = fullfile(directory,'channel selection.txt');
channel_list = readtable(file_loc);

disp('LOADING animal 1')
animal_1_ch = channel_list.Var2(ismember(channel_list.Var1, animal_info{1}));
if exist([directory, '\LFP_' ,animal_info{1},'.mat'], "file")==2
    load ([directory, '\LFP_' ,animal_info{1},'.mat'], 'LFP')
    LFP_animal_1_PAG = LFP;
    clear LFP
else
    file_pointer                = fopen([directory, '\continuous_' ,animal_info{1},'.dat'], 'r');
    LFP_animal_1_PAG             = fread(file_pointer,'int16');

    LFP_animal_1_PAG             = reshape(LFP_animal_1_PAG, 384, numel(LFP_animal_1_PAG)/384);
end
lfp_animal_1 = double(LFP_animal_1_PAG(animal_1_ch,:));
clear LFP_animal_1_PAG

disp('LOADING animal 2')
animal_2_ch = channel_list.Var2(ismember(channel_list.Var1, animal_info{2}));
if exist([directory, '\LFP_' ,animal_info{2},'.mat'], "file")==2
    load ([directory, '\LFP_' ,animal_info{2},'.mat'], 'LFP')
    LFP_animal_2_PAG = LFP;
    clear LFP
else
    file_pointer                = fopen([directory, '\continuous_' ,animal_info{2},'.dat'], 'r');
    LFP_animal_2_PAG             = fread(file_pointer,'int16');
    LFP_animal_2_PAG             = reshape(LFP_animal_2_PAG, 384, numel(LFP_animal_2_PAG)/384);
end
lfp_animal_2 = double(LFP_animal_2_PAG(animal_2_ch,:));
clear LFP_animal_2_PAG

% play_song([],[],[])
%% Estimate Playbouts

max_lag_sec = max(cross_correlogram_lenth);
max_lag = round(max_lag_sec * sr);

config.Behavior         = Behavior;
config.repeated_animal  = animal_1;
config.animal_types     = animal_types        ;
config.play_behaviors   = play_behaviors      ;
config.beh_bin          = bin_size             ;
config.conv_length      = conv_length;
config.behavior_window  = 0;


[play_bouts_table]      = play_bout(config);



play_session = Behavior{ismember(Behavior.Type, 'Partners session'), {'Start', 'End'}};
if isempty(play_session)
    play_session = [min(Behavior.Start)-5 max(Behavior.End)+5];
end





%% same as before but now using power spectrum isntead of hilbert transform
disp('Estimating Spectrogram and mean delta power')


delta_animal_1 =   filtfilt(Hd_freq.Coefficients, 1,lfp_animal_1);
delta_animal_2 =   filtfilt(Hd_freq.Coefficients, 1,lfp_animal_2);
lfp_time       = (1:size(delta_animal_1,2))/sr;


% play_song([],[],[])


 %%

 cross_corr = nan(size(play_bouts_table,1), range(cross_correlogram_lenth)*sr +1);


for j=1:size(play_bouts_table, 1)

    beg_time = play_bouts_table(j,1);

    [~, loc] = min(abs(lfp_time-beg_time));

    index2exctact = round((loc+cross_correlogram_lenth(1)*sr):(loc+cross_correlogram_lenth(2)*sr));
    possible_index = ismember(index2exctact,1:numel(lfp_time));

   delta_this_pb_1 =  (delta_animal_1(index2exctact(possible_index)));
   delta_this_pb_2 =  (delta_animal_2(index2exctact(possible_index)));



   lfp1_filt = zscore(delta_this_pb_1);
   lfp2_filt = zscore(delta_this_pb_2);


 % Cross-correlation

 [xc,~] = xcorr(lfp1_filt, ...
     lfp2_filt, ...
     max_lag, ...
     'coeff');

    cross_corr(j,:) = xc;


end





%% estimating local (time convolved)


%%
cc_structure = struct();

% Cross-correlogram structure

cc_structure = struct();
% Metadata

cc_structure.time_range                = time_range;
cc_structure.cross_correlogram_length  = cross_correlogram_lenth;

cc_structure.sr                        = sr;
cc_structure.max_lag_sec               = max_lag_sec;
cc_structure.max_lag                   = max_lag;

cc_structure.play_behaviors            = play_behaviors;

% Filter info

cc_structure.filter                    = Hd_freq;

% Behavior/session info

cc_structure.Behavior                  = Behavior;
cc_structure.play_bouts_table          = play_bouts_table;
cc_structure.play_session              = play_session;

% Time axes

cc_structure.lfp_time                  = lfp_time;

cc_structure.cross_corr_lags_sec = ...
    (-max_lag:max_lag) ./ sr;

% Filtered signals

cc_structure.delta_animal_1            = delta_animal_1;
cc_structure.delta_animal_2            = delta_animal_2;

% Cross-correlation results

cc_structure.cross_corr                = cross_corr;

% Summary statistics

cc_structure.mean_cross_corr = ...
    mean(cross_corr,1,'omitnan');

cc_structure.sem_cross_corr = ...
    std(cross_corr,[],1,'omitnan') ./ ...
    sqrt(sum(~isnan(cross_corr),1));






end