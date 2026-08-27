function mi_structure = GENERATE_STRUCUTRE_animal_synch(directory)
%% set parameters

%time range for the delta psth
time_range = [-5 10];
%lfp sampling rate
sr = 2500;
%params mutual information
n_bins_local        = 10;
n_bins_global       = 50;
t_size              = 200;
nIter               = 10000;  % number of randomizations
%params playbout
% non_play_behaviors   = {'PounceI',  'Bite'};
play_behaviors      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};
% play_behaviors      = non_play_behaviors;
bin_size            = 0.01;
conv_length         = 1;

%params spectrogram
wind_length         = 1;
wind_overlap        = .990;
f                   = .1:.1:6;
freq_pow_range      = [1 4];
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


index_play_session_range = round((play_session(1)-5)*sr:(play_session(2)+5)*sr);

[pow_spectrogram,~,spect_time]  = spectrogram(lfp_animal_1(index_play_session_range),wind_length*2500, wind_overlap*2500, f,2500);
spect_time                     = spect_time+play_session(1)-5;
pow_spectrogram = abs(pow_spectrogram);
mean_pow_animal_1 = mean(pow_spectrogram(f >=freq_pow_range(1) & f<=freq_pow_range(2),:));


[pow_spectrogram,~,~]  = spectrogram(lfp_animal_2(index_play_session_range),wind_length*2500, wind_overlap*2500, f,2500);
pow_spectrogram = abs(pow_spectrogram);
mean_pow_animal_2 = mean(pow_spectrogram(f >=freq_pow_range(1) & f<=freq_pow_range(2),:));
% play_song([],[],[])


 %%

 spect_sr = round(1/(wind_length-wind_overlap));

psth_time =( time_range(1)*spect_sr:time_range(2)*spect_sr)/spect_sr;
baseline = psth_time<-1;
n_samples = round((range(time_range)*spect_sr) + 1);
psth_onset = nan(2,size(play_bouts_table, 1),n_samples);
psth_offset = nan(2,size(play_bouts_table, 1),n_samples);

psth_zscored = nan(2,size(play_bouts_table, 1),n_samples);
for j=1:size(play_bouts_table, 1)

    beg_time = play_bouts_table(j,1);

    [~, loc] = min(abs(spect_time-beg_time));

    index2exctact = round((loc+time_range(1)*spect_sr):(loc+time_range(2)*spect_sr));
    possible_index = ismember(index2exctact,1:numel(spect_time));
    psth_onset(1,j,possible_index) = (mean_pow_animal_1(index2exctact(possible_index)));
    psth_onset(2,j,possible_index) = (mean_pow_animal_2(index2exctact(possible_index)));

    psth_zscored(1,j,:)  = ( psth_onset(1,j,:) - mean( psth_onset(1,j,baseline), 'omitmissing'))/ std( psth_onset(1,j,baseline), 'omitmissing');
    psth_zscored(2,j,:)  = ( psth_onset(2,j,:) - mean( psth_onset(2,j,baseline), 'omitmissing'))/ std( psth_onset(2,j,baseline), 'omitmissing');

    end_time = play_bouts_table(j,2);

    [~, loc] = min(abs(spect_time-end_time));

    index2exctact   = round((loc+time_range(1)*spect_sr):(loc+time_range(2)*spect_sr));
    possible_index  = ismember(index2exctact,1:numel(spect_time));
    psth_offset(1,j,possible_index) = (mean_pow_animal_1(index2exctact(possible_index)));
    psth_offset(2,j,possible_index) = (mean_pow_animal_2(index2exctact(possible_index)));

end


%% estimating global MI to psth onset (zscored)
global_time_range = [-2 4];
time4global = psth_time>global_time_range(1)  & psth_time<global_time_range(2);

I1 = squeeze(psth_zscored(1,:,time4global));
I2 = squeeze(psth_zscored(2,:,time4global));

mi_global = image_mutual_info(I1, I2, 50);

% generate_rand_distribution
MI_rand_global = nan(nIter,1);
for k = 1:nIter
    I2_shuffled = I2(randperm(size(I2,1)),:);
    MI_rand_global(k) = image_mutual_info(I1,   I2_shuffled, n_bins_global);
end

I1 = squeeze(psth_zscored(1,:,:));
I2 = squeeze(psth_zscored(2,:,:));

%% estimating local (time convolved)
[n,T] = size(I1);
mi_time             = psth_time(1:T-t_size) + .5*t_size/spect_sr;
play_onset = round(find(mi_time==0)- t_size/2);
% generate_rand_distribution
MI_t_rand = zeros(1,nIter);
for k = 1:nIter
    I2_shuffled = I2(randperm(size(I2,1)),(1:t_size));
    MI_t_rand(k) = image_mutual_info(I1(:,(1:t_size)),   I2_shuffled, n_bins_local);
end

MI_t_rand_play  = zeros(1,nIter);
for k = 1:nIter
    I2_shuffled = I2(randperm(size(I2,1)),(1:t_size)+play_onset);
    MI_t_rand_play(k) = image_mutual_info(I1(:,(1:t_size)+play_onset),   I2_shuffled, n_bins_local);
end




mi_t = nan(T-t_size,1);
mi_t_pctl =mi_t;
mi_t_by_trial = nan(n,T-t_size,1);
for t=1:T-t_size

    index = t:(t+t_size-1);
    mi_t(t) = image_mutual_info(I1(:,index), I2(:,index), n_bins_local);
    mi_t_pctl(t) = sum(MI_t_rand>mi_t(t))/nIter;

    for trial=1:n
        if ~(all(isnan(I1(trial,index))) ||  all(isnan(I2(trial,index))))
            mi_t_by_trial(trial,t) =  image_mutual_info(I1(trial,index), I2(trial,index), n_bins_local);
        end
    end
end

%% estimate mi properties per play bout
pb_lengths          = diff(play_bouts_table')';

mi_time             = psth_time(1:T-t_size) + .5*t_size/spect_sr;
baseline_index      = mi_time<0;
mi_t_by_trial_bc    = mi_t_by_trial;
mean_mi             = nan(n,1);
max_mi              = nan(n,1);
before_after_mi     = nan(n,2);
start_index         = mi_time>-1 & mi_time<0;

for  trial=1:n
        this_mi =   mi_t_by_trial(trial,:) ;


        mean_mi(trial)              = mean(mi_t_by_trial(trial,mi_time>0 & mi_time<pb_lengths(trial)));
        max_mi(trial)               = max(mi_t_by_trial(trial,mi_time>0 & mi_time<pb_lengths(trial)));
        mi_t_by_trial_bc(trial,:)   = (this_mi -mean(this_mi(baseline_index)))/std(this_mi(baseline_index));

        before_after_mi(trial,1)    = mean(mi_t_by_trial(trial,start_index))-mean_mi(trial);
        end_index                   = mi_time-pb_lengths(trial)>0 & mi_time-pb_lengths(trial)<1;
        before_after_mi(trial,2)    = mean(mi_t_by_trial(trial,end_index), 'omitmissing')-mean_mi(trial);


end


%%
mi_structure = [];
mi_structure.mi_global          = mi_global;
mi_structure.MI_rand_global     = MI_rand_global;
mi_structure.mi_t_pctl          = mi_t_pctl;
mi_structure.MI_t_rand          = MI_t_rand;
mi_structure.mi_t               = mi_t;
mi_structure.mi_t_by_trial      = mi_t_by_trial;
mi_structure.mi_t_by_trial_bc   = mi_t_by_trial_bc;
mi_structure.mean_mi            = mean_mi;
mi_structure.max_mi             = max_mi;
mi_structure.psth_time          = psth_time;
mi_structure.psth_zscored       = psth_zscored;
mi_structure.psth_onset         = psth_onset;
mi_structure.psth_offset        = psth_offset;
mi_structure.play_bouts_table   = play_bouts_table;
mi_structure.t_size             = t_size;
mi_structure.spect_sr           = spect_sr;
mi_structure.mi_time            = mi_time;
mi_structure.before_after_mi    = before_after_mi;
mi_structure.MI_t_rand_play     = MI_t_rand_play;
mi_structure.global_time_range  = global_time_range;
mi_structure.play_behaviors     = play_behaviors;
mi_structure.Behavior           = Behavior;







end