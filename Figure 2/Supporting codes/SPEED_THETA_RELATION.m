function psth_struct = SPEED_THETA_RELATION(animal_code, pt,freq_range,f,wind_length,wind_overlap )
% GENERATE_THETA_PSTH
% Computes freq-band (6–12 Hz) LFP power peri-event time histograms (PSTHs)
% around social play bouts and key behaviors. Outputs a struct with aligned data.

varNames2store = { 'model_data','lm',...
     'wind_length', 'wind_overlap', ...
    'f', 'freq_range', 'mid_PAG_channel' ...
     'play_behaviors', 'Behavior' ...
     'play_bouts_table','lm'};
%% ---- DEFINE PARAMETERS FOR SPECTROGRAM

% Spectrogram parameters
hist_range      = [-20 20];       % peri-event window for onset/offset (s)






%% -------------------- SETUP AND PATHS --------------------
% Define behavior types for playbout 
play_behaviors = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};

% Define data directories for loading data
synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\Area_limits_GoodLooking.xlsx';
npx_folder          = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\NPX raw data';
trakcking_folder    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Traking backups';
behavior_data       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Behavior backups';

% Extract animal/session code information from animal_code

animal_code_params = strsplit(animal_code, ' ');
animal_batch       = animal_code_params{1};
date               = animal_code_params{2};
repeated_animal    = animal_code_params{3};

%% -------------------- LOAD SYNCHRONIZATION MODEL --------------------
% Load mapping between video time and neural time
load([synch_directory,'\', animal_code, '\synch_model_video2NPX.mat'])
load([synch_directory,'\', animal_code, '\synch_model_audio2NPX.mat'])

%% --------------------------- LOAD TRAKING AND ANIMAL VARIABLES for partner = pt --------


load([trakcking_folder,'\',animal_code,' P', num2str(pt), ' traking_structure.mat']) 


traking_time        = predict(synch_model_video2NPX, (traking_structure.frames2stract/30)')';
full_traking_time   = (traking_structure.frames2stract(1):traking_structure.frames2stract(end))/30;
full_traking_time   = predict(synch_model_video2NPX, full_traking_time')';
animal_pos          = nan(numel(full_traking_time),2);
partner_pos         = nan(numel(full_traking_time),2);
animal_pos(:,1)     = interp1(traking_time, smoothdata(traking_structure.animal_pos(:,1), 'loess',5), full_traking_time,'cubic');
animal_pos(:,2)     = interp1(traking_time,smoothdata(traking_structure.animal_pos(:,2), 'loess',5), full_traking_time,'cubic');


partner_pos(:,1)    = interp1(traking_time, smoothdata(traking_structure.partner_pos(:,1), 'loess',5), full_traking_time,'cubic');
partner_pos(:,2)    = interp1(traking_time,smoothdata(traking_structure.partner_pos(:,2), 'loess',5), full_traking_time,'cubic');


animal_dist         = sqrt(sum((partner_pos-animal_pos).*(partner_pos-animal_pos),2));
animal_speed        = sqrt(sum(diff(animal_pos).*diff(animal_pos),2));
animal_accel        = abs(diff(animal_speed));
animal_rel_speed    = diff(animal_dist);

%% ----------- LOAD BEHAVIOR DATA AND OBTAIN PLAYBOUTS


Behavior_file =[behavior_data,'\', animal_code,'.txt'];%load behavior data

Behavior                            = readtable(Behavior_file);
Behavior(:,2)                       = [];
Behavior.Properties.VariableNames   = {'Animal', 'Start', 'End', 'Length', 'Type'};


bin_size                = 0.01;
conv_length             = 1;
Behavior.Type2          = Behavior.Type;

Behavior.Type2(ismember(Behavior.Type2, {'Pounce_A', 'Pounce_B'}))      = {'Pounce'}; %% Merging behaviors to Type2
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'}))    = {'PounceI'};
Behavior.Type2(ismember( Behavior.Type2,''))                            = {'Other'};
Behavior(ismember(Behavior.Animal, 'Reversal'),:)                       = [];
Behavior(ismember(Behavior.Animal, 'Session_structure'),:)              = [];

animal_types                                                            = unique(Behavior.Animal);
disp('Animal this session')
disp(animal_types')


Behavior.Start          = predict(synch_model_video2NPX, Behavior.Start);
Behavior.End            = predict(synch_model_video2NPX, Behavior.End);

config.Behavior         = Behavior;
config.repeated_animal  = repeated_animal;
config.animal_types     = animal_types        ;
config.play_behaviors   = play_behaviors      ;
config.beh_bin          = bin_size             ;
config.conv_length      = conv_length;
config.behavior_window  = 0;


[play_bouts_table]      = play_bout(config);
%% -------------------- LOAD LFP DATA --------------------
disp('LOADING LFP')
if exist([npx_folder,'\',animal_code,'\LFP_PAG.mat'], 'file')==2
    % Preprocessed LFP file exists
    NPX_Type = 2;
    load([npx_folder,'\',animal_code,'\LFP_PAG.mat'], 'LFP')
elseif exist([npx_folder,'\',animal_code,'\LFP_PAG.dat'], 'file')==2
    % Load raw binary LFP file
    NPX_Type = 1;
    file_pointer = fopen([npx_folder,'\',animal_code,'\LFP_PAG.dat'], 'r');
    LFP = fread(file_pointer,'int16');
    LFP = reshape(LFP, 384, numel(LFP)/384);
end
disp('LFP LOADED')


%% -------------------- SELECT PAG CHANNEL(S) --------------------
disp('Loading Channel Map')
hard_coded_x_coords = [8 40;258 290; 508 540; 758 790];
area_limit = readtable(area_limit_table);

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

%% -------------------- COMPUTE SPECTROGRAM & THETA POWER --------------------

% Extract LFP from mid-PAG channel(s)
PAG_LFP = double(LFP(mid_PAG_channel,:));
clear LFP
LFP_time = (1:size(PAG_LFP,2))/2500; % sample rate: 2.5 kHz

% Identify freq indices in spectrogram
f_index = f>=freq_range(1) & f<=freq_range(2);

% Only compute spectrogram over LFP range covering all play bouts
range2exctract      = LFP_time>=min(traking_time)+hist_range(1) & LFP_time<=max(traking_time)+hist_range(2);


% Initialize output struct
psth_struct = [];

for ch_n=1:numel(mid_PAG_channel)
    disp('Estimating spectrogram')
    % Compute spectrogram for selected LFP channel
    [pow_spectrogram,~,spect_time]  = spectrogram(PAG_LFP(ch_n,range2exctract),wind_length*2500, floor(wind_overlap*2500), f,2500);
    spect_time = spect_time+min(traking_time)+hist_range(1) - .5*wind_length;
    
        
    disp('Estimating Onset and Offset PSTHs')
    animal_speed_spect_time     = interp1(full_traking_time(1:end-1)', animal_speed,spect_time');
    animal_speed_spect_time = (animal_speed_spect_time- mean(animal_speed_spect_time, 'omitmissing'))/std(animal_speed_spect_time, 'omitmissing');
    animal_accel_spect_time     = interp1(full_traking_time(2:end-1)',animal_accel,spect_time');
    animal_accel_spect_time = (animal_accel_spect_time- mean(animal_accel_spect_time, 'omitmissing'))/std(animal_accel_spect_time, 'omitmissing'); 

    animal_dist_spect_time    = interp1(full_traking_time,animal_dist,spect_time');
    animal_dist_spect_time = (animal_dist_spect_time- mean(animal_dist_spect_time, 'omitmissing'))/std(animal_dist_spect_time, 'omitmissing'); 

    animal_rel_speed_spect_time    = interp1(full_traking_time(1:end-1)',animal_rel_speed,spect_time');
    animal_rel_speed_spect_time = (animal_dist_spect_time- mean(animal_rel_speed_spect_time, 'omitmissing'))/std(animal_rel_speed_spect_time, 'omitmissing'); 

    play_bout_index = any(spect_time>=play_bouts_table(:,1) & spect_time<=play_bouts_table(:,2));
    
    pow_spectrogram = abs(pow_spectrogram);

    % Extract freq-band power (log-scaled, smoothed)
    freq_pow = mean(log10(pow_spectrogram(f_index,:)));
    freq_pow = movmean(freq_pow,1/max(freq_range));
    
    model_data = table(animal_speed_spect_time, animal_accel_spect_time, categorical(play_bout_index'),animal_dist_spect_time,animal_rel_speed_spect_time, freq_pow', 'VariableNames',{'Speed','Acc','Play','Distance','RelativeSpeed', 'Power'} );
    lm = fitlm(model_data, 'Power ~ Speed + Acc + Distance + RelativeSpeed  +Play');
    %% -------------------- STORE RESULTS IN STRUCT --------------------
    if ch_n==1
        % First channel: initialize struct
        psth_struct = struct();
        for i = 1:numel(varNames2store)
            name = varNames2store{i};
            psth_struct.(name) = eval(name);
        end
       
    else
        % Additional channels: store as array of structs
       for i = 1:numel(varNames2store)
            name = varNames2store{i};
            psth_struct(ch_n).(name) = eval(name);
        end
     
    end
end
end
