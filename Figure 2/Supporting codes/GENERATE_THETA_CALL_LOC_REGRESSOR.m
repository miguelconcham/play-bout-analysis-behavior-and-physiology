function psth_struct = GENERATE_THETA_CALL_LOC_REGRESSOR(animal_code, pt, f, freq_range, wind_length,spect_bin_size )
% GENERATE_THETA_PSTH
% Computes freq-band (6–12 Hz) LFP power peri-event time histograms (PSTHs)
% around social play bouts and key behaviors. Outputs a struct with aligned data.
% kalman bool (if true will apply kalman estiamtes for speed and
% acceleration)

do_kalman       = false;
area2analyze    = 'PAG';
varNames2store = { 'play_bout_onset_low_freq','play_bout_onset_all_low_freq',...
'play_bout_onset_high_freq','play_bout_onset_all_high_freq',...  
    'regressors','variable_names', ...
    'hist_range', 'wind_length','spect_bin_size', ...
    'f', 'freq_range', 'mid_PAG_channel' ...
     'play_behaviors', 'CallStats', 'Behavior' ...
     'play_bouts_table'};

variable_names = {
    'animal_angle', 'animal_angle_speed', 'animal_angle_acc', ...
    'animal_speed', 'animal_accel', 'animal_speed_kalman', 'animal_accel_kalman',...
    'partner_angle', 'partner_angle_speed', 'partner_angle_acc', ...
    'partner_speed', 'partner_accel', 'partner_speed_kalman', 'partner_accel_kalman', ...
    'relative_distance',  'relative_angle', ...
    'relative_angle_speed', 'relative_angle_acc', 'relative_speed', 'relative_acc', ...
    'relative_acc_kalman', 'relative_speed_kalman'};

% Define groups for interpolation rules
if do_kalman
    group1 = {'animal_speed', 'partner_speed', 'relative_speed', 'relative_angle_speed', 'animal_angle_speed', 'partner_angle_speed'};
    group2 = {'animal_accel', 'animal_angle_acc', 'partner_accel', 'partner_angle_acc', 'relative_acc', 'relative_angle_acc'};
else
    group1 = {'animal_speed', 'partner_speed', 'relative_speed', 'relative_angle_speed', 'animal_angle_speed', 'partner_angle_speed', ...
        'relative_speed_kalman','partner_speed_kalman','animal_speed_kalman'};
    group2 = {'animal_accel', 'animal_angle_acc', 'partner_accel', 'partner_angle_acc', 'relative_acc', 'relative_angle_acc', ...
        'relative_acc_kalman','animal_accel_kalman','partner_accel_kalman'};
end

% Parameters for Kalman Filter 
% (read demo in ...Codes\Kalman Filter Package\Examples\KF_Demo2.m)
if do_kalman
    MAlen = 30; % moving average filter length (1 second in our case)
    N = 15; % lookback window length for covariance estimation (half a second in our case)
end
%% ---- DEFINE PARAMETERS FOR SPECTROGRAM

% Spectrogram parameters
hist_range      = [-20 20];       % peri-event window for onset/offset (s)
low_f           = f{1};
low_freq_range  = freq_range{1};
high_f          = f{2};
high_freq_range = freq_range{2};

% Spectrogram parameters
low_wind_length     = wind_length{1};    % 250 ms
low_wind_overlap    = low_wind_length-spect_bin_size;   % 249 ms overlap


high_wind_length     = wind_length{2};    % 250 ms
high_wind_overlap    = high_wind_length-spect_bin_size;





%% -------------------- SETUP AND PATHS --------------------
% Define behavior types for playbout 
play_behaviors = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};

% Define data directories for loading data
call_folder         =  '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\CallDetectionBackup';
synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Synch data';
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Area_limits_GoodLooking.xlsx';
npx_folder          = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
trakcking_folder    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Traking backups';
behavior_data       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Behavior backups';

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
full_tracking_time   = (traking_structure.frames2stract(1):traking_structure.frames2stract(end))/30;
full_tracking_time   = predict(synch_model_video2NPX, full_tracking_time')';

animal_pos          = nan(numel(full_tracking_time),2);
partner_pos         = nan(numel(full_tracking_time),2);
animal_pos(:,1)     = interp1(traking_time, smoothdata(traking_structure.animal_pos(:,1), 'loess',5), full_tracking_time,'cubic');
animal_pos(:,2)     = interp1(traking_time,smoothdata(traking_structure.animal_pos(:,2), 'loess',5), full_tracking_time,'cubic');
animal_velocity     = diff(animal_pos);
animal_angle        = cart2pol(animal_pos(:,1), animal_pos(:,2));
animal_angle_speed  = angdiff(animal_angle);
animal_angle_acc    = angdiff(animal_angle_speed);

animal_speed        = sqrt(sum(diff(animal_pos).*diff(animal_pos),2));
animal_accel        = abs(diff(animal_speed));

% Perform Square Root Kalman filtering


if do_kalman
    start_pos           = animal_pos(1,:);
    z                   = animal_pos' - start_pos';
    [animal_kf_chat, ~]       = SquareRootKalmanFilter_CA(z,MAlen,N);
    animal_speed_kalman = sqrt(sum(animal_kf_chat([3 4],:).*animal_kf_chat([3 4],:)));
    animal_accel_kalman = sqrt(sum(animal_kf_chat([5 6],:).*animal_kf_chat([5 6],:)));
else
    animal_speed_kalman = animal_speed*NaN;
    animal_accel_kalman = animal_accel*NaN;
end

   




partner_pos(:,1)    = interp1(traking_time, smoothdata(traking_structure.partner_pos(:,1), 'loess',5), full_tracking_time,'cubic');
partner_pos(:,2)    = interp1(traking_time,smoothdata(traking_structure.partner_pos(:,2), 'loess',5), full_tracking_time,'cubic');
partner_velocity     = diff(partner_pos);
partner_angle       = cart2pol(partner_pos(:,1), partner_pos(:,2));
partner_angle_speed = angdiff(partner_angle);
partner_angle_acc   = angdiff(partner_angle_speed);

partner_speed       = sqrt(sum(diff(partner_pos).*diff(partner_pos),2));
partner_accel       = abs(diff(partner_speed));


% Perform Square Root Kalman filtering
if do_kalman
    start_pos               = partner_pos(1,:);
    z                       = partner_pos' - start_pos';
    [partner_kf_chat, ~]          = SquareRootKalmanFilter_CA(z,MAlen,N);
    partner_speed_kalman    = sqrt(sum(partner_kf_chat([3 4],:).*partner_kf_chat([3 4],:)));
    partner_accel_kalman    = sqrt(sum(partner_kf_chat([5 6],:).*partner_kf_chat([5 6],:)));

else
    partner_speed_kalman = partner_speed*NaN;
    partner_accel_kalman = partner_accel*NaN;

end


relative_distance       = sqrt(sum((partner_pos-animal_pos).*(partner_pos-animal_pos),2));
relative_pos            = animal_pos-partner_pos;
relative_velocity       = animal_velocity-partner_velocity;
relative_distance       = sqrt(sum(relative_pos.*relative_pos,2));
relative_angle          = cart2pol(relative_pos(:,1), relative_pos(:,2));
relative_angle_speed    =  cart2pol(relative_velocity(:,1), relative_velocity(:,2));
relative_angle_acc      =  angdiff(relative_angle_speed);

relative_speed          = diff(relative_distance);
relative_acc            = diff(relative_speed);
relative_speed          = sqrt(sum(relative_speed.*relative_speed,2));
relative_acc            = sqrt(sum(relative_acc.*relative_acc,2));
if do_kalman
    relative_vel_kalman     = partner_kf_chat([3 4],:) - animal_kf_chat([3 4],:);
    relative_acc_kalman     = partner_kf_chat([5 6],:)  - animal_kf_chat([5 6],:);
    relative_speed_kalman   = sqrt(sum(relative_vel_kalman.*relative_vel_kalman));
    relative_acc_kalman     = sqrt(sum(relative_acc_kalman.*relative_acc_kalman));

else
    relative_speed_kalman = partner_speed_kalman*NaN;
    relative_acc_kalman = partner_accel_kalman*NaN;
end





vars_raw = struct();
for iVar = 1:length(variable_names)
    varName = variable_names{iVar};
    % Assume these variables are in workspace
    vars_raw.(varName) = eval(varName);
end



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

%% Load calls
CallStats = readtable([call_folder, '\', animal_code, '_Stats.xlsx']) ;
CallStats.Properties.VariableNames = cellfun(@(x) strrep(x, '_', ''),CallStats.Properties.VariableNames, 'UniformOutput',false );
CallStats.BeginTimes = predict(synch_model_audio2NPX,CallStats.BeginTimes);
CallStats.EndTimes = predict(synch_model_audio2NPX,CallStats.EndTimes);
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

load([npx_folder,'\',animal_code,'\chann_map_',area2analyze,'.mat'], 'xcoords', 'ycoords','chanMap')
% Build animal identifier for area selection
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal,animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName,this_animal),:);
 figure('units','normalized','outerposition',[0 0 .2 1]);

if NPX_Type == 1

    if ~ismember(192, chanMap)

        pos_191 = find(chanMap==191);
        pos_193 = find(chanMap==193);

        if pos_193 == pos_191+1

            xcoords = [xcoords;NaN];
            xcoords(pos_193+1:end) = xcoords(pos_193:end-1);
            xcoords(pos_193) = 43;
            ycoords = [ycoords;NaN];
            ycoords(pos_193+1:end) = ycoords(pos_193:end-1);
            ycoords(pos_193) = 1900;
            chanMap = [chanMap;NaN];
            chanMap(pos_193+1:end) = chanMap(pos_193:end-1);
            chanMap(pos_193) = 192;
        else
            disp('Inconsistent ChannelMap')
            return
        end
    end
    plot(xcoords,ycoords, 'k.'); hold on
    % Raw LFP: select channel range for LPAG region

    Y_Range = area_limit{ismember(area_limit.area, {'LPAG'}), {'ProbeNum','depth_start', 'depth_end'}};
    this_indexes = ycoords>=Y_Range(2) & ycoords<=Y_Range(3);
    all_locs = [xcoords(this_indexes) ycoords(this_indexes)];
    plot(all_locs(:,1),all_locs(:,2), 'r.')
    mean_loc = mean(all_locs);
    [~, closest_channel]= min(sum(abs([xcoords ycoords]-repmat(mean_loc,numel(ycoords),1)),2));
    plot(xcoords(closest_channel), ycoords(closest_channel), 'xb')
    mid_PAG_channel = chanMap(closest_channel);
    title([this_animal, ' Probe#', num2str(Y_Range(1))])
else
    % Preprocessed: use chann_map_PAG.mat to locate mid-PAG channel
    plot(xcoords,ycoords, 'k.'); hold on
    Y_Range = area_limit{ismember(area_limit.area, {'LPAG'}), {'ProbeNum','depth_start', 'depth_end'}};

    mid_PAG_channel = nan(size(Y_Range,1),1);

    for j=1:size(Y_Range,1)
        this_indexes = ycoords>=Y_Range(j,2) & ycoords<=Y_Range(j,3) & ismember(xcoords,hard_coded_x_coords(Y_Range(j,1),:));
        all_locs = [xcoords(this_indexes) ycoords(this_indexes)];
        plot(all_locs(:,1),all_locs(:,2), 'r.')
        mean_loc = mean(all_locs);
        [~, closest_channel]= min(sum(abs([xcoords ycoords]-repmat(mean_loc,numel(ycoords),1)),2));
        plot(xcoords(closest_channel), ycoords(closest_channel), 'xb')
        mid_PAG_channel(j) = chanMap(closest_channel);
        title([this_animal, ' Probe#', num2str(Y_Range(1,:))])
    end
end

pause(.1)
%% -------------------- COMPUTE SPECTROGRAM & THETA POWER --------------------

% Extract LFP from mid-PAG channel(s)
PAG_LFP = double(LFP(mid_PAG_channel,:));
clear LFP
LFP_time = (1:size(PAG_LFP,2))/2500; % sample rate: 2.5 kHz

% Identify freq indices in spectrogram
f_index_low     = low_f>=low_freq_range(1) & low_f<=low_freq_range(2);
f_index_high    = high_f>=high_freq_range(1) & high_f<=high_freq_range(2);
% Only compute spectrogram over LFP range covering all play bouts
range2exctract      = LFP_time>=min(traking_time)+hist_range(1) & LFP_time<=max(traking_time)+hist_range(2);
play_bouts_table    = play_bouts_table(play_bouts_table(:,1)>=min(traking_time) & play_bouts_table(:,2)<=max(traking_time),:);

slef_onset_offset = [Behavior.Start(ismember(Behavior.Animal,repeated_animal)) Behavior.End(ismember(Behavior.Animal,repeated_animal))];
other_onset_offset = [Behavior.Start(~ismember(Behavior.Animal,repeated_animal)) Behavior.End(~ismember(Behavior.Animal,repeated_animal))];

% Initialize output struct
psth_struct = [];

for ch_n=1:numel(mid_PAG_channel)
    disp('Estimating spectrogram')
    % Compute spectrogram for selected LFP channel
    [pow_spectrogram_low_freq,~,spect_time]  = spectrogram(PAG_LFP(ch_n,range2exctract),low_wind_length*2500, floor(low_wind_overlap*2500), low_f,2500);
    spect_time = spect_time+min(traking_time)+hist_range(1) - .5*low_wind_length;


    [pow_spectrogram_high_freq,~,spect_time_high_freq]  = spectrogram(PAG_LFP(ch_n,range2exctract),high_wind_length*2500, floor(high_wind_overlap*2500), high_f,2500);
    pow_spectrogram_high_freq = abs(pow_spectrogram_high_freq);
    pow_spectrogram_high_freq_interp = nan(size(pow_spectrogram_high_freq,1),numel(spect_time));
    spect_time_high_freq = spect_time_high_freq+min(traking_time)+hist_range(1) - .5*high_wind_length;

    for fn=1:size(pow_spectrogram_high_freq,1)
        pow_spectrogram_high_freq_interp(fn,:) = interp1(spect_time_high_freq,pow_spectrogram_high_freq(fn,:),spect_time);
    end
    pow_spectrogram_high_freq = pow_spectrogram_high_freq_interp;
    clear pow_spectrogram_high_freq_interp
    real_bin_size = mean(diff(spect_time));




    disp('Estimating Onset and Offset PSTHs')
    
    vars_spect_time = struct();

    for iVar = 1:length(variable_names)
        varName = variable_names{iVar};
        rawVar = vars_raw.(varName);

        if ismember(varName, group1)
            % group 1: interp1(full_tracking_time(1:end-1)', VARIABLE, spect_time')
            vars_spect_time.(varName) = interp1(full_tracking_time(1:end-1)', rawVar, spect_time');

        elseif ismember(varName, group2)
            % group 2: interp1(full_tracking_time(2:end-1)', VARIABLE, spect_time')
            vars_spect_time.(varName) = interp1(full_tracking_time(2:end-1)', rawVar, spect_time');

        elseif strcmp(varName, 'partner_pos') || strcmp(varName, 'relative_pos')
            % Special cases for 2D pos?
            % If your pos variables are 2D matrices, interpolate each dimension separately
            % Assuming rawVar is Nx2 and full_tracking_time is Nx1
            vars_spect_time.(varName)(:,1) = interp1(full_tracking_time, rawVar(:,1), spect_time');
            vars_spect_time.(varName)(:,2) = interp1(full_tracking_time, rawVar(:,2), spect_time');

        else
            % default: interp1(full_tracking_time, VARIABLE, spect_time')
            vars_spect_time.(varName) = interp1(full_tracking_time, rawVar, spect_time');
        end
    end


    pow_spectrogram_low_freq = abs(pow_spectrogram_low_freq);
   
    % Extract freq-band power (log-scaled, smoothed)
    freq_pow_low = mean(log10(pow_spectrogram_low_freq(f_index_low,:)));
    freq_pow_low = movmean(freq_pow_low,1/(max(low_freq_range)*real_bin_size));

    freq_pow_high = mean(log10(pow_spectrogram_high_freq(f_index_high,:)));
    freq_pow_high = movmean(freq_pow_high,1/(max(high_freq_range)*real_bin_size));

    %% -------------------- COMPUTE PERI-EVENT psth power by frequency and behaviors --------------------
    % Allocate arrays for different peri-event alignments
   
    play_bout_onset_low_freq                = nan(size(play_bouts_table,1), round(range(hist_range)/spect_bin_size));
    play_bout_onset_all_low_freq            = nan(size(play_bouts_table,1),numel(low_f), round(range(hist_range)/spect_bin_size));

    play_bout_onset_high_freq                = nan(size(play_bouts_table,1), round(range(hist_range)/spect_bin_size));
    play_bout_onset_all_high_freq            = nan(size(play_bouts_table,1),numel(high_f), round(range(hist_range)/spect_bin_size));

    %Create all needed regressors (play bout, call, speed, acceleration,
    %self and other)
    n_bins = round(range(hist_range) / spect_bin_size); % number of bins in PSTH
    n_playbouts = size(play_bouts_table, 1);

    % Pre-allocate regressors container (structure)
    regressors = struct();

    % Pre-allocate all regressors to NaNs (play_bout and call handled separately)
    for iVar = 1:length(variable_names)
        varName = variable_names{iVar};
        regressors.([varName '_onset']) = nan(n_playbouts, n_bins);
      
    end
    play_bout_onset_regressor = nan(n_playbouts, n_bins);
    call_onset_regressor = nan(n_playbouts, n_bins);
    self_onset_regressor               = nan(n_playbouts, n_bins);
    other_onset_regressor  = nan(n_playbouts, n_bins);
    % Loop over play bouts
    for pb_n=1:size(play_bouts_table,1)
        % Get bout start/end
        play_bout_start     = play_bouts_table(pb_n,1);
        play_bout_end       = play_bouts_table(pb_n,2);
        [~,loc_start]       = min(abs(spect_time-play_bout_start));
        [~,loc_end]         = min(abs(spect_time-play_bout_end));       

        % --- Align to bout onset ---
        entire_range_onset                                      = round(loc_start+(hist_range(1)/spect_bin_size)):(loc_start+(hist_range(2)/spect_bin_size));
        if numel(entire_range_onset)>round(range(hist_range)/spect_bin_size)
            entire_range_onset(1) = [];
        end
        allowed_index_freq                                      = ismember(entire_range_onset,1:size(freq_pow_low,2));
        play_bout_onset_low_freq(pb_n,allowed_index_freq)       = freq_pow_low(entire_range_onset(allowed_index_freq));
        play_bout_onset_all_low_freq(pb_n,:,allowed_index_freq)  = pow_spectrogram_low_freq(:,entire_range_onset(allowed_index_freq));
        play_bout_onset_high_freq(pb_n,allowed_index_freq)       = freq_pow_high(entire_range_onset(allowed_index_freq));
        play_bout_onset_all_high_freq(pb_n,:,allowed_index_freq)  = pow_spectrogram_high_freq(:,entire_range_onset(allowed_index_freq));

        current_time_onset                                      = nan(1,size(call_onset_regressor,2));
        current_time_onset(allowed_index_freq)                  = spect_time(entire_range_onset(allowed_index_freq));
        play_bout_onset_regressor(pb_n,:)                       = any(current_time_onset>=play_bouts_table(:,1) & current_time_onset<=play_bouts_table(:,2),1);
        call_onset_regressor(pb_n,:)                            = any(current_time_onset>=CallStats.BeginTimes & current_time_onset<=CallStats.EndTimes,1);
        

                         
        
        self_onset_regressor(pb_n,:)                            = any(current_time_onset>=slef_onset_offset(:,1) & current_time_onset<=slef_onset_offset(:,2),1);
        other_onset_regressor(pb_n,:)                           = any(current_time_onset>=other_onset_offset(:,1) & current_time_onset<=other_onset_offset(:,2),1);
        for iVar = 1:length(variable_names)
            varName = variable_names{iVar};
            spect_var = vars_spect_time.(varName);

            % Assign onset regressor for this play bout and bins
            % Example: regressors.animal_speed_onset(pb_n, allowed_index_freq) = animal_speed_spect_time(entire_range_onset(allowed_index_freq));
            regressors.([varName '_onset'])(pb_n, allowed_index_freq) = spect_var(entire_range_onset(allowed_index_freq));

            % Similarly for offset if you have offset indices available:
            % regressors.([varName '_offset'])(pb_n, allowed_index_freq) = spect_var(entire_range_offset(allowed_index_freq));
        end
    
    end
     
    
    regressors.play_bout_onset_regressor    = play_bout_onset_regressor;
    regressors.call_onset_regressor         = call_onset_regressor; 
    regressors.self_onset_regressor         = self_onset_regressor ;  
    regressors.other_onset_regressor        = other_onset_regressor;

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
