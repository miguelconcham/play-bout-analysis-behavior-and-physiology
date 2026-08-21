function psth_struct = GENERATE_THETA_PSTH_CALLS(animal_code,f,freq_range, wind_params)
% GENERATE_THETA_PSTH
% Computes freq-band (6–12 Hz) LFP power peri-event time histograms (PSTHs)
% around social play bouts and key behaviors. Outputs a struct with aligned data.
%% define parameters

% Spectrogram parameters
hist_range      = [-2 2];       % peri-event window for onset/offset (s)

% Spectrogram parameters
wind_length     = wind_params(1);    % 250 ms
wind_overlap    = wind_params(2);   % 249 ms overlap
spect_bin_size  = wind_length-wind_overlap;

% f          = 4:.1:15;      % frequency range for spectrogram
% freq_range = [6 12];       % freq band



%% -------------------- SETUP AND PATHS --------------------
% Define behavior types of interest
play_behaviors = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};
area2analyze = 'PAG';
% Define data directories for synchronization, anatomical limits, and behavior
call_folder         =  '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\CallDetectionBackup';
synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\Area_limits_GoodLooking.xlsx';
npx_folder          = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\NPX raw data';

% Extract animal/session code information from current directory

animal_code_params = strsplit(animal_code, ' ');
animal_batch       = animal_code_params{1};
date               = animal_code_params{2};
repeated_animal    = animal_code_params{3};

%% -------------------- LOAD SYNCHRONIZATION MODEL --------------------
% Load mapping between video time and neural time
load([synch_directory,'\', animal_code, '\synch_model_video2NPX.mat'])
load([synch_directory,'\', animal_code, '\synch_model_audio2NPX.mat'])



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
    % Preprocessed: use ChannelMap.mat to locate mid-PAG channel
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
f_index = f>=freq_range(1) & f<=freq_range(2);

% Only compute spectrogram over LFP range covering all play bouts
range2exctract = LFP_time>=min(CallStats.BeginTimes)+hist_range(1) & LFP_time<=max(CallStats.EndTimes)+hist_range(2);

% Initialize output struct
psth_struct = [];

for ch_n=1:numel(mid_PAG_channel)
    disp('Estimating spectrogram')
    % Compute spectrogram for selected LFP channel
    [pow_spectrogram,~,spect_time]  = spectrogram(PAG_LFP(ch_n,range2exctract),wind_length*2500, floor(wind_overlap*2500), f,2500);
    spect_time = spect_time+min(CallStats.BeginTimes)+hist_range(1);
    pow_spectrogram = abs(pow_spectrogram);

    % Extract freq-band power (log-scaled, smoothed)
    freq_pow = mean(log10(pow_spectrogram(f_index,:)));
    freq_pow = movmean(freq_pow,1/max(freq_range));

    %% -------------------- COMPUTE PERI-EVENT THETA PSTHs --------------------
    % Allocate arrays for different peri-event alignments
    call_onset              = nan(size(CallStats,1), round(range(hist_range)/spect_bin_size) +1);
    call_offset             = nan(size(CallStats,1), round(range(hist_range)/spect_bin_size)+1);
    call_onset_regressor    = nan(size(CallStats,1), round(range(hist_range)/spect_bin_size)+1);
    call_offset_regressor   = nan(size(CallStats,1), round(range(hist_range)/spect_bin_size)+1);
  
    % Loop over play bouts
    for pb_n=1:size(CallStats,1)
        % Get bout start/end
        play_bout_start = CallStats.BeginTimes(pb_n);
        play_bout_end   = CallStats.EndTimes(pb_n);
        [~,loc_start] = min(abs(spect_time-play_bout_start));
        [~,loc_end]   = min(abs(spect_time-play_bout_end));
       

        % --- Align to bout onset ---
        entire_range = round(loc_start+(hist_range(1)/spect_bin_size)+1):(loc_start+(hist_range(2)/spect_bin_size));

        allowed_index_freq = ismember(entire_range,1:size(freq_pow,2));
        call_onset(pb_n,allowed_index_freq) = freq_pow(entire_range(allowed_index_freq));
        current_time = nan(1,size(call_onset,2));
        current_time(allowed_index_freq) = spect_time(entire_range(allowed_index_freq));
        call_onset_regressor(pb_n,:) =any(current_time>=CallStats.BeginTimes & current_time<=CallStats.EndTimes,1);

        % --- Align to bout offset ---
        entire_range = round(loc_end+(hist_range(1)/spect_bin_size):loc_end+(hist_range(2)/spect_bin_size));
        allowed_index_freq = ismember(entire_range,1:size(freq_pow,2));
        call_offset(pb_n,allowed_index_freq) = freq_pow(entire_range(allowed_index_freq));
        current_time = nan(1,size(call_offset,2));
        current_time(allowed_index_freq) = spect_time(entire_range(allowed_index_freq));
        call_offset_regressor(pb_n,:) =any(current_time>=CallStats.BeginTimes & current_time<=CallStats.EndTimes,1);
       

    end
     

    %% -------------------- STORE RESULTS IN STRUCT --------------------
    if ch_n==1
        % First channel: initialize struct
        psth_struct.call_onset              = call_onset;
        psth_struct.call_offset             = call_offset;
        psth_struct.call_onset_regressor    = call_onset_regressor;
        psth_struct.call_offset_regressor   = call_offset_regressor;
        psth_struct.CallStats               = CallStats;
        psth_struct.hist_range              = hist_range;
        psth_struct.wind_length             = wind_length;
        psth_struct.wind_overlap            = wind_overlap;
        psth_struct.ch                      = mid_PAG_channel(ch_n);
        psth_struct.f                       = f;
        psth_struct.freq_range              = freq_range;
    else
        % Additional channels: store as array of structs
        psth_struct(ch_n).call_onset              = call_onset;
        psth_struct(ch_n).call_offset             = call_offset;
        psth_struct(ch_n).call_onset_regressor    = call_onset_regressor;
        psth_struct(ch_n).call_offset_regressor   = call_offset_regressor;
        psth_struct(ch_n).CallStats               = CallStats;
        psth_struct(ch_n).hist_range              = hist_range;
        psth_struct(ch_n).wind_length             = wind_length;
        psth_struct(ch_n).wind_overlap            = wind_overlap;
        psth_struct(ch_n).ch                      = mid_PAG_channel(ch_n);
        psth_struct(ch_n).f                       = f;
        psth_struct(ch_n).freq_range              = freq_range;
     
    end
end
end
