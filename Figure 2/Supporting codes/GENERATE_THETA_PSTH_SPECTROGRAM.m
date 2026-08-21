function psth_struct = GENERATE_THETA_PSTH_SPECTROGRAM(current_dir,wind_length,wind_overlap,min_separation,f,freq_pow_range )

play_behaviors      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};
% play_behaviors      = {'PounceI', 'Bite'};
% non_play_behaviors      = {'Grooming', 'PounceI','Rearing', 'Sniffing','Scratching', 'Bite'};
% exploratory_behaviors   = {'Grooming','Rearing', 'Sniffing','Scratching'};
% play_behaviors          = exploratory_behaviors;
synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\Area_limits_GoodLooking.xlsx';
behavior_data       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Behavior backups';
% npx_raw_data = 
animal_code         = strsplit(current_dir, '\');
animal_code         = animal_code{end};
area2analyze        = 'PAG';
animal_code_params  = strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
date                = animal_code_params{2};
repeated_animal     = animal_code_params{3};
%% define parameters



%% load synch from synch folder
load([synch_directory,'\', animal_code, '\synch_model_video2NPX.mat'])

%% 2 Load hmm data from HMM folder


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

animal_types            = unique(Behavior.Animal);

animal_types(ismember(animal_types,'Session_structure'))                =[];

Behavior.Start          = predict(synch_model_video2NPX, Behavior.Start);
Behavior.End            = predict(synch_model_video2NPX, Behavior.End);

partner_names           = animal_types;
partner_names(ismember(animal_types, repeated_animal))                  = [];

config.Behavior         = Behavior;
config.repeated_animal  = repeated_animal;
config.animal_types     = animal_types        ;
config.play_behaviors   = play_behaviors      ;
config.beh_bin          = bin_size             ;
config.conv_length      = conv_length;
config.behavior_window  = 0;


[play_bouts_table]      = play_bout(config);




%%  load lfp from current dir
disp('LOADING LFP')
if exist([current_dir,'\','LFP_PAG.mat'], 'file')==2

    NPX_Type        = 2;
    load([current_dir,'\','LFP_PAG.mat'], 'LFP')
elseif exist([current_dir,'\','LFP_PAG.dat'], 'file')==2
    NPX_Type        = 1;
    file_pointer    = fopen([current_dir,'\','LFP_PAG.dat'], 'r');
    LFP             = fread(file_pointer,'int16');
    LFP             = reshape(LFP, 384, numel(LFP)/384);
end


disp('LFP LOADED')
%% -------------------- SELECT PAG CHANNEL(S) --------------------
disp('Loading Channel Map')
hard_coded_x_coords = [8 40;258 290; 508 540; 758 790];
area_limit = readtable(area_limit_table);

load([current_dir,'\chann_map_',area2analyze,'.mat'], 'xcoords', 'ycoords','chanMap')
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
%% obtain_psth
hist_range      = [-20 20];
range_time_wrap = [-5 5];

n_bins_time_wrap = range(range_time_wrap)*100;

PAG_LFP         = double(LFP(mid_PAG_channel,:));
clear LFP
LFP_time        = (1:size(PAG_LFP,2))/2500;

spect_bin_size = wind_length-wind_overlap;

range2exctract  = LFP_time>=min(play_bouts_table(:))+hist_range(1) & LFP_time<=max(play_bouts_table(:))+hist_range(2);
psth_struct     = [];
for ch_n=1:numel(mid_PAG_channel)
    disp('Estimating spectrogram')
    [pow_spectrogram,~,spect_time]  = spectrogram(PAG_LFP(ch_n,range2exctract),wind_length*2500, wind_overlap*2500, f,2500);
    spect_time                      = spect_time+min(play_bouts_table(:))+hist_range(1);

    pow_spectrogram = abs(pow_spectrogram);
    disp('ready')

    play_bout_onset          = nan(size(play_bouts_table,1), size(pow_spectrogram,1),round(range(hist_range)/spect_bin_size));
    play_bout_offset         = nan(size(play_bouts_table,1), size(pow_spectrogram,1), round(range(hist_range)/spect_bin_size));
    play_bout_tw_this        = nan(size(play_bouts_table,1), size(pow_spectrogram,1), round(range(range_time_wrap)/spect_bin_size)  ...
         + n_bins_time_wrap);

    for pb_n=1:size(play_bouts_table,1)

        play_bout_start     = play_bouts_table(pb_n,1);
        play_bout_end       = play_bouts_table(pb_n,2);
        [~,loc_start]       = min(abs(spect_time-play_bout_start));
        [~,loc_end]         = min(abs(spect_time-play_bout_end));

        entire_range        = round(loc_start+(hist_range(1)/spect_bin_size):loc_start+(hist_range(2)/spect_bin_size));
        allowed_index_fres = ismember(entire_range,1:size(pow_spectrogram,2));
        play_bout_onset(pb_n,:,allowed_index_fres) = pow_spectrogram(:,entire_range(allowed_index_fres));

        entire_range        = round(loc_end+(hist_range(1)/spect_bin_size):loc_end+(hist_range(2)/spect_bin_size));
        allowed_index_fres = ismember(entire_range,1:size(pow_spectrogram,2));
        play_bout_offset(pb_n,:,allowed_index_fres) = pow_spectrogram(:,entire_range(allowed_index_fres));

        pre_time            = loc_start+ round(range_time_wrap(1)/spect_bin_size):loc_start-1;
        post_time           = loc_end+1:loc_end+round(range_time_wrap(2)/spect_bin_size);
        in_between_time     = loc_start:loc_end;
        time_wrapped_freq_pow  = interp1(in_between_time, pow_spectrogram(:,in_between_time)', linspace(loc_start,loc_end,n_bins_time_wrap));
        play_bout_tw_this(pb_n,:,:) =  [pow_spectrogram(:,pre_time),time_wrapped_freq_pow', pow_spectrogram(:,post_time)];

    end
    
    if ch_n==1
        psth_struct.play_bout_onset         = play_bout_onset;
        psth_struct.play_bout_offset        = play_bout_offset;
        psth_struct.play_bout_tw_this       = play_bout_tw_this;
        psth_struct.hist_range              = hist_range;
        psth_struct.range_time_wrap         = range_time_wrap;
        psth_struct.n_bins_time_wrap        = n_bins_time_wrap;
        psth_struct.wind_length             = wind_length;
        psth_struct.wind_overlap            = wind_overlap;
        psth_struct.n_bins_time_wrap        = n_bins_time_wrap;
        psth_struct.ch                      = mid_PAG_channel(ch_n);
        psth_struct.play_bouts_table        = play_bouts_table;
        psth_struct.Behavior                = Behavior;
        psth_struct.f                       = f;
        psth_struct.freq_pow_range          = freq_pow_range;
        psth_struct.min_separation          = min_separation;
        psth_struct.pow_spectrogram         = pow_spectrogram;
        psth_struct.spect_time              = spect_time;


    else
        psth_struct(ch_n).play_bout_onset           = play_bout_onset;
        psth_struct(ch_n).play_bout_offset          = play_bout_offset;
        psth_struct(ch_n).play_bout_tw_this         = play_bout_tw_this;
        psth_struct(ch_n).hist_range                = hist_range;
        psth_struct(ch_n).range_time_wrap           = range_time_wrap;
        psth_struct(ch_n).n_bins_time_wrap          = n_bins_time_wrap;
        psth_struct(ch_n).wind_length               = wind_length;
        psth_struct(ch_n).wind_overlap              = wind_overlap;
        psth_struct(ch_n).n_bins_time_wrap          = n_bins_time_wrap;
        psth_struct(ch_n).ch                        = mid_PAG_channel(ch_n);
        psth_struct(ch_n).play_bouts_table          = play_bouts_table;
        psth_struct(ch_n).Behavior                  = Behavior;
        psth_struct(ch_n).f                         = f;
        psth_struct(ch_n).freq_pow_range             = freq_pow_range;
        psth_struct(ch_n).min_separation            = min_separation;
        psth_struct(ch_n).pow_spectrogram           = pow_spectrogram;
    end

end
end
