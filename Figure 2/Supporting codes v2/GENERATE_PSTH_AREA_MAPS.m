function psth_struct = GENERATE_PSTH_AREA_MAPS(current_dir, wind_length, wind_overlap, min_separation, f, freq_pow_range)

%% Parameters
play_behaviors    = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD'};
synch_directory   = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
chan_map_folder   = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\StarndarChannMap';
area_limit_table  = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\Area_limits_GoodLooking.xlsx';
behavior_data     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Behavior backups';

animal_code        = strsplit(current_dir, '\');
animal_code        = animal_code{end};
animal_code_params = strsplit(animal_code, ' ');
animal_batch       = animal_code_params{1};
repeated_animal    = animal_code_params{3};

%% Load synchronization model
load([synch_directory, '\', animal_code, '\synch_model_video2NPX.mat'])

%% Load and preprocess behavior
Behavior_file = [behavior_data, '\', animal_code, '.txt'];
Behavior = readtable(Behavior_file);
Behavior(:, 2) = [];
Behavior.Properties.VariableNames = {'Animal', 'Start', 'End', 'Length', 'Type'};

bin_size       = 0.01;
conv_length    = 1;
Behavior.Type2 = Behavior.Type;
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_A', 'Pounce_B'}))   = {'Pounce'};
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
Behavior.Type2(ismember(Behavior.Type2, ''))                           = {'Other'};
Behavior(ismember(Behavior.Animal, 'Reversal'), :) = [];

animal_types = unique(Behavior.Animal);
animal_types(ismember(animal_types, 'Session_structure')) = [];

Behavior.Start = predict(synch_model_video2NPX, Behavior.Start);
Behavior.End   = predict(synch_model_video2NPX, Behavior.End);

config.Behavior        = Behavior;
config.repeated_animal = repeated_animal;
config.animal_types    = animal_types;
config.play_behaviors  = play_behaviors;
config.beh_bin         = bin_size;
config.conv_length     = conv_length;
config.behavior_window = 0;

[play_bouts_table] = play_bout(config);

%% Load LFP
disp('LOADING LFP')
if exist([current_dir, '\', 'LFP_PAG.mat'], 'file') == 2
    NPX_Type = 2;
    load([current_dir, '\', 'LFP_PAG.mat'], 'LFP')
elseif exist([current_dir, '\', 'LFP_PAG.dat'], 'file') == 2
    NPX_Type = 1;
    file_pointer = fopen([current_dir, '\', 'LFP_PAG.dat'], 'r');
    LFP = fread(file_pointer, 'int16');
    LFP = reshape(LFP, 384, numel(LFP) / 384);
end
disp('LFP LOADED')

%% Build channel map and area labels
disp('Loading Channel Map')
areas_by_channel    = cell(384, 1);
channel_map         = nan(384, 2);
hard_coded_x_coords = [8 40; 258 290; 508 540; 758 790];
area_limit = readtable(area_limit_table);

if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName, this_animal), :);

if NPX_Type == 1
    load([chan_map_folder, '\neuropixPhase3A_kilosortChanMap.mat'], 'xcoords', 'ycoords', 'chanMap')
    for ch_n = 1:384
        ch = chanMap(ch_n) + 1;
        channel_map(ch, 1) = xcoords(ch_n);
        channel_map(ch, 2) = ycoords(ch_n);
        areas_by_channel{ch} = area_limit.area{ycoords(ch_n) >= area_limit.depth_start & ycoords(ch_n) < area_limit.depth_end + 1 & ismember(area_limit.Probe_Area, 'PAG')};
    end
else
    load([current_dir, '\ChannelMap.mat'], 'xcoords', 'ycoords', 'chanMap')
    for ch_n = 1:384
        probe_n = find(any(ismember(hard_coded_x_coords, xcoords(ch_n)), 2));
        ch = chanMap(ch_n) + 1;
        channel_map(ch, 1) = xcoords(ch_n);
        channel_map(ch, 2) = ycoords(ch_n);
        areas_by_channel{ch} = area_limit.area{ycoords(ch_n) >= area_limit.depth_start & ycoords(ch_n) < area_limit.depth_end + 1 & area_limit.ProbeNum == probe_n};
    end
end

%% Compute band-limited power PSTH for all channels
hist_range       = [-20 20];
range_time_wrap  = [-5 5];
n_bins_time_wrap = range(range_time_wrap) * 100;
spect_bin_size   = wind_length - wind_overlap;

PAG_LFP  = double(LFP);
clear LFP
LFP_time = (1:size(PAG_LFP, 2)) / 2500;

range2exctract = LFP_time >= min(play_bouts_table(:)) + hist_range(1) & LFP_time <= max(play_bouts_table(:)) + hist_range(2);
n_time_bins   = round(range(hist_range) / spect_bin_size);
play_bout_onset = nan(384, size(play_bouts_table, 1), n_time_bins);

for ch_n = 1:size(PAG_LFP, 1)
    if mod(ch_n, 25) == 1
        disp(['Processing channel #', num2str(ch_n)])
    end

    [pow_spectrogram, ~, spect_time] = spectrogram(PAG_LFP(ch_n, range2exctract), wind_length * 2500, wind_overlap * 2500, f, 2500);
    spect_time = spect_time + min(play_bouts_table(:)) + hist_range(1);
    pow_spectrogram = abs(pow_spectrogram);

    freq_pow = mean(log10(pow_spectrogram(f >= freq_pow_range(1) & f <= freq_pow_range(2), :)));
    freq_pow = movmean(freq_pow, 1 / max(freq_pow_range));

    for pb_n = 1:size(play_bouts_table, 1)
        play_bout_start = play_bouts_table(pb_n, 1);
        [~, loc_start]  = min(abs(spect_time - play_bout_start));

        entire_range = round(loc_start + (hist_range(1) / spect_bin_size):loc_start + (hist_range(2) / spect_bin_size));
        allowed_index_fres = ismember(entire_range, 1:size(freq_pow, 2));
        play_bout_onset(ch_n, pb_n, allowed_index_fres) = freq_pow(entire_range(allowed_index_fres));
    end
end

psth_struct.play_bout_onset   = play_bout_onset;
psth_struct.play_bouts_table  = play_bouts_table;
psth_struct.channel_map       = channel_map;
psth_struct.areas_by_channel  = areas_by_channel;
psth_struct.hist_range        = hist_range;
psth_struct.wind_length       = wind_length;
psth_struct.wind_overlap      = wind_overlap;
psth_struct.range_time_wrap   = range_time_wrap;
psth_struct.n_bins_time_wrap  = n_bins_time_wrap;
end
