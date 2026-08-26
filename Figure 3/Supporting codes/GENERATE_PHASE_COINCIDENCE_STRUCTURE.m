function coindicence_phase_struct = GENERATE_PHASE_COINCIDENCE_STRUCTURE(npx_data_dir,Hd_freq1,bin_size_freq , hist_range,time_precision,areas2compare)


% synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Synch data';
synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Synch data';
area2analyze        = 'PAG';
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Area_limits_GoodLooking.xlsx';
behavior_data       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Behavior backups';
play_behaviors      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};
% non_play_behaviors  = {'Grooming', 'PounceI','Rearing', 'Sniffing','Scratching', 'Bite'};

% npx_raw_data =
animal_code         = strsplit(npx_data_dir, '\');
animal_code         = animal_code{end};
animal_code_params  = strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
repeated_animal     = animal_code_params{3};


%% define parameters
n_rand              = 1000;
psth_range_freq    = round(1.25*[-1 1]./min(Hd_freq1.CutoffFrequency1),2);

edges_freq = psth_range_freq(1):bin_size_freq:psth_range_freq(2);


%% load synch from synch folder
load([synch_directory,'\', animal_code, '\synch_model_video2NPX.mat'], 'synch_model_video2NPX')

%% 2 Load beahvior data from behavior folder


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

Behavior.Start          = predict(synch_model_video2NPX, Behavior.Start);
Behavior.End            = predict(synch_model_video2NPX, Behavior.End);


partner_sessions = Behavior(strcmp(Behavior.Type2, 'Partners session'),:);



%% 2 Create play bout array

bin_size                = 0.01;
conv_length             = 1;
animal_types            = unique(Behavior.Animal);

animal_types(ismember(animal_types,'Session_structure'))                =[];




config.Behavior         = Behavior;
config.repeated_animal  = repeated_animal;
config.animal_types     = animal_types        ;
config.play_behaviors   = play_behaviors      ;
config.beh_bin          = bin_size             ;
config.conv_length      = conv_length;
config.behavior_window  = 0;


[play_bouts_table]      = play_bout(config);




%% determining NPX type
hard_coded_x_coords_NPX2 = [8 40;258 290; 508 540; 758 790];
load([npx_data_dir,'\','chann_map_', area2analyze, '.mat'], 'chanMap', 'xcoords', 'ycoords')
if any(ismember(xcoords, hard_coded_x_coords_NPX2))
    NPX_Type        = 2;
else
    NPX_Type        = 1;
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
end

%% Create (load) channel map
disp('Loading Channel Map')
areas_by_channel = cell(384,1);
channel_map      = nan(384,2);

area_limit = readtable(area_limit_table);

% Build animal identifier for area selection
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal,animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName,this_animal),:);

if NPX_Type == 1


    for ch_n=1:384
        ch = chanMap(ch_n);
        channel_map(ch,1) = xcoords(ch_n);
        channel_map(ch,2) = ycoords(ch_n);
        areas_by_channel{ch} = area_limit.area{ycoords(ch_n)>=area_limit.depth_start &  ycoords(ch_n)<area_limit.depth_end+1 & ismember(area_limit.Probe_Area, area2analyze) };
    end
else



    for ch_n=1:384
        probe_n = find(any(ismember(hard_coded_x_coords_NPX2,xcoords(ch_n)),2));
        ch = chanMap(ch_n);
        channel_map(ch,1) = xcoords(ch_n);
        channel_map(ch,2) = ycoords(ch_n);
        areas_by_channel{ch} = area_limit.area{ycoords(ch_n)>=area_limit.depth_start &  ycoords(ch_n)<area_limit.depth_end+1 & area_limit.ProbeNum==probe_n & ismember(area_limit.Probe_Area, 'PAG')};
    end

end


%% load spikes
spike_times             = double(readNPY([npx_data_dir,'\spike_times_', area2analyze, '.npy']))/30000;
spike_clusters          = readNPY([npx_data_dir,'\spike_clusters_', area2analyze, '.npy']);
cluster_info            = readtable([npx_data_dir,'\cluster_info_', area2analyze, '.tsv'] ,"FileType","text",'Delimiter', '\t');
cluster_info            = cluster_info(ismember(cluster_info.group,{'mua', 'good'}),:);
these_neurons_areas     = areas_by_channel(cluster_info.ch+1);
cluster_info.area       = these_neurons_areas;
cluster_info.ch         = cluster_info.ch+1;


disp(['Selecting ' , num2str(sum(ismember(cluster_info.area,areas2compare))), ' out of ' , num2str(size(cluster_info,1)), ' clusters'])
cluster_info =       cluster_info(ismember(cluster_info.area,areas2compare),:);


channels_with_spikes    = unique(cluster_info.ch)';

%%  load lfp from current dir
disp('LOADING LFP')
if NPX_Type==2
    load([npx_data_dir,'\','LFP_', area2analyze, '.mat'], 'LFP')
elseif  NPX_Type        == 1

    file_pointer    = fopen([npx_data_dir,'\','LFP_',area2analyze,'.dat'], 'r');
    LFP             = fread(file_pointer,'int16');
    LFP             = reshape(LFP, 384, numel(LFP)/384);
end

disp('LFP LOADED')

%% prealocating varaibles, estimating phase for channels with spikes


PAG_LFP                         = double(LFP(channels_with_spikes,:));
filtered_signal_all_channels    = PAG_LFP;
amplitud_data_all_channels      = PAG_LFP;
phase_data_all_channels         = PAG_LFP;
clear LFP
sr_LFP = 2500;
figure('units','normalized','outerposition',[0 0 1 1]);
sp_n =1;
for nch=1:size(PAG_LFP,1)

    if sp_n>12
        sp_n = 1;
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    subplot(3,4,sp_n)



    filtered_signal_all_channels(nch,:) = filtfilt(Hd_freq1.Coefficients, 1, filtered_signal_all_channels(nch,:));
    hiblert_data_this_channel           = hilbert(filtered_signal_all_channels(nch,:));
    phase_data_this_channel             = angle(hiblert_data_this_channel);
    amplitud_data_all_channels(nch,:)   = abs(hiblert_data_this_channel);

    std_amp                             = std(amplitud_data_all_channels(nch,:));
    [~,max_locs_freq]                   = findpeaks(filtered_signal_all_channels(nch,:), 'MinPeakProminence',.5*std_amp, 'MinPeakDistance', sr_LFP/(Hd_freq1.CutoffFrequency2  )) ;
    [ ~,min_locs]                       = findpeaks(-filtered_signal_all_channels(nch,:), 'MinPeakProminence',.5*std_amp, 'MinPeakDistance', sr_LFP/(Hd_freq1.CutoffFrequency2  )) ;


    phase_data_this_channel              = circular_uniformize(phase_data_this_channel);


    original_distribution   = phase_data_this_channel(max_locs_freq);
    mean_angle_original     = angle(mean(exp(1i*original_distribution)));
    mean_angle_original     = mod( mean_angle_original+2*pi,2*pi);
    phase_data_this_channel = mod(phase_data_this_channel  - mean_angle_original + 5*pi , 2*pi) - pi; %% centering step

    phase_data_all_channels(nch,:) = phase_data_this_channel;
    polarhistogram(phase_data_this_channel, -pi:(pi/36):pi, 'EdgeColor','none', 'Normalization','percentage')
    hold on
    polarhistogram(phase_data_this_channel(max_locs_freq), -pi:(pi/36):pi, 'EdgeColor','none', 'Normalization','percentage')
    polarhistogram(phase_data_this_channel(min_locs), -pi:(pi/36):pi, 'EdgeColor','none', 'Normalization','percentage')
    title(channels_with_spikes(nch))
    pause(.1)
    sp_n = sp_n+1;


end
%%

synch_comb                      = nchoosek(1:size(cluster_info,1), 2);

%%

% --- Step 1: pre alocated arrays ---
lfp_time                        = (1:size(PAG_LFP,2))/sr_LFP;

%% Estimating phases of coincidences

%% Preallocate output
% sync_event_data{n_comb} = [spike1_rel, spike2_rel, lag, phase11, phase12, phase21, phase22, period]
sync_event_data = cell(size(synch_comb,1),2);

%% Optional speedup: cache spike trains for all clusters
cluster_spike_times = cell(height(cluster_info),1);
for c = 1:height(cluster_info)
    cluster_spike_times{c} = spike_times(spike_clusters == cluster_info.cluster_id(c));
end

%% Main loop over neuron pairs
for n_comb = 1:size(synch_comb,1)

    if mod(n_comb, 25) == 1
        disp([num2str(n_comb), ' out of ', num2str(size(synch_comb,1))])
    end

    clust1 = synch_comb(n_comb,1);
    clust2 = synch_comb(n_comb,2);

    spikes_cluster_1 = cluster_spike_times{clust1};
    spikes_cluster_2 = cluster_spike_times{clust2};

    % Get correct phase channel for each cluster
    ch1 = cluster_info.ch(clust1);
    ch2 = cluster_info.ch(clust2);

    idx_ch1 = (channels_with_spikes == ch1);
    idx_ch2 = (channels_with_spikes == ch2);

    phase_data_ch_1 = phase_data_all_channels(idx_ch1,:);
    phase_data_ch_2 = phase_data_all_channels(idx_ch2,:);

    % Store all synchronous events for this pair here
    pair_events_2 = [];
    pair_events_1 = [];

    %% Loop over play bouts
    for pb_n = 1:size(play_bouts_table,1)

        pb_start = play_bouts_table(pb_n,1);

        this_pb_spikes_1 = spikes_cluster_1( ...
            spikes_cluster_1 >= pb_start + hist_range(1) & ...
            spikes_cluster_1 <= pb_start + hist_range(2));

        this_pb_spikes_2 = spikes_cluster_2( ...
            spikes_cluster_2 >= pb_start + hist_range(1) & ...
            spikes_cluster_2 <= pb_start + hist_range(2));

        if ~isempty(this_pb_spikes_1) && ~isempty(this_pb_spikes_2)

            [s1_index, s2_index] = find_synchronous_spikes_twopointer( ...
                this_pb_spikes_1, this_pb_spikes_2, time_precision);

            if any(s1_index) && any(s2_index)

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % NEURON 1 synchronous spikes
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sync_spikes_1 = this_pb_spikes_1(s1_index);
                spike1_rel    = sync_spikes_1 - pb_start;

                phase11 = interp1(lfp_time, phase_data_ch_1, sync_spikes_1, 'linear', NaN);
                phase21 = interp1(lfp_time, phase_data_ch_2, sync_spikes_1, 'linear', NaN);

                valid_1 = ~isnan(phase11) & ~isnan(phase21);

                spike1_rel = spike1_rel(valid_1);
                phase11    = phase11(valid_1);
                phase21    = phase21(valid_1);

                if ~isempty(spike1_rel)
                    period1 = ones(size(spike1_rel));
                    period1(spike1_rel < 0) = 2;

                    these_events_1 = [ ...
                        spike1_rel(:), ...
                        phase11(:), ...
                        phase21(:), ...
                        period1(:)];

                    pair_events_1 = [pair_events_1; these_events_1];
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % NEURON 2 synchronous spikes
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sync_spikes_2 = this_pb_spikes_2(s2_index);
                spike2_rel    = sync_spikes_2 - pb_start;

                phase12 = interp1(lfp_time, phase_data_ch_1, sync_spikes_2, 'linear', NaN);
                phase22 = interp1(lfp_time, phase_data_ch_2, sync_spikes_2, 'linear', NaN);

                valid_2 = ~isnan(phase12) & ~isnan(phase22);

                spike2_rel = spike2_rel(valid_2);
                phase12    = phase12(valid_2);
                phase22    = phase22(valid_2);

                if ~isempty(spike2_rel)
                    period2 = ones(size(spike2_rel));
                    period2(spike2_rel < 0) = 2;

                    these_events_2 = [ ...
                        spike2_rel(:), ...
                        phase12(:), ...
                        phase22(:), ...
                        period2(:)];

                    pair_events_2 = [pair_events_2; these_events_2];
                end

            end
        end
    end

    % Save this pair's event matrix
    sync_event_data{n_comb,1} = pair_events_1;
    sync_event_data{n_comb,2} = pair_events_2;

end

%% storing section
coindicence_phase_struct.sync_event_data        = sync_event_data;
coindicence_phase_struct.sync_event_columns     = {'spike1_rel','spike2_rel','lag','phase11','phase12','phase21','phase22','period'};
coindicence_phase_struct.cluster_spike_times    = cluster_spike_times;
coindicence_phase_struct.synch_comb             = synch_comb;
coindicence_phase_struct.cluster_info           = cluster_info;
coindicence_phase_struct.sync_event_data        = sync_event_data;           % Synchronous spike arrays per neuron
coindicence_phase_struct.sync_event_columns     = {'spike_rel','phase_on_ch1','phase_on_ch2','period'}; % column names
coindicence_phase_struct.cluster_spike_times    = cluster_spike_times;      % raw spike times for all clusters
coindicence_phase_struct.cluster_info           = cluster_info;             % cluster metadata (area, channel, etc.)
coindicence_phase_struct.channels_with_spikes   = channels_with_spikes;     % channels included in analysis
coindicence_phase_struct.hist_range             = hist_range;               % relative window around play bouts
coindicence_phase_struct.time_precision         = time_precision;           % coincidence detection precision
coindicence_phase_struct.play_bouts_table       = play_bouts_table;         % play bout start times



    function [s1_index, s2_index] = find_synchronous_spikes_twopointer(s1, s2, tol)

        n1 = numel(s1);
        n2 = numel(s2);

        s1_index = false(n1,1);
        s2_index = false(n2,1);

        j_start = 1;

        for i = 1:n1

            % Move lower bound forward until s2 is within left tolerance
            while j_start <= n2 && s2(j_start) < s1(i) - tol
                j_start = j_start + 1;
            end

            j = j_start;

            % Mark all s2 spikes within tolerance of s1(i)
            while j <= n2 && s2(j) <= s1(i) + tol
                s1_index(i) = true;
                s2_index(j) = true;
                j = j + 1;
            end
        end
    end
end
