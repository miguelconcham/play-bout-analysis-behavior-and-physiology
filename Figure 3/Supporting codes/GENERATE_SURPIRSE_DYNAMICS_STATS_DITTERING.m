function synch_structure = GENERATE_SURPIRSE_DYNAMICS_STATS_DITTERING(npx_data_dir, bin_size, hist_range, time_precision, areas2compare,behaviors4playbout)


% synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
area2analyse        = 'PAG';
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\Area_limits_GoodLooking.xlsx';
behavior_data       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Behavior backups';

animal_code         = strsplit(npx_data_dir, '\');
animal_code         = animal_code{end};
animal_code_params  = strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
repeated_animal     = animal_code_params{3};
nperm               = 2000;
smoth_wind_sec      = 0.05;
dittering_range    = 0.009;


% Make permutation/jitter reproducible per animal/session.
% Without this, regenerated *.mat files will differ run-to-run.
rng(sum(double(animal_code)), 'twister');


area2analyze = 'PAG';
%% define parameters


psth_edges = hist_range(1):bin_size:hist_range(2);


%% load synch from synch folder
load([synch_directory,'\', animal_code, '\synch_model_video2NPX.mat'], 'synch_model_video2NPX')

%% 2 Load beahvior data from behavior folder


Behavior_file =[behavior_data,'\', animal_code,'.txt'];%load behavior data

Behavior                            = readtable(Behavior_file);
Behavior(:,2)                       = [];
Behavior.Properties.VariableNames   = {'Animal', 'Start', 'End', 'Length', 'Type'};



Behavior.Type2          = Behavior.Type;

Behavior.Type2(ismember(Behavior.Type2, {'Pounce_A', 'Pounce_B'}))      = {'Pounce'}; %% Merging behaviors to Type2
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'}))    = {'PounceI'};
Behavior.Type2(ismember( Behavior.Type2,''))                            = {'Other'};
Behavior(ismember(Behavior.Animal, 'Reversal'),:)                       = [];

Behavior.Start          = predict(synch_model_video2NPX, Behavior.Start);
Behavior.End            = predict(synch_model_video2NPX, Behavior.End);


partner_sessions = Behavior(strcmp(Behavior.Type2, 'Partners session'),:);



%% 2 Create play bout array

conv_length             = 1;
play_bout_bins_size     = 0.01;

animal_types            = unique(Behavior.Animal);

animal_types(ismember(animal_types,'Session_structure'))                =[];




config.Behavior         = Behavior;
config.repeated_animal  = repeated_animal;
config.animal_types     = animal_types        ;
config.play_behaviors   = behaviors4playbout      ;
config.beh_bin          = play_bout_bins_size             ;
config.conv_length      = conv_length;
config.behavior_window  = 0;


[play_bouts_table]      = play_bout(config);


%
% config.Behavior         = Behavior;
% config.repeated_animal  = repeated_animal;
% config.animal_types     = animal_types        ;
% config.play_behaviors   = non_play_behaviors      ;
% config.beh_bin          = bin_size             ;
% config.conv_length      = conv_length;
% config.behavior_window  = 0;
%
% [non_play_bouts_table]      = play_bout(config);


%% determining NPX type
hard_coded_x_coords_NPX2 = [8 40;258 290; 508 540; 758 790];
load([npx_data_dir,'\','chann_map_', area2analyse, '.mat'], 'chanMap', 'xcoords', 'ycoords')
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
        areas_by_channel{ch} = area_limit.area{ycoords(ch_n)>=area_limit.depth_start &  ycoords(ch_n)<area_limit.depth_end+1 & ismember(area_limit.Probe_Area, area2analyse) };
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



%% obtain_psth





% --- Step 1: pre alocated arrays ---



synch_comb = nchoosek(1:size(cluster_info,1), 2);


synch_spikes_histogram          = nan(2,size(synch_comb,1),numel(psth_edges)-1);
synch_spikes_pctl               = nan(2,size(synch_comb,1),numel(psth_edges)-1);


for n_comb = 1:size(synch_spikes_histogram,2)


    if mod(n_comb, 25) == 1
        disp([num2str(n_comb), ' out of ', num2str(size(synch_spikes_histogram,2))])
    end



    clust1  = synch_comb(n_comb,1);
    clust2  = synch_comb(n_comb,2);

    spikes_cluseter_1 = spike_times(spike_clusters==cluster_info.cluster_id(clust1));
    spikes_cluseter_2 = spike_times(spike_clusters==cluster_info.cluster_id(clust2));

    synch_spikes_histogram_n1 = nan(size(play_bouts_table,1),numel(psth_edges)-1);
    synch_spikes_histogram_n2 = nan(size(play_bouts_table,1),numel(psth_edges)-1);

    synch_spikes_histogram_shifted_n1 = nan(size(play_bouts_table,1),numel(psth_edges)-1,nperm);
    synch_spikes_histogram_shifted_n2 = nan(size(play_bouts_table,1),numel(psth_edges)-1,nperm);


    for pb_n=1:size(play_bouts_table,1)

        pb_start            = play_bouts_table(pb_n,1);
        this_pb_spikes_1    = spikes_cluseter_1(spikes_cluseter_1>=pb_start+hist_range(1) & spikes_cluseter_1<=pb_start+hist_range(2));
        this_pb_spikes_2    = spikes_cluseter_2(spikes_cluseter_2>=pb_start+hist_range(1) & spikes_cluseter_2<=pb_start+hist_range(2));

        if ~isempty(this_pb_spikes_2) && ~isempty(this_pb_spikes_1)
                [s1_index, s2_index] = find_synchronous_spikes_twopointer( ...
                this_pb_spikes_1, this_pb_spikes_2, time_precision);
           
            synch_spikes_histogram_n1(pb_n,:) = histcounts(this_pb_spikes_1(s1_index)-pb_start,psth_edges);
            synch_spikes_histogram_n2(pb_n,:) = histcounts(this_pb_spikes_2(s2_index)-pb_start,psth_edges);

        end

    end

    for perm_n=1:nperm
        for pb_n=1:size(play_bouts_table,1)

            pb_start            = play_bouts_table(pb_n,1);
            this_pb_spikes_1    = spikes_cluseter_1(spikes_cluseter_1>=pb_start+hist_range(1) & spikes_cluseter_1<=pb_start+hist_range(2));
            this_pb_spikes_2    = spikes_cluseter_2(spikes_cluseter_2>=pb_start+hist_range(1) & spikes_cluseter_2<=pb_start+hist_range(2));
            if rand<0.5
                this_pb_spikes_shifted_1    = this_pb_spikes_1 + 2*(rand(size(this_pb_spikes_1)) -.5)*dittering_range ;
                this_pb_spikes_shifted_1    = mod(this_pb_spikes_shifted_1-(pb_start+hist_range(1)),range(hist_range))+(pb_start+hist_range(1));
                this_pb_spikes_shifted_2    = this_pb_spikes_2;
            else
                this_pb_spikes_shifted_2    = this_pb_spikes_2 +2*(rand(size(this_pb_spikes_2))-.5)*dittering_range;
                this_pb_spikes_shifted_2    = mod(this_pb_spikes_shifted_2-(pb_start+hist_range(1)),range(hist_range))+(pb_start+hist_range(1));
                this_pb_spikes_shifted_1    = this_pb_spikes_1;
            end

            if ~isempty(this_pb_spikes_2) && ~isempty(this_pb_spikes_1)


                [s1_index_shift, s2_index_shift] = find_synchronous_spikes_twopointer( ...
                    this_pb_spikes_shifted_1, this_pb_spikes_shifted_2, time_precision);

                synch_spikes_histogram_shifted_n1(pb_n,:,perm_n) = histcounts(this_pb_spikes_shifted_1(s1_index_shift)-pb_start,psth_edges);
                synch_spikes_histogram_shifted_n2(pb_n,:,perm_n) = histcounts(this_pb_spikes_shifted_2(s2_index_shift)-pb_start,psth_edges);
            end

        end
    end

    synch_spikes_histogram(1,n_comb,:) = mean(synch_spikes_histogram_n1, 'omitmissing');
    synch_spikes_histogram(2,n_comb,:) = mean(synch_spikes_histogram_n2, 'omitmissing');

    real_mean = mean(synch_spikes_histogram_n1, 'omitmissing');
    real_mean = movmean(real_mean, smoth_wind_sec/bin_size);
    surroageted_mean = squeeze(mean(synch_spikes_histogram_shifted_n1,1, 'omitmissing'));
    for perm_n=1:nperm
        surroageted_mean(:,perm_n) = movmean( surroageted_mean(:,perm_n) , smoth_wind_sec/bin_size);
    end
    pctl_Val = (1 + sum(surroageted_mean' >= real_mean, 1)) / (nperm + 1);
    synch_spikes_pctl(1,n_comb,:) = pctl_Val;

    real_mean = mean(synch_spikes_histogram_n2, 'omitmissing');
    real_mean = movmean(real_mean, smoth_wind_sec/bin_size);
    surroageted_mean = squeeze(mean(synch_spikes_histogram_shifted_n2,1, 'omitmissing'));
    for perm_n=1:nperm
        surroageted_mean(:,perm_n) = movmean( surroageted_mean(:,perm_n) , smoth_wind_sec/bin_size);
    end
    pctl_Val = (1 + sum(surroageted_mean' >= real_mean, 1)) / (nperm + 1);
    synch_spikes_pctl(2,n_comb,:) = pctl_Val;


end

synch_structure.synch_spikes_histogram          = synch_spikes_histogram;
synch_structure.synch_spikes_pctl               = synch_spikes_pctl;
synch_structure.psth_edges                      = psth_edges;
synch_structure.these_neurons_areas             = these_neurons_areas;
synch_structure.cluster_info                    = cluster_info;
synch_structure.synch_comb                      = synch_comb;
synch_structure.nperm                           = nperm;
synch_structure.smoth_wind_sec                  = smoth_wind_sec;
synch_structure.dittering_range                = dittering_range;
synch_structure.play_behaviors                 = behaviors4playbout;
synch_structure.dittering_range                = dittering_range;



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