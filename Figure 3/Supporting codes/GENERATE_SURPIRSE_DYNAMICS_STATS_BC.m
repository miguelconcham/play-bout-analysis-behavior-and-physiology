function synch_structure = GENERATE_SURPIRSE_DYNAMICS_STATS_BC(npx_data_dir, bin_size, hist_range, time_precision, areas2compare)


% synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
area2analyse        = 'PAG';
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\Area_limits_GoodLooking.xlsx';
behavior_data       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Behavior backups';
play_behaviors      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};
% non_play_behaviors  = {'Grooming', 'PounceI','Rearing', 'Sniffing','Scratching', 'Bite'};
% play_behaviors      = non_play_behaviors;
% npx_raw_data = 
animal_code         = strsplit(npx_data_dir, '\');
animal_code         = animal_code{end};
animal_code_params  = strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
repeated_animal     = animal_code_params{3};
nperm               = 2000;
smoth_wind_sec      = bin_size;
jitter_multiplier  = 1;


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
%% select only channels within areas2compare

chanel_selection = find(ismember(areas_by_channel,areas2compare));

%% load spikes
spike_times             = double(readNPY([npx_data_dir,'\spike_times_', area2analyze, '.npy']))/30000;
spike_clusters          = readNPY([npx_data_dir,'\spike_clusters_', area2analyze, '.npy']);
cluster_info            = readtable([npx_data_dir,'\cluster_info_', area2analyze, '.tsv'] ,"FileType","text",'Delimiter', '\t');
cluster_info            = cluster_info(ismember(cluster_info.group,{'mua', 'good'}),:);
these_neurons_areas     = areas_by_channel(cluster_info.ch+1);
cluster_info.area       = these_neurons_areas;
cluster_info.ch         = cluster_info.ch+1;
disp(['Selecting ' , num2str(sum(ismember(cluster_info.ch,chanel_selection))), ' out of ' , num2str(size(cluster_info,1)), ' clusters'])
cluster_info = cluster_info(ismember(cluster_info.ch,chanel_selection),:);



%% obtain_psth





% --- Step 1: pre alocated arrays ---



synch_comb = nchoosek(1:size(cluster_info,1), 2);
time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));
baseline_index = time_centers>0;
synch_spikes_histogram          = nan(2,size(synch_comb,1),numel(psth_edges)-1);
synch_spikes_pctl               = nan(2,size(synch_comb,1),numel(psth_edges)-1);


for n_comb = 1:size(synch_spikes_histogram,2)

    clust1  = synch_comb(n_comb,1);
    clust2  = synch_comb(n_comb,2);

    spikes_cluseter_1 = spike_times(spike_clusters==cluster_info.cluster_id(clust1));
    spikes_cluseter_2 = spike_times(spike_clusters==cluster_info.cluster_id(clust2));

    synch_spikes_histogram_n1 = nan(size(play_bouts_table,1),numel(psth_edges)-1);
    synch_spikes_histogram_n2 = nan(size(play_bouts_table,1),numel(psth_edges)-1);

    synch_spikes_histogram_shifted_n1 = nan(size(play_bouts_table,1),numel(psth_edges)-1,nperm);
    synch_spikes_histogram_shifted_n2 = nan(size(play_bouts_table,1),numel(psth_edges)-1,nperm);
    pb_lengths = diff(play_bouts_table');
    [~, length_order] = sort(pb_lengths);

    play_bouts_table = play_bouts_table(length_order,:);

    for pb_n=1:size(play_bouts_table,1)

        pb_start            = play_bouts_table(pb_n,1);
        this_pb_spikes_1    = spikes_cluseter_1(spikes_cluseter_1>=pb_start+hist_range(1) & spikes_cluseter_1<=pb_start+hist_range(2));
        this_pb_spikes_2    = spikes_cluseter_2(spikes_cluseter_2>=pb_start+hist_range(1) & spikes_cluseter_2<=pb_start+hist_range(2));

        if ~isempty(this_pb_spikes_2) && ~isempty(this_pb_spikes_1)
            s1_index = any(abs(this_pb_spikes_1-this_pb_spikes_2')<time_precision,2);
            s2_index = any(abs(this_pb_spikes_1-this_pb_spikes_2')<time_precision,1);

            synch_spikes_histogram_n1(pb_n,:) = histcounts(this_pb_spikes_1(s1_index)-pb_start,psth_edges);
            synch_spikes_histogram_n2(pb_n,:) = histcounts(this_pb_spikes_2(s2_index)-pb_start,psth_edges);
           
        end

    end


    for perm_n=1:nperm
    for pb_n=1:size(play_bouts_table,1)

        pb_start            = play_bouts_table(pb_n,1);
        this_pb_spikes_1    = spikes_cluseter_1(spikes_cluseter_1>=pb_start+hist_range(1) & spikes_cluseter_1<=pb_start+hist_range(2));
        this_pb_spikes_2    = spikes_cluseter_2(spikes_cluseter_2>=pb_start+hist_range(1) & spikes_cluseter_2<=pb_start+hist_range(2));
        direction  = (rand(size(this_pb_spikes_1))>.5)*2 -1;
        this_pb_spikes_shifted_1    = this_pb_spikes_1 + (rand(size(this_pb_spikes_1))*jitter_multiplier*time_precision + time_precision).*direction;
        direction  = (rand(size(this_pb_spikes_2))>.5)*2 -1;
        this_pb_spikes_shifted_2    = this_pb_spikes_2  + (rand(size(this_pb_spikes_2))*jitter_multiplier*time_precision + time_precision).*direction;
        this_pb_spikes_shifted_1    = mod(this_pb_spikes_shifted_1-(pb_start+hist_range(1)),range(hist_range))+(pb_start+hist_range(1));
        this_pb_spikes_shifted_2    = mod(this_pb_spikes_shifted_2-(pb_start+hist_range(1)),range(hist_range))+(pb_start+hist_range(1));
        
        if ~isempty(this_pb_spikes_2) && ~isempty(this_pb_spikes_1)           

            s1_index = any(abs(this_pb_spikes_1-this_pb_spikes_shifted_2')<time_precision,2);
            s2_index = any(abs(this_pb_spikes_shifted_1-this_pb_spikes_2')<time_precision,1);

            synch_spikes_histogram_shifted_n1(pb_n,:,perm_n) = histcounts(this_pb_spikes_1(s1_index)-pb_start,psth_edges);
            synch_spikes_histogram_shifted_n2(pb_n,:,perm_n) = histcounts(this_pb_spikes_2(s2_index)-pb_start,psth_edges);
        end

    end
    end

    synch_spikes_histogram(1,n_comb,:) = mean(synch_spikes_histogram_n1, 'omitmissing');
    synch_spikes_histogram(2,n_comb,:) = mean(synch_spikes_histogram_n2, 'omitmissing');        

    real_mean = mean(synch_spikes_histogram_n1, 'omitmissing');
    real_mean = movmean(real_mean, smoth_wind_sec/bin_size);
    real_mean = (real_mean - mean(real_mean(baseline_index)))/std(real_mean(baseline_index));
    surroageted_mean = squeeze(mean(synch_spikes_histogram_shifted_n1,1, 'omitmissing'));
    for perm_n=1:nperm
        surroageted_mean(:,perm_n) = movmean( surroageted_mean(:,perm_n) , smoth_wind_sec/bin_size);
        surroageted_mean(:,perm_n) = (surroageted_mean(:,perm_n) - mean(surroageted_mean(baseline_index,perm_n)))/std(surroageted_mean(baseline_index,perm_n));

    end
    pctl_Val = mean(surroageted_mean'>=real_mean);
    synch_spikes_pctl(1,n_comb,:) = pctl_Val;

    real_mean = mean(synch_spikes_histogram_n2, 'omitmissing');
    real_mean = movmean(real_mean, smoth_wind_sec/bin_size);
    real_mean = (real_mean - mean(real_mean(baseline_index)))/std(real_mean(baseline_index));
    surroageted_mean = squeeze(mean(synch_spikes_histogram_shifted_n2,1, 'omitmissing'));
    for perm_n=1:nperm
        surroageted_mean(:,perm_n) = movmean( surroageted_mean(:,perm_n) , smoth_wind_sec/bin_size);
        surroageted_mean(:,perm_n) = (surroageted_mean(:,perm_n) - mean(surroageted_mean(baseline_index,perm_n)))/std(surroageted_mean(baseline_index,perm_n));
    end
    pctl_Val = mean(surroageted_mean'>=real_mean);
    synch_spikes_pctl(2,n_comb,:) = pctl_Val;

    if mod(n_comb, 500) == 0
        disp([num2str(n_comb), ' out of ', num2str(size(synch_spikes_histogram,2))])
    end

end

synch_structure.synch_spikes_histogram          = synch_spikes_histogram;
synch_structure.synch_spikes_pctl               = synch_spikes_pctl;
synch_structure.psth_edges                      = psth_edges;
synch_structure.these_neurons_areas             = these_neurons_areas;
synch_structure.cluster_info                    = cluster_info;
synch_structure.synch_comb                      = synch_comb;


end