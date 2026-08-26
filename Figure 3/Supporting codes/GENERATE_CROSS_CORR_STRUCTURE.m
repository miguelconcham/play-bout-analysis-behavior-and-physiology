function cross_corr_structure = GENERATE_CROSS_CORR_STRUCTURE(npx_data_dir, bin_size, hist_range )


% synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Synch data';
synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Synch data';
area2analyse        = 'PAG';
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Area_limits_GoodLooking.xlsx';
behavior_data       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Behavior backups';
play_behaviors      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};
non_play_behaviors  = {'Grooming', 'PounceI','Rearing', 'Sniffing','Scratching', 'Bite'};

% npx_raw_data = 
animal_code         = strsplit(npx_data_dir, '\');
animal_code         = animal_code{end};
animal_code_params  = strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
repeated_animal     = animal_code_params{3};

area2analyze = 'PAG';
%% define parameters
n_rand              = 1000;

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



config.Behavior         = Behavior;
config.repeated_animal  = repeated_animal;
config.animal_types     = animal_types        ;
config.play_behaviors   = non_play_behaviors      ;
config.beh_bin          = bin_size             ;
config.conv_length      = conv_length;
config.behavior_window  = 0;

[non_play_bouts_table]      = play_bout(config);


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


%% obtain_psth





% --- Step 1: pre alocated arrays ---



cross_corr_comb = nchoosek(1:size(cluster_info,1), 2);


all_cross_corr_entire_session   = nan(size(cross_corr_comb,1),numel(psth_edges)-1);
all_cross_corr_play             = zeros(size(cross_corr_comb,1),numel(psth_edges)-1);
all_cross_corr_non_play         = zeros(size(cross_corr_comb,1),numel(psth_edges)-1);
all_cross_corr_non_play_pctl    = zeros(size(cross_corr_comb,1),numel(psth_edges)-1);
all_cross_corr_play_pctl    = zeros(size(cross_corr_comb,1),numel(psth_edges)-1);

disp('ESTIMATING CROSS CORRELATIONS')
for n_comb = 1:size(all_cross_corr_entire_session,1)

    clust1  = cross_corr_comb(n_comb,1);
    clust2  = cross_corr_comb(n_comb,2);

    spikes_cluseter_1 = spike_times(spike_clusters==cluster_info.cluster_id(clust1));
    spikes_cluseter_2 = spike_times(spike_clusters==cluster_info.cluster_id(clust2));


    tsOffsets = crosscorrelogram(spikes_cluseter_1, spikes_cluseter_2, hist_range);

    all_cross_corr_entire_session(n_comb,:) = histcounts(tsOffsets,psth_edges);



    for pb_n=1:size(play_bouts_table,1)

        pb_start = play_bouts_table(pb_n,1);
        % pb_end   = play_bouts_table(pb_n,2);
        pb_end = pb_start;
        this_pb_spikes_1    = spikes_cluseter_1(spikes_cluseter_1>=pb_start+hist_range(1) & spikes_cluseter_1<=pb_end+hist_range(2));
        this_pb_spikes_2    = spikes_cluseter_2(spikes_cluseter_2>=pb_start+hist_range(1) & spikes_cluseter_2<=pb_end+hist_range(2));
        if ~isempty(this_pb_spikes_1) &&  ~isempty(this_pb_spikes_2)
            tsOffsets           = crosscorrelogram(this_pb_spikes_2, this_pb_spikes_1, hist_range);
            all_cross_corr_play(n_comb,:) = all_cross_corr_play(n_comb,:)+ histcounts(tsOffsets,psth_edges);
        end
    end
    % shrink to  %% ---------- Permutation cross-correlation across bouts ----------
    N = size(play_bouts_table,1);
    all_pairs = nchoosek(1:N, 2);  % size M x 2, M = N*(N-1)/2

    M = nchoosek(N,2);      % total number of unique pairs
    % max_permutations = min(10, nchoosek(M, N));
    max_permutations = 10;

    perm_crosscorr_play = zeros(max_permutations, length(psth_edges)-1);

    used_sets = containers.Map(); % track unique sets

    perm_count = 0;
    attempts = 0;
    max_attempts = max_permutations * 50;

    while perm_count < max_permutations && attempts < max_attempts
        attempts = attempts + 1;

        % randomly pick N rows from all_pairs
        idx = randperm(M, N);
        chosen_pairs = all_pairs(idx, :);

        % sort each row and flatten to string to check uniqueness
        sorted_pairs = sort(chosen_pairs, 2);
        key = sprintf('%d,', sorted_pairs'); % flatten all numbers

        if isKey(used_sets, key)
            continue
        end

        used_sets(key) = true;

        % compute cross-correlation for this permutation
        this_perm_hist = zeros(1, length(psth_edges)-1);
        for k = 1:N
            i = chosen_pairs(k,1);
            j = chosen_pairs(k,2);

            % Bout i (neuron 1)
            pb_start_i = play_bouts_table(i,1);
            % pb_end_i   = play_bouts_table(i,2); 
            pb_end_i = pb_start_i;
            spikes1 = spikes_cluseter_1(spikes_cluseter_1 >= pb_start_i + hist_range(1) & ...
                spikes_cluseter_1 <= pb_end_i   + hist_range(2));

            % Bout j (neuron 2)
            pb_start_j = play_bouts_table(j,1);
            % pb_end_j   = play_bouts_table(j,2);
            pb_end_j = pb_start_j;
            spikes2 = spikes_cluseter_2(spikes_cluseter_2 >= pb_start_j + hist_range(1) & ...
                spikes_cluseter_2 <= pb_end_j   + hist_range(2));

            if isempty(spikes1) || isempty(spikes2)
                continue
            end

            % align spikes
            spikes1_al = spikes1 - pb_start_i;
            spikes2_al = spikes2 - pb_start_j;

            % cross-correlogram
            tsOffsets_perm = crosscorrelogram(spikes2_al, spikes1_al, hist_range);
            this_perm_hist = this_perm_hist + histcounts(tsOffsets_perm, psth_edges);
        end

        perm_count = perm_count + 1;
        perm_crosscorr_play(perm_count,:) = this_perm_hist;
    end

    % shrink to actual size
    perm_crosscorr_play = perm_crosscorr_play(1:perm_count,:);

    all_cross_corr_play_pctl(n_comb,:) = mean(perm_crosscorr_play)/N;
    all_cross_corr_play(n_comb,:) = all_cross_corr_play(n_comb,:)/N;
    if all(this_perm_hist == 0)
        continue  % skip permutations with no valid spikes at all
    end


    for pb_n=1:size(non_play_bouts_table,1)
        pb_start = non_play_bouts_table(pb_n,1);
        % pb_end   = non_play_bouts_table(pb_n,2);
        pb_end = pb_start;
        this_pb_spikes_1    = spikes_cluseter_1(spikes_cluseter_1>=pb_start +hist_range(1) & spikes_cluseter_1<=pb_end +hist_range(2));
        this_pb_spikes_2    = spikes_cluseter_2(spikes_cluseter_2>=pb_start +hist_range(1) & spikes_cluseter_2<=pb_end +hist_range(2));
        if ~isempty(this_pb_spikes_1) &&  ~isempty(this_pb_spikes_2)
            tsOffsets           = crosscorrelogram(this_pb_spikes_2, this_pb_spikes_1, hist_range);
            all_cross_corr_non_play(n_comb,:) =all_cross_corr_non_play(n_comb,:)+ histcounts(tsOffsets,psth_edges);
        end
    end
    
    
    % shrink to  %% ---------- Permutation cross-correlation across bouts ----------
    N = size(non_play_bouts_table,1);
    all_pairs = nchoosek(1:N, 2);  % size M x 2, M = N*(N-1)/2

    M = nchoosek(N,2);      % total number of unique pairs
    % max_permutations = min(10, nchoosek(M, N));
    max_permutations = 10;
    perm_crosscorr_non_play = zeros(max_permutations, length(psth_edges)-1);

    used_sets = containers.Map(); % track unique sets

    perm_count = 0;
    attempts = 0;
    max_attempts = max_permutations * 50;

    while perm_count < max_permutations && attempts < max_attempts
        attempts = attempts + 1;

        % randomly pick N rows from all_pairs
        idx = randperm(M, N);
        chosen_pairs = all_pairs(idx, :);

        % sort each row and flatten to string to check uniqueness
        sorted_pairs = sort(chosen_pairs, 2);
        key = sprintf('%d,', sorted_pairs'); % flatten all numbers

        if isKey(used_sets, key)
            continue
        end

        used_sets(key) = true;

        % compute cross-correlation for this permutation
        this_perm_hist = zeros(1, length(psth_edges)-1);
        for k = 1:N
            i = chosen_pairs(k,1);
            j = chosen_pairs(k,2);

            % Bout i (neuron 1)
            pb_start_i = non_play_bouts_table(i,1);
            % pb_end_i   = non_play_bouts_table(i,2);
            pb_end_i = pb_start_i;
            spikes1 = spikes_cluseter_1(spikes_cluseter_1 >= pb_start_i + hist_range(1) & ...
                spikes_cluseter_1 <= pb_end_i   + hist_range(2));

            % Bout j (neuron 2)
            pb_start_j = non_play_bouts_table(j,1);
            % pb_end_j   = non_play_bouts_table(j,2);
            pb_end_j = pb_start_j;
            spikes2 = spikes_cluseter_2(spikes_cluseter_2 >= pb_start_j + hist_range(1) & ...
                spikes_cluseter_2 <= pb_end_j   + hist_range(2));

            if isempty(spikes1) || isempty(spikes2)
                continue
            end

            % align spikes
            spikes1_al = spikes1 - pb_start_i;
            spikes2_al = spikes2 - pb_start_j;

            % cross-correlogram
            tsOffsets_perm = crosscorrelogram(spikes2_al, spikes1_al, hist_range);
            this_perm_hist = this_perm_hist + histcounts(tsOffsets_perm, psth_edges);
        end

        perm_count = perm_count + 1;
        perm_crosscorr_non_play(perm_count,:) = this_perm_hist;
    end

    % shrink to actual size
    perm_crosscorr_non_play = perm_crosscorr_non_play(1:perm_count,:);
  

    all_cross_corr_non_play_pctl(n_comb,:) = mean(perm_crosscorr_non_play)/N;
    all_cross_corr_non_play(n_comb,:)  = all_cross_corr_non_play(n_comb,:) /N;
    if all(this_perm_hist == 0)
        continue  % skip permutations with no valid spikes at all
    end

end

cross_corr_structure.all_cross_corr_entire_session  = all_cross_corr_entire_session;
cross_corr_structure.all_cross_corr_play            = all_cross_corr_play;
cross_corr_structure.all_cross_corr_non_play        = all_cross_corr_non_play;
cross_corr_structure.all_cross_corr_non_play_pctl   = all_cross_corr_non_play_pctl;
cross_corr_structure.all_cross_corr_play_pctl       = all_cross_corr_play_pctl;
cross_corr_structure.these_neurons_areas            = these_neurons_areas;
cross_corr_structure.cluster_info                   = cluster_info;
cross_corr_structure.cross_corr_comb                = cross_corr_comb;
end

    

