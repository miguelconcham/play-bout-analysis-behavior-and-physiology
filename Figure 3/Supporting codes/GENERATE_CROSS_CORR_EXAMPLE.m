function [] = GENERATE_CROSS_CORR_EXAMPLE(npx_data_dir, bin_size, hist_range, id_list ,time_precision)


% synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
area2analyse        = 'PAG';
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\Area_limits_GoodLooking.xlsx';
behavior_data       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Behavior backups';
play_behaviors      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};
non_play_behaviors  = {'Grooming', 'PounceI','Rearing', 'Sniffing','Scratching', 'Bite'};

% npx_raw_data = 
animal_code         = strsplit(npx_data_dir, '\');
animal_code         = animal_code{end};
animal_code_params  = strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
repeated_animal     = animal_code_params{3};
n_sec_to_shfift = [1 2];
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



spike_times     = spike_times(ismember(spike_clusters,id_list));
spike_clusters  = spike_clusters(ismember(spike_clusters,id_list));

cluster_info = cluster_info(ismember(cluster_info.cluster_id,id_list),:)


these_neurons_areas     = areas_by_channel(cluster_info.ch+1);
cluster_info.area       = these_neurons_areas;
cluster_info.ch         = cluster_info.ch+1;


%% obtain_psth





% --- Step 1: pre alocated arrays ---



cross_corr_comb = nchoosek(1:size(cluster_info,1), 2);


all_cross_corr_entire_session   = nan(size(cross_corr_comb,1),numel(psth_edges)-1);
all_cross_corr_play             = zeros(size(cross_corr_comb,1),numel(psth_edges)-1);


for n_comb = 1:size(all_cross_corr_entire_session,1)

    clust1  = cross_corr_comb(n_comb,1);
    clust2  = cross_corr_comb(n_comb,2);

    spikes_cluseter_1 = spike_times(spike_clusters==cluster_info.cluster_id(clust1));
    spikes_cluseter_2 = spike_times(spike_clusters==cluster_info.cluster_id(clust2));


    % tsOffsets                               = crosscorrelogram(spikes_cluseter_1, spikes_cluseter_2, hist_range);
    % all_cross_corr_entire_session(n_comb,:) = histcounts(tsOffsets,psth_edges);
    % 
    % 
    % randomzied_cross_corr = nan(n_rand,numel(all_cross_corr_entire_session));
    % 
    % for nr= 1:n_rand
    %     tsOffsets = crosscorrelogram(spikes_cluseter_1, mod(spikes_cluseter_2+rand*range(n_sec_to_shfift) + n_sec_to_shfift(1),max(spikes_cluseter_2)), hist_range);
    %     randomzied_cross_corr(nr,:) = histcounts(tsOffsets,psth_edges);
    % end
    % 
    % pctl = mean(movmean(randomzied_cross_corr,10,2)<movmean(all_cross_corr_entire_session,10));
    % figure
    % plot(time_centers,pctl)
%%
    x_lim = [-1 2];
    figure
    subplot(3,1,1)
    hold on
  
   synch_spikes_histogram_n1 = nan(size(play_bouts_table,1),numel(psth_edges)-1);
     synch_spikes_histogram_n2 = nan(size(play_bouts_table,1),numel(psth_edges)-1);
    num_coincidence_before_after = nan(size(play_bouts_table,1),2);
     pb_lengths = diff(play_bouts_table');
     [~, length_order] = sort(pb_lengths);

     play_bouts_table = play_bouts_table(length_order,:);
     ploting_pb_n = 1;
    for pb_n=1:size(play_bouts_table,1)

        pb_start = play_bouts_table(pb_n,1);
        pb_end   = play_bouts_table(pb_n,2);
        pb_end = min(pb_end, pb_start+hist_range(2));
      
        this_pb_spikes_1    = spikes_cluseter_1(spikes_cluseter_1>=pb_start+hist_range(1) & spikes_cluseter_1<=pb_end+hist_range(2));
        this_pb_spikes_2    = spikes_cluseter_2(spikes_cluseter_2>=pb_start+hist_range(1) & spikes_cluseter_2<=pb_end+hist_range(2));

        
        s1_index = any(abs(this_pb_spikes_1-this_pb_spikes_2')<time_precision,2);
         s2_index = any(abs(this_pb_spikes_1-this_pb_spikes_2')<time_precision,1);
        if numel(this_pb_spikes_1)>0 & numel(this_pb_spikes_2)>0


         fill([pb_start pb_end pb_end pb_start]-pb_start, [0 0 1 1]+ploting_pb_n, 'k', 'FaceAlpha',.1, 'EdgeColor','none')

         others_pb_in_ragne = (play_bouts_table(:,1) <pb_start+hist_range(2) & play_bouts_table(:,1)>pb_start+hist_range(1)) | ...
                                (play_bouts_table(:,2) <pb_start+hist_range(2) & play_bouts_table(:,2)>pb_start+hist_range(1));
         others_pb_in_ragne(pb_n) = 0;

         if any(others_pb_in_ragne)
             list_pb = find(others_pb_in_ragne);

             for o_pb_n = 1:numel(list_pb)

                 o_pb_start = play_bouts_table(o_pb_n,1);
                 o_pb_end   = play_bouts_table(o_pb_n,2);
                 o_pb_end = min(o_pb_end, pb_start+hist_range(2));

                    fill([o_pb_start o_pb_end o_pb_end o_pb_start]-pb_start, [0 0 1 1]+ploting_pb_n, 'k', 'FaceAlpha',.1, 'EdgeColor','none')
             end
         end



        plot([this_pb_spikes_1,this_pb_spikes_1]'-pb_start, [this_pb_spikes_1*0+ploting_pb_n,this_pb_spikes_1*0+ploting_pb_n+.5]' , 'k')
        plot([this_pb_spikes_2,this_pb_spikes_2]'-pb_start, [this_pb_spikes_2*0+ploting_pb_n,this_pb_spikes_2*0+ploting_pb_n+.5]' +.5, 'b')


         plot([this_pb_spikes_1(s1_index),this_pb_spikes_1(s1_index)]'-pb_start, [this_pb_spikes_1(s1_index)*0+ploting_pb_n,this_pb_spikes_1(s1_index)*0+ploting_pb_n+.5]' , 'r')
         plot([this_pb_spikes_2(s2_index),this_pb_spikes_2(s2_index)]'-pb_start, [this_pb_spikes_2(s2_index)*0+ploting_pb_n,this_pb_spikes_2(s2_index)*0+ploting_pb_n+.5]' +.5 , 'r')

         fill([pb_start pb_end pb_end pb_start]-pb_start, [0 1 1 0]+ploting_pb_n, 'k', 'FaceAlpha',.1, 'EdgeColor','none')
         ploting_pb_n = ploting_pb_n+1;
        end
      
         synch_spikes_histogram_n1(pb_n,:) = histcounts(this_pb_spikes_1(s1_index)-pb_start,psth_edges);
         synch_spikes_histogram_n2(pb_n,:) = histcounts(this_pb_spikes_2(s2_index)-pb_start,psth_edges);


       
        if ~isempty(this_pb_spikes_1) &&  ~isempty(this_pb_spikes_2)
            tsOffsets           = crosscorrelogram(this_pb_spikes_2, this_pb_spikes_1, hist_range);
            all_cross_corr_play(n_comb,:) = all_cross_corr_play(n_comb,:)+ histcounts(tsOffsets,psth_edges);
        end
    end
    axis tight
    xlim(x_lim)
     title('original synch')
     subplot(3,1,3)
     time_centers =.5*(psth_edges(1:end-1) + psth_edges(2:end));

     plot(time_centers, movmean(mean(synch_spikes_histogram_n1),10)/bin_size, 'k')
     hold on
     plot(time_centers, movmean(mean(synch_spikes_histogram_n2),10)/bin_size, 'b')
     xlim(x_lim)


 
      
    subplot(3,1,2)
    hold on
 
    desynch_spikes_histogram_n1 = nan(size(play_bouts_table,1),numel(psth_edges)-1);
     desynch_spikes_histogram_n2 = nan(size(play_bouts_table,1),numel(psth_edges)-1);
     pb_lengths = diff(play_bouts_table');
     [~, length_order] = sort(pb_lengths);

     play_bouts_table = play_bouts_table(length_order,:);
     anderes_play_bouts_table = play_bouts_table(randsample(size(play_bouts_table,1),size(play_bouts_table,1), false),:);
     ploting_pb_n = 1;
    for pb_n=1:size(play_bouts_table,1)

        pb_start = play_bouts_table(pb_n,1);
        pb_end   = play_bouts_table(pb_n,2);
        pb_end = min(pb_end, pb_start+hist_range(2));

        this_pb_spikes_1    = spikes_cluseter_1(spikes_cluseter_1>=pb_start+hist_range(1) & spikes_cluseter_1<=pb_start+hist_range(2));

        anderes_pb_start = anderes_play_bouts_table(pb_n,1);

        this_pb_spikes_2    = spikes_cluseter_2(spikes_cluseter_2>=anderes_pb_start+hist_range(1) & spikes_cluseter_2<=anderes_pb_start+hist_range(2));

        s1_index = any(abs((this_pb_spikes_1-pb_start)-(this_pb_spikes_2-anderes_pb_start)')<time_precision,2);
        s2_index = any(abs((this_pb_spikes_1-pb_start)-(this_pb_spikes_2-anderes_pb_start)')<time_precision,1);


        if numel(this_pb_spikes_1)>0 & numel(this_pb_spikes_2)>0

            fill([pb_start pb_end pb_end pb_start]-pb_start, [0 0 1 1]+ploting_pb_n, 'k', 'FaceAlpha',.1, 'EdgeColor','none')

            others_pb_in_ragne = (play_bouts_table(:,1) <pb_start+hist_range(2) & play_bouts_table(:,1)>pb_start+hist_range(1)) | ...
                (play_bouts_table(:,2) <pb_start+hist_range(2) & play_bouts_table(:,2)>pb_start+hist_range(1));
            others_pb_in_ragne(pb_n) = 0;

            if any(others_pb_in_ragne)
                list_pb = find(others_pb_in_ragne);

                for o_pb_n = 1:numel(list_pb)

                    o_pb_start = play_bouts_table(o_pb_n,1);
                    o_pb_end   = play_bouts_table(o_pb_n,2);
                    o_pb_end = min(o_pb_end, pb_start+hist_range(2));

                    fill([o_pb_start o_pb_end o_pb_end o_pb_start]-pb_start, [0 0 1 1]+ploting_pb_n, 'k', 'FaceAlpha',.1, 'EdgeColor','none')
                end
            end



            plot([this_pb_spikes_1,this_pb_spikes_1]'-pb_start, [this_pb_spikes_1*0+ploting_pb_n,this_pb_spikes_1*0+ploting_pb_n+.5]' , 'k')
            plot([this_pb_spikes_2,this_pb_spikes_2]'-anderes_pb_start, [this_pb_spikes_2*0+ploting_pb_n,this_pb_spikes_2*0+ploting_pb_n+.5]' +.5, 'b')


            plot([this_pb_spikes_1(s1_index),this_pb_spikes_1(s1_index)]'-pb_start, [this_pb_spikes_1(s1_index)*0+ploting_pb_n,this_pb_spikes_1(s1_index)*0+ploting_pb_n+.5]' , 'r')
            plot([this_pb_spikes_2(s2_index),this_pb_spikes_2(s2_index)]'-anderes_pb_start, [this_pb_spikes_2(s2_index)*0+ploting_pb_n,this_pb_spikes_2(s2_index)*0+ploting_pb_n+.5]' +.5, 'r')

            fill([pb_start pb_end pb_end pb_start]-pb_start, [0 1 1 0]+ploting_pb_n, 'k', 'FaceAlpha',.1, 'EdgeColor','none')
            ploting_pb_n = ploting_pb_n+1;
        end
        desynch_spikes_histogram_n1(pb_n,:) = histcounts(this_pb_spikes_1(s1_index)-pb_start,psth_edges);
        desynch_spikes_histogram_n2(pb_n,:) = histcounts(this_pb_spikes_2(s2_index)-anderes_pb_start,psth_edges);



        if ~isempty(this_pb_spikes_1) &&  ~isempty(this_pb_spikes_2)
            tsOffsets           = crosscorrelogram(this_pb_spikes_2, this_pb_spikes_1, hist_range);
            all_cross_corr_play(n_comb,:) = all_cross_corr_play(n_comb,:)+ histcounts(tsOffsets,psth_edges);
        end

    end
    axis tight
    xlim(x_lim)
  
   subplot(3,1,3)
   time_centers =.5*(psth_edges(1:end-1) + psth_edges(2:end));

   plot(time_centers, movmean(mean(desynch_spikes_histogram_n1),10)/bin_size, 'k:')
   hold on
   plot(time_centers, movmean(mean(desynch_spikes_histogram_n2),10)/bin_size, ':b')
   xlim(x_lim)
   title([animal_code, ' ID', num2str(id_list(1))])
   ylabel({num2str(id_list(2)), 'Synch events'})
   legend({'Synch event neuron 1','Synch event neuron 1','Synch event neuron 1 non matched playbout','Synch event neuron 2 non matched playbout'})


%%
end


end