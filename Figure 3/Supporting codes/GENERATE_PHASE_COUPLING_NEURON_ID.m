function phase_struct = GENERATE_PHASE_COUPLING_NEURON_ID(current_dir,Hd_freq1,bin_size_freq, id_list, plot_bool )


data_root           = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
synch_directory     = [data_root, '\Synch data'];
play_behaviors      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};
chan_map_folder     = [data_root, '\NPX data\StarndarChannMap'];
area_limit_table    = [data_root, '\Area_limits_GoodLooking.xlsx'];
behavior_data       = [data_root, '\Behavior backups'];
 area2analyse = 'PAG';
animal_code         = strsplit(current_dir, '\');
animal_code         = animal_code{end};
animal_code_params  = strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
repeated_animal     = animal_code_params{3};

area2analyze = 'PAG';
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

%% determining NPX type
hard_coded_x_coords_NPX2 = [8 40;258 290; 508 540; 758 790];
load([current_dir,'\','chann_map_', area2analyse, '.mat'], 'chanMap', 'xcoords', 'ycoords')
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

LFP = double(LFP);
disp('LFP LOADED')

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
      areas_by_channel{ch} = area_limit.area{ycoords(ch_n)>=area_limit.depth_start &  ycoords(ch_n)<area_limit.depth_end+1 & area_limit.ProbeNum==probe_n & ismember(area_limit.Probe_Area, area2analyse)};
  end

end   
%% load spikes
spike_times             = double(readNPY([current_dir,'\spike_times_', area2analyze, '.npy']))/30000;
spike_clusters          = readNPY([current_dir,'\spike_clusters_', area2analyze, '.npy']);
cluster_info            = readtable([current_dir,'\cluster_info_', area2analyze, '.tsv'] ,"FileType","text",'Delimiter', '\t');
cluster_info            = cluster_info(ismember(cluster_info.group,{'mua', 'good'}),:);
these_neurons_areas     = areas_by_channel(cluster_info.ch+1);
cluster_info.area       = these_neurons_areas;
cluster_info.ch         = cluster_info.ch+1;

%% find intended neurons


neuron_index = ismember(cluster_info.cluster_id,id_list);
cluster_info = cluster_info(neuron_index,:);
channel_list = unique(cluster_info.ch);






PAG_LFP         = double(LFP);
clear LFP
sr_LFP = 2500;


% --- Step 1: real data binning ---
lfp_time                        = (1:size(PAG_LFP,2))/sr_LFP;
session_phase_stats             = nan(size(partner_sessions,1), size(cluster_info,1),5);
entire_recording_phase_stats    = nan(1, size(cluster_info,1),5);
session_psth                    = nan(size(partner_sessions,1), size(cluster_info,1),numel(edges_freq)-1);
entire_recording_psth           = nan(1, size(cluster_info,1),numel(edges_freq)-1);
assigned                        = false( size(cluster_info,1),1);
all_psth                        = cell( size(cluster_info,1),size(partner_sessions,1)+1);
all_loc_times                   = cell( size(cluster_info,1),size(partner_sessions,1)+1);
all_phases                      = cell( size(cluster_info,1),size(partner_sessions,1)+1);
all_lfp                         = cell(size(cluster_info,1),4);
all_spikes                      = cell(size(cluster_info,1),1);
spikes4raster                   = cell(size(cluster_info,1),1);

for ch_n=channel_list

    disp(['Processing ch #', num2str(ch_n)])


    filtered_signal_freq    = filtfilt(Hd_freq1.Coefficients, 1, PAG_LFP(ch_n,:));
    hiblert_data            = hilbert(filtered_signal_freq);
    phase_data              = angle(hiblert_data);
    amplitud_data           = abs(hiblert_data);

    std_amp                 = std(amplitud_data);
    [~,max_locs_freq]       = findpeaks(filtered_signal_freq, 'MinPeakProminence',.5*std_amp, 'MinPeakDistance', sr_LFP/(Hd_freq1.CutoffFrequency2  )) ;

    original_distribution   = phase_data(max_locs_freq);
    mean_angle_original     = angle(mean(exp(1i*original_distribution)));
    mean_angle_original     = mod( mean_angle_original+2*pi,2*pi);
    phase_data              = mod(phase_data  - mean_angle_original + 5*pi , 2*pi) - pi; %% centering step

    phase_data              = circular_uniformize(phase_data);


    neurons_this_channel  = find(ismember(cluster_info.ch, ch_n))';
    disp([num2str(numel(neurons_this_channel)), ' neuron(s) in this channel'])
    disp([num2str(sum(assigned==0)), ' neuron(s) to be analyzed'])

    for nn=neurons_this_channel
        all_lfp{nn,1} = PAG_LFP(ch_n,:);
        all_lfp{nn,2} = filtered_signal_freq;
        all_lfp{nn,3} = lfp_time;
        all_lfp{nn,4} = phase_data;

        this_id     = cluster_info.cluster_id(nn);
        this_spikes = spike_times(spike_clusters==this_id);
        all_spikes{nn} = this_spikes;

        this_neuron_phases           = interp1(lfp_time, phase_data,this_spikes);
        all_phases{nn,1}             = this_neuron_phases;

        %estimating entire recording  MVL, p_val, prefered angle, mean rate
        %and the corresdponign psth
        this_neuron_r               = circ_r(this_neuron_phases);
        this_neuron_prefered_angle  = circ_mean(this_neuron_phases);

        rand_r            = nan(n_rand,1);
        for nr=1:n_rand
            rand_r(nr) =circ_r((randsample(phase_data,numel(this_neuron_phases), false))');
        end

        this_neuron_p       = sum(rand_r>this_neuron_r)/n_rand;
        [~,~, mean_rate] = kernel_rate(this_spikes, lfp_time(1), lfp_time(end), 0.01, 0.1);


        if isempty(this_neuron_prefered_angle)
            this_neuron_prefered_angle = NaN;
        end
        if isempty(this_neuron_p)
            this_neuron_p = NaN;
        end
        if isempty(this_neuron_r)
            this_neuron_r = NaN;
        end
        if isempty(mean_rate)
            mean_rate = NaN;
        end

        entire_recording_phase_stats(1,nn,:) =  [this_neuron_prefered_angle this_neuron_r this_neuron_p mean_rate this_id];
        entire_recording_psth_this_neuron  = nan(numel(max_locs_freq),numel(edges_freq)-1);
    
        if plot_bool

            figure
            hold on
        end
        spikes4raster{nn,1} = [];
        for j=1:numel(max_locs_freq)
            this_peak_time                  = lfp_time(max_locs_freq(j));
            spikes_this_peak                = this_spikes(this_spikes>=this_peak_time+psth_range_freq(1) & this_spikes<=this_peak_time+psth_range_freq(2))-this_peak_time;
            spikes_this_peak = reshape(spikes_this_peak, numel(spikes_this_peak),1);
            spikes4raster{nn,1} = [spikes4raster{nn,1};[spikes_this_peak spikes_this_peak*0+j]];
            entire_recording_psth_this_neuron(j,:)    = histcounts(spikes_this_peak,edges_freq);
            if plot_bool
                plot(spikes_this_peak, spikes_this_peak*0 + j, 'k.')
            end
        end
        if plot_bool
           xlim(psth_range_freq)
           title(nn)
           pause(.1)
        end

        all_psth{nn,1}          = entire_recording_psth_this_neuron;
        all_loc_times{nn,1}     = lfp_time(max_locs_freq);

        entire_recording_psth(1,nn,:) = mean(entire_recording_psth_this_neuron);
        assigned(nn) = true;

        if all(isnan(mean(entire_recording_psth_this_neuron)))
            disp('entire psth is nan??')
        end



        for sn=1:size(partner_sessions,1)

            this_session_time       = find(lfp_time>=partner_sessions.Start(sn) & lfp_time<=partner_sessions.End(sn));
            this_session_peaks      = max_locs_freq(ismember(max_locs_freq,this_session_time));


            %estimating each partner session recording  MVL, p_val,
            % prefered angle, mean rate and the corresdponign psth
            psth_this_session_freq  = nan(numel(this_session_peaks),numel(edges_freq)-1);

             spikes4raster{nn,sn+1} = [];
             for j=1:numel(this_session_peaks)
                 this_peak_time                  = lfp_time(this_session_peaks(j));
                 spikes_this_peak                = this_spikes(this_spikes>=this_peak_time+psth_range_freq(1) & this_spikes<=this_peak_time+psth_range_freq(2))-this_peak_time;
                 spikes_this_peak = reshape(spikes_this_peak, numel(spikes_this_peak),1);
                 spikes4raster{nn,sn+1} = [spikes4raster{nn,sn+1};[spikes_this_peak spikes_this_peak*0+j]];

                 psth_this_session_freq(j,:)    = histcounts(spikes_this_peak,edges_freq);
             end
            all_psth{nn,1+sn}          = psth_this_session_freq;
            all_loc_times{nn,1+sn}     = lfp_time(this_session_peaks);
            session_psth(sn,nn,:)       = mean(psth_this_session_freq);


            this_session_spikes         = this_spikes(this_spikes>=partner_sessions.Start(sn) & this_spikes<=partner_sessions.End(sn));

            [~,~, mean_rate] = kernel_rate(this_session_spikes, partner_sessions.Start(sn), partner_sessions.End(sn), 0.01, 0.1);
            this_neuron_phases_freq     = interp1(lfp_time, phase_data,this_session_spikes);
            this_session_phases_freq    = phase_data(lfp_time>=partner_sessions.Start(sn) & lfp_time<=partner_sessions.End(sn));

            all_phases{nn,1+sn}         = this_neuron_phases_freq;
            this_neuron_r               = circ_r(this_neuron_phases_freq);
            this_neuron_prefered_angle  = circ_mean(this_neuron_phases_freq);

            rand_r = nan(n_rand,1);

            for nr=1:n_rand

                rand_r(nr) =circ_r((randsample(this_session_phases_freq,numel(this_neuron_phases_freq), false))');
            end

            this_neuron_p = sum(rand_r>this_neuron_r)/n_rand;

            if isempty(this_neuron_prefered_angle)
                this_neuron_prefered_angle = NaN;
            end
            if isempty(this_neuron_p)
                this_neuron_p = NaN;
            end
            if isempty(this_neuron_r)
                this_neuron_r = NaN;
            end
            if isempty(mean_rate)
                mean_rate = NaN;
            end

            session_phase_stats(sn,nn,:) = [this_neuron_prefered_angle this_neuron_r this_neuron_p mean_rate this_id];
        end
    end


end
    
if any(~assigned)

    disp('Not all neurons assigned')
end
phase_struct.cluster_info                       = cluster_info;
phase_struct.session_psth                       = session_psth;
phase_struct.session_phase_stats                = session_phase_stats;
phase_struct.entire_recording_phase_stats       = entire_recording_phase_stats;
phase_struct.entire_recording_psth              = entire_recording_psth;
phase_struct.session_phase_stats                = session_phase_stats;
phase_struct.psth_range_freq                    = psth_range_freq;
phase_struct.partner_sessions                   = partner_sessions;
phase_struct.edges_freq                         = edges_freq;
phase_struct.areas_by_channel                   = areas_by_channel;
phase_struct.channel_map                        = channel_map;
phase_struct.all_psth                           = all_psth;
phase_struct.all_loc_times                      = all_loc_times;
phase_struct.all_phases                         = all_phases;
phase_struct.all_lfp                            = all_lfp;
phase_struct.all_spikes                         = all_spikes;
phase_struct.spikes4raster                      = spikes4raster;
end
