function phase_struct = GENERATE_PHASE_COUPLING_STRUCTURE(npx_data_dir,Hd_freq1,bin_size_freq )


data_root           = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
synch_directory     = [data_root, '\Synch data'];
area2analyse        = 'PAG';
area_limit_table    = [data_root, '\Area_limits_GoodLooking.xlsx'];
behavior_data       = [data_root, '\Behavior backups'];
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
psth_alignment = 'peak'; %or peak
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

pre_play_bouts = build_balanced_periods(play_bouts_table);


config.Behavior         = Behavior;
config.repeated_animal  = repeated_animal;
config.animal_types     = animal_types        ;
config.play_behaviors   = non_play_behaviors      ;
config.beh_bin          = bin_size             ;
config.conv_length      = conv_length;
config.behavior_window  = 0;

[non_play_bouts_table]      = play_bout(config);
pre_non_play_bouts          = build_balanced_periods(non_play_bouts_table);


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
channels_with_spikes    = unique(cluster_info.ch)';

%%  load lfp from current dir
disp('LOADING LFP')
if NPX_Type==2
    load([npx_data_dir,'\','LFP_', area2analyse, '.mat'], 'LFP')
elseif  NPX_Type        == 1

    file_pointer    = fopen([npx_data_dir,'\','LFP_',area2analyse,'.dat'], 'r');
    LFP             = fread(file_pointer,'int16');
    LFP             = reshape(LFP, 384, numel(LFP)/384);
end

disp('LFP LOADED')

%% obtain_psth



PAG_LFP         = double(LFP);
clear LFP
sr_LFP = 2500;


% --- Step 1: pre alocated arrays ---
lfp_time                        = (1:size(PAG_LFP,2))/sr_LFP;
session_phase_stats             = nan(size(partner_sessions,1), size(cluster_info,1),7);
entire_recording_phase_stats    = nan(1, size(cluster_info,1),7);
play_phase_stats                = nan(1, size(cluster_info,1),7);
pre_play_phase_stats            = nan(1, size(cluster_info,1),7);
non_play_phase_stats            = nan(1, size(cluster_info,1),7);
pre_non_play_phase_stats        = nan(1, size(cluster_info,1),7);

session_psth                    = nan(size(partner_sessions,1), size(cluster_info,1),numel(edges_freq)-1);
entire_recording_psth           = nan(1, size(cluster_info,1),numel(edges_freq)-1);
play_psth                       = nan(size(partner_sessions,1), size(cluster_info,1),numel(edges_freq)-1);
pre_play_psth                   = nan(1, size(cluster_info,1),numel(edges_freq)-1);
non_play_psth                    = nan(size(partner_sessions,1), size(cluster_info,1),numel(edges_freq)-1);
pre_non_play_psth               = nan(1, size(cluster_info,1),numel(edges_freq)-1);
assigned                        = false( size(cluster_info,1),1);
for ch_n=channels_with_spikes
    
    disp(['Processing ch #', num2str(ch_n)])
   

    filtered_signal_freq    = filtfilt(Hd_freq1.Coefficients, 1, PAG_LFP(ch_n,:));
    hiblert_data            = hilbert(filtered_signal_freq);        
    phase_data              = angle(hiblert_data);
    amplitud_data           = abs(hiblert_data);

    std_amp                 = std(amplitud_data);

    if strcmp(psth_alignment,'trough')
        [~,peaksforalignment]       = findpeaks(-filtered_signal_freq, 'MinPeakProminence',.5*std_amp, 'MinPeakDistance', sr_LFP/(Hd_freq1.CutoffFrequency2  )) ;
    else
        [~,peaksforalignment]       = findpeaks(filtered_signal_freq, 'MinPeakProminence',.5*std_amp, 'MinPeakDistance', sr_LFP/(Hd_freq1.CutoffFrequency2  )) ;
    end



    
    phase_data              = circular_uniformize(phase_data);


    original_distribution   = phase_data(peaksforalignment);       
    mean_angle_original     = angle(mean(exp(1i*original_distribution)));
    mean_angle_original     = mod( mean_angle_original+2*pi,2*pi);  
    phase_data              = mod(phase_data  - mean_angle_original + 5*pi , 2*pi) - pi; %% centering step

    

    
    neurons_this_channel  = find(ismember(cluster_info.ch, ch_n))';
    disp([num2str(numel(neurons_this_channel)), ' neuron(s) in this channel'])
    disp([num2str(sum(assigned==0)), ' neuron(s) to be analyzed'])

    for nn=neurons_this_channel
        this_id     = cluster_info.cluster_id(nn);
        this_spikes = spike_times(spike_clusters==this_id);     

        
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %estimating entire recording  MVL, PPC,, p_val, prefered angle,
        % mean rate and the corresdponign psth 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        this_neuron_phases           = interp1(lfp_time, phase_data,this_spikes);        
        this_neuron_mvl                 = circ_r(this_neuron_phases);
        this_neuron_prefered_angle      = circ_mean(this_neuron_phases);
        this_neuron_ppc                 = compute_ppc_fast(this_neuron_phases);
        tic
        rand_r            = nan(n_rand,1);
        rand_ppc          = nan(n_rand,1);
        for nr=1:n_rand
            rand_r(nr) =circ_r((randsample(phase_data,numel(this_neuron_phases), false))');
            rand_ppc(nr) = compute_ppc_fast((randsample(phase_data,numel(this_neuron_phases), false))');
        end
        toc
        
        this_neuron_mvl_p       = sum(rand_r>this_neuron_mvl)/n_rand;
        this_neuron_ppc_p       = sum(rand_ppc>this_neuron_ppc)/n_rand;
        [dt,lambda, mean_rate]        = kernel_rate(this_spikes, lfp_time(1), lfp_time(end), 0.01, 0.1);


        if isempty(this_neuron_prefered_angle),this_neuron_prefered_angle = NaN; end
        if isempty(this_neuron_ppc), this_neuron_ppc = NaN; end
        if isempty(this_neuron_ppc_p), this_neuron_ppc_p = NaN; end
        if isempty(this_neuron_mvl_p),this_neuron_mvl_p = NaN; end
        if isempty(this_neuron_mvl),this_neuron_mvl = NaN;end
        if isempty(mean_rate),mean_rate = NaN;end

        entire_recording_phase_stats(1,nn,:) =  [this_neuron_prefered_angle this_neuron_mvl this_neuron_mvl_p this_neuron_ppc this_neuron_ppc_p mean_rate this_id];

        entire_recording_psth_this_neuron  = nan(numel(peaksforalignment),numel(edges_freq)-1);
        
        for j=1:numel(peaksforalignment)
            this_peak_time                  = lfp_time(peaksforalignment(j));
            spikes_this_peak                = this_spikes(this_spikes>=this_peak_time+psth_range_freq(1) & this_spikes<=this_peak_time+psth_range_freq(2))-this_peak_time;
            entire_recording_psth_this_neuron(j,:)    = histcounts(spikes_this_peak,edges_freq);
        end
        

        entire_recording_psth(1,nn,:) = mean(entire_recording_psth_this_neuron);
        assigned(nn) = true;

        if all(isnan(mean(entire_recording_psth_this_neuron)))
        disp('entire psth is nan??')
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% NOw ESTIAMTE EACH PARTNER PHASE LOCKING %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for sn=1:size(partner_sessions,1)

            this_session_time       = find(lfp_time>=partner_sessions.Start(sn) & lfp_time<=partner_sessions.End(sn));
            this_session_peaks      = peaksforalignment(ismember(peaksforalignment,this_session_time));
            

            %estimating each partner session recording  MVL, p_val,
            % prefered angle, mean rate and the corresdponign psth
            psth_this_session_freq  = nan(numel(this_session_peaks),numel(edges_freq)-1);

            for j=1:numel(this_session_peaks)    
                this_peak_time                  = lfp_time(this_session_peaks(j));         
                spikes_this_peak                = this_spikes(this_spikes>=this_peak_time+psth_range_freq(1) & this_spikes<=this_peak_time+psth_range_freq(2))-this_peak_time;
                psth_this_session_freq(j,:)    = histcounts(spikes_this_peak,edges_freq);
            end
            session_psth(sn,nn,:) = mean(psth_this_session_freq);


            this_session_spikes         = this_spikes(this_spikes>=partner_sessions.Start(sn) & this_spikes<=partner_sessions.End(sn));


            [~,~, mean_rate] = kernel_rate(this_session_spikes, partner_sessions.Start(sn), partner_sessions.End(sn), 0.01, 0.1);
            
            this_neuron_phases_freq     = interp1(lfp_time, phase_data,this_session_spikes);
            this_session_phases_freq    = phase_data(lfp_time>=partner_sessions.Start(sn) & lfp_time<=partner_sessions.End(sn));
            this_neuron_mvl             = circ_r(this_neuron_phases_freq);
            this_neuron_ppc             = compute_ppc_fast(this_neuron_phases_freq);


            this_neuron_prefered_angle  = circ_mean(this_neuron_phases_freq);

            rand_r      = nan(n_rand,1);
            rand_ppc    = nan(n_rand,1);

            for nr=1:n_rand

                rand_r(nr) =circ_r((randsample(this_session_phases_freq,numel(this_neuron_phases_freq), false))');
                rand_ppc(nr) =compute_ppc_fast((randsample(this_session_phases_freq,numel(this_neuron_phases_freq), false))');
            end
        
            this_neuron_mvl_p = sum(rand_r>this_neuron_mvl)/n_rand;
            this_neuron_ppc_p       = sum(rand_ppc>this_neuron_ppc)/n_rand;

            if isempty(this_neuron_prefered_angle),this_neuron_prefered_angle = NaN; end
            if isempty(this_neuron_ppc), this_neuron_ppc = NaN; end
            if isempty(this_neuron_ppc_p), this_neuron_ppc_p = NaN; end
            if isempty(this_neuron_mvl_p),this_neuron_mvl_p = NaN; end
            if isempty(this_neuron_mvl), this_neuron_mvl = NaN;end
            if isempty(mean_rate), mean_rate = NaN; end

            session_phase_stats(sn,nn,:) = [this_neuron_prefered_angle this_neuron_mvl this_neuron_mvl_p this_neuron_ppc this_neuron_ppc_p mean_rate this_id];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% NOW ESTIAMTE EACH PLAYBOUT PHASE LOCKING %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        mask = false(size(this_spikes));
        for b = 1:size(play_bouts_table,1)
            mask = mask | (this_spikes > play_bouts_table(b,1) & this_spikes < play_bouts_table(b,2));
        end
        play_spikes     = this_spikes(mask);

         mask = false(size(this_spikes));
        for b = 1:size(pre_play_bouts,1)
            mask = mask | (this_spikes > pre_play_bouts(b,1) & this_spikes < pre_play_bouts(b,2));
        end
        pre_paly_spikes = this_spikes(mask);

        play_phases     = interp1(lfp_time, phase_data,play_spikes);
        pre_paly_phases = interp1(lfp_time, phase_data,pre_paly_spikes);

        %estimating entire recording  MVL, p_val, prefered angle, mean rate
        %and the corresdponign psth
        play_mvl                 = circ_r(play_phases);
        play_prefered_angle      = circ_mean(play_phases);
        play_ppc                 = compute_ppc_fast(play_phases);
        tic
        rand_r            = nan(n_rand,1);
        rand_ppc          = nan(n_rand,1);
        for nr=1:n_rand
            rand_r(nr) =circ_r((randsample(phase_data,numel(play_phases), false))');
            rand_ppc(nr) = compute_ppc_fast((randsample(phase_data,numel(play_phases), false))');
        end
        toc

        play_mvl_p       = sum(rand_r>play_mvl)/n_rand;
        play_ppc_p       = sum(rand_ppc>play_ppc)/n_rand;
        mask = false(size(dt));
        for b = 1:size(play_bouts_table,1)
            mask = mask | (dt > play_bouts_table(b,1) & dt < play_bouts_table(b,2));
        end
        play_rate = mean(lambda( mask));


        pre_play_mvl                 = circ_r(pre_paly_phases);
        pre_play_prefered_angle      = circ_mean(pre_paly_phases);
        pre_play_ppc                 = compute_ppc_fast(pre_paly_phases);
        tic
        rand_r            = nan(n_rand,1);
        rand_ppc          = nan(n_rand,1);
        for nr=1:n_rand
            rand_r(nr) =circ_r((randsample(phase_data,numel(pre_paly_phases), false))');
            rand_ppc(nr) = compute_ppc_fast((randsample(phase_data,numel(pre_paly_phases), false))');
        end
        toc
        pre_play_mvl_p       = sum(rand_r>pre_play_mvl)/n_rand;
        pre_play_ppc_p       = sum(rand_ppc>pre_play_ppc)/n_rand;
         
        mask = false(size(dt));
        for b = 1:size(pre_play_bouts,1)
            mask = mask | (dt > pre_play_bouts(b,1) & dt < pre_play_bouts(b,2));
        end
        pre_play_rate        = mean(lambda(mask));



        if isempty(pre_play_prefered_angle), pre_play_prefered_angle = NaN; end
        if isempty(pre_play_ppc),pre_play_ppc = NaN; end
        if isempty(pre_play_ppc_p),pre_play_ppc_p = NaN; end
        if isempty(pre_play_mvl_p), pre_play_mvl_p = NaN; end
        if isempty(pre_play_mvl), pre_play_mvl = NaN; end
        if isempty(pre_play_rate), pre_play_rate = NaN; end

        if isempty(play_prefered_angle), play_prefered_angle = NaN; end
        if isempty(play_ppc),play_ppc = NaN; end
        if isempty(play_ppc_p),play_ppc_p = NaN; end
        if isempty(play_mvl_p), play_mvl_p = NaN; end
        if isempty(play_mvl), play_mvl = NaN; end
        if isempty(play_rate), play_rate = NaN; end


        play_phase_stats(1,nn,:)        = [    play_prefered_angle     play_mvl     play_mvl_p     play_ppc     play_ppc_p     play_rate this_id];
        pre_play_phase_stats(1,nn,:)    = [pre_play_prefered_angle pre_play_mvl pre_play_mvl_p pre_play_ppc pre_play_ppc_p pre_play_rate this_id];
        
        mask = false(size(lfp_time));
        for b = 1:size(play_bouts_table,1)
            mask = mask | (lfp_time > play_bouts_table(b,1) & lfp_time < play_bouts_table(b,2));
        end
        play_time                       = find(mask);
        play_peaks                      = peaksforalignment(ismember(peaksforalignment,play_time));

        psth_this_session_freq  = nan(numel(play_peaks),numel(edges_freq)-1);
        for j=1:numel(play_peaks)
            this_peak_time                  = lfp_time(play_peaks(j));
            spikes_this_peak                = this_spikes(this_spikes>=this_peak_time+psth_range_freq(1) & this_spikes<=this_peak_time+psth_range_freq(2))-this_peak_time;
            psth_this_session_freq(j,:)    = histcounts(spikes_this_peak,edges_freq);
        end
        play_psth (sn,nn,:) =  mean(psth_this_session_freq);

    
        mask = false(size(lfp_time));
        for b = 1:size(pre_play_bouts,1)
            mask = mask | (lfp_time > pre_play_bouts(b,1) & lfp_time < pre_play_bouts(b,2));
        end
        pre_play_time       = find(mask);
        pre_play_peaks      = peaksforalignment(ismember(peaksforalignment,pre_play_time));
        psth_this_session_freq  = nan(numel(pre_play_peaks),numel(edges_freq)-1);
        for j=1:numel(pre_play_peaks)
            this_peak_time                  = lfp_time(pre_play_peaks(j));
            spikes_this_peak                = this_spikes(this_spikes>=this_peak_time+psth_range_freq(1) & this_spikes<=this_peak_time+psth_range_freq(2))-this_peak_time;
            psth_this_session_freq(j,:)    = histcounts(spikes_this_peak,edges_freq);
        end
        pre_play_psth (sn,nn,:) =  mean(psth_this_session_freq);



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% NOW ESTIAMTE EACH NON-PLAYBOUT PHASE LOCKING %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

        mask = false(size(this_spikes));
        for b = 1:size(non_play_bouts_table,1)
            mask = mask | (this_spikes > non_play_bouts_table(b,1) & this_spikes < non_play_bouts_table(b,2));
        end
        non_play_spikes     = this_spikes(mask);

         mask = false(size(this_spikes));
        for b = 1:size(pre_non_play_bouts,1)
            mask = mask | (this_spikes > pre_non_play_bouts(b,1) & this_spikes < pre_non_play_bouts(b,2));
        end
        pre_paly_spikes = this_spikes(mask);

        non_play_phases     = interp1(lfp_time, phase_data,non_play_spikes);
        pre_paly_phases = interp1(lfp_time, phase_data,pre_paly_spikes);

        %estimating entire recording  MVL, p_val, prefered angle, mean rate
        %and the corresdponign psth
        non_play_mvl                 = circ_r(non_play_phases);
        non_play_prefered_angle      = circ_mean(non_play_phases);
        non_play_ppc                 = compute_ppc_fast(non_play_phases);
        tic
        rand_r            = nan(n_rand,1);
        rand_ppc          = nan(n_rand,1);
        for nr=1:n_rand
            rand_r(nr) =circ_r((randsample(phase_data,numel(non_play_phases), false))');
            rand_ppc(nr) = compute_ppc_fast((randsample(phase_data,numel(non_play_phases), false))');
        end
        toc

        non_play_mvl_p       = sum(rand_r>non_play_mvl)/n_rand;
        non_play_ppc_p       = sum(rand_ppc>non_play_ppc)/n_rand;
        mask = false(size(dt));
        for b = 1:size(non_play_bouts_table,1)
            mask = mask | (dt > non_play_bouts_table(b,1) & dt < non_play_bouts_table(b,2));
        end
        non_play_rate = mean(lambda( mask));


        pre_non_play_mvl                 = circ_r(pre_paly_phases);
        pre_non_play_prefered_angle      = circ_mean(pre_paly_phases);
        pre_non_play_ppc                 = compute_ppc_fast(pre_paly_phases);
        tic
        rand_r            = nan(n_rand,1);
        rand_ppc          = nan(n_rand,1);
        for nr=1:n_rand
            rand_r(nr) =circ_r((randsample(phase_data,numel(pre_paly_phases), false))');
            rand_ppc(nr) = compute_ppc_fast((randsample(phase_data,numel(pre_paly_phases), false))');
        end
        toc
        pre_non_play_mvl_p       = sum(rand_r>pre_non_play_mvl)/n_rand;
        pre_non_play_ppc_p       = sum(rand_ppc>pre_non_play_ppc)/n_rand;
         
        mask = false(size(dt));
        for b = 1:size(pre_non_play_bouts,1)
            mask = mask | (dt > pre_non_play_bouts(b,1) & dt < pre_non_play_bouts(b,2));
        end
        pre_non_play_rate        = mean(lambda(mask));



        if isempty(pre_non_play_prefered_angle), pre_non_play_prefered_angle = NaN; end
        if isempty(pre_non_play_ppc),pre_non_play_ppc = NaN; end
        if isempty(pre_non_play_ppc_p),pre_non_play_ppc_p = NaN; end
        if isempty(pre_non_play_mvl_p), pre_non_play_mvl_p = NaN; end
        if isempty(pre_non_play_mvl), pre_non_play_mvl = NaN; end
        if isempty(pre_non_play_rate), pre_non_play_rate = NaN; end

        if isempty(non_play_prefered_angle), non_play_prefered_angle = NaN; end
        if isempty(non_play_ppc),non_play_ppc = NaN; end
        if isempty(non_play_ppc_p),non_play_ppc_p = NaN; end
        if isempty(non_play_mvl_p), non_play_mvl_p = NaN; end
        if isempty(non_play_mvl), non_play_mvl = NaN; end
        if isempty(non_play_rate), non_play_rate = NaN; end


        non_play_phase_stats(1,nn,:)        = [    non_play_prefered_angle     non_play_mvl     non_play_mvl_p     non_play_ppc     non_play_ppc_p     non_play_rate this_id];
        pre_non_play_phase_stats(1,nn,:)    = [pre_non_play_prefered_angle pre_non_play_mvl pre_non_play_mvl_p pre_non_play_ppc pre_non_play_ppc_p pre_non_play_rate this_id];
        
        mask = false(size(lfp_time));
        for b = 1:size(non_play_bouts_table,1)
            mask = mask | (lfp_time > non_play_bouts_table(b,1) & lfp_time < non_play_bouts_table(b,2));
        end
        non_play_time                       = find(mask);
        non_play_peaks                      = peaksforalignment(ismember(peaksforalignment,non_play_time));

        psth_this_session_freq  = nan(numel(non_play_peaks),numel(edges_freq)-1);
        for j=1:numel(non_play_peaks)
            this_peak_time                  = lfp_time(non_play_peaks(j));
            spikes_this_peak                = this_spikes(this_spikes>=this_peak_time+psth_range_freq(1) & this_spikes<=this_peak_time+psth_range_freq(2))-this_peak_time;
            psth_this_session_freq(j,:)    = histcounts(spikes_this_peak,edges_freq);
        end
        non_play_psth (sn,nn,:) =  mean(psth_this_session_freq);

    
        mask = false(size(lfp_time));
        for b = 1:size(pre_non_play_bouts,1)
            mask = mask | (lfp_time > pre_non_play_bouts(b,1) & lfp_time < pre_non_play_bouts(b,2));
        end
        pre_non_play_time       = find(mask);
        pre_non_play_peaks      = peaksforalignment(ismember(peaksforalignment,pre_non_play_time));
        psth_this_session_freq  = nan(numel(pre_non_play_peaks),numel(edges_freq)-1);
        for j=1:numel(pre_non_play_peaks)
            this_peak_time                  = lfp_time(pre_non_play_peaks(j));
            spikes_this_peak                = this_spikes(this_spikes>=this_peak_time+psth_range_freq(1) & this_spikes<=this_peak_time+psth_range_freq(2))-this_peak_time;
            psth_this_session_freq(j,:)    = histcounts(spikes_this_peak,edges_freq);
        end
        pre_non_play_psth (sn,nn,:) =  mean(psth_this_session_freq);

    end


end
    
if any(~assigned)

    disp('Not all neurons assigned')
end
phase_struct.cluster_info                       = cluster_info;
%store partners session data

phase_struct.session_phase_stats                = session_phase_stats;
phase_struct.session_psth                       = session_psth;


%store entire session data
phase_struct.entire_recording_phase_stats       = entire_recording_phase_stats;
phase_struct.entire_recording_psth              = entire_recording_psth;

%store  play data
phase_struct.play_phase_stats                   = play_phase_stats;
phase_struct.play_psth                          = play_psth;

%store pre-play data
phase_struct.pre_play_phase_stats               = pre_play_phase_stats;
phase_struct.pre_play_psth                      = pre_play_psth;


%store  non-play data
phase_struct.non_play_phase_stats                   = non_play_phase_stats;
phase_struct.non_play_psth                          = non_play_psth;

%store pre-non-play data
phase_struct.pre_non_play_phase_stats               = pre_non_play_phase_stats;
phase_struct.pre_non_play_psth                      = pre_non_play_psth;



phase_struct.psth_range_freq                    = psth_range_freq;
phase_struct.partner_sessions                   = partner_sessions;
phase_struct.edges_freq                         = edges_freq;
phase_struct.areas_by_channel                   = areas_by_channel;
phase_struct.channel_map                        = channel_map;
phase_struct.psth_alignment                     = psth_alignment;

end
