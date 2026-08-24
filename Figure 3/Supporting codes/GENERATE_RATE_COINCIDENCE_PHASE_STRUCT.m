function synch_structure = GENERATE_RATE_COINCIDENCE_PHASE_STRUCT( ...
    npx_data_dir, Hd_freq1, bin_size, hist_range, time_precision, areas2compare)

%% =========================
% PATHS / META
%% =========================

% synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
area2analyze        = 'PAG';
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\Area_limits_GoodLooking.xlsx';
behavior_data       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Behavior backups';
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

psth_edges = psth_range_freq(1):bin_size:psth_range_freq(2);


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





%% 2 Create play bout array

bin_size_pb                = 0.01;
conv_length             = 1;
animal_types            = unique(Behavior.Animal);

animal_types(ismember(animal_types,'Session_structure'))                =[];




config.Behavior         = Behavior;
config.repeated_animal  = repeated_animal;
config.animal_types     = animal_types        ;
config.play_behaviors   = play_behaviors      ;
config.beh_bin          = bin_size_pb             ;
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


%% =========================
% LFP + PHASE (YOUR REAL METHOD)
%% =========================



sr_LFP = 2500;
lfp_time = (1:size(LFP,2))/sr_LFP;

% only channels used by spikes
channels_with_spikes = unique(cluster_info.ch)';

PAG_LFP = double(LFP(channels_with_spikes,:));
clear LFP
sp_n =1;
figure('units','normalized','outerposition',[0 0 1 1]);
pause(.1)
filtered_signal_all_channels    = PAG_LFP;
amplitud_data_all_channels      = nan(size(PAG_LFP));
phase_data_all_channels         = nan(size(PAG_LFP));
phase_data = nan(size(PAG_LFP));
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

%% =========================
% PREALLOC OUTPUT
%% =========================
nPairs = nchoosek(size(cluster_info,1),2);
synch_comb = nchoosek(1:size(cluster_info,1),2);
nBouts = size(play_bouts_table,1);
nBins = numel(psth_edges)-1;

rate_n1   = nan(nPairs,nBouts,nBins);
rate_n2   = nan(nPairs,nBouts,nBins);

coin_n1   = nan(nPairs,nBouts,nBins);
coin_n2   = nan(nPairs,nBouts,nBins);

phase_n1  = nan(nPairs,nBouts,nBins);
phase_n2  = nan(nPairs,nBouts,nBins);

plv_n1    = nan(nPairs,nBouts,nBins);
plv_n2    = nan(nPairs,nBouts,nBins);

cluster_spike_times = cell(size(cluster_info,1),1);

for c=1:size(cluster_info,1)
    cluster_spike_times{c} = spike_times(spike_clusters == cluster_info.cluster_id(c));
end

%% =========================
% MAIN LOOP
%% =========================

for p=1:nPairs

    c1 = synch_comb(p,1);
    c2 = synch_comb(p,2);

    s1 = cluster_spike_times{c1};
    s2 = cluster_spike_times{c2};

    ch1 = cluster_info.ch(c1);
    ch2 = cluster_info.ch(c2);

    i_ch1 = find(channels_with_spikes == ch1);
    i_ch2 = find(channels_with_spikes == ch2);

    ph1 = phase_data(i_ch1,:);
    ph2 = phase_data(i_ch2,:);

    for pb=1:nBouts

        pb_start = play_bouts_table(pb,1);

        s1_pb = s1(s1>=pb_start+hist_range(1) & s1<=pb_start+hist_range(2));
        s2_pb = s2(s2>=pb_start+hist_range(1) & s2<=pb_start+hist_range(2));

        % =====================
        % FIRING RATE
        % =====================
        rate_n1(p,pb,:) = histcounts(s1_pb-pb_start,psth_edges)/bin_size;
        rate_n2(p,pb,:) = histcounts(s2_pb-pb_start,psth_edges)/bin_size;

        if isempty(s1_pb) || isempty(s2_pb)
            continue;
        end

        % =====================
        % COINCIDENCE
        % =====================
        [i1,i2] = find_synchronous_spikes_twopointer(s1_pb,s2_pb,time_precision);

        csp1 = s1_pb(i1);
        csp2 = s2_pb(i2);

        coin_n1(p,pb,:) = histcounts(csp1-pb_start,psth_edges)/bin_size;
        coin_n2(p,pb,:) = histcounts(csp2-pb_start,psth_edges)/bin_size;

        % =====================
        % PHASE PER BIN
        % =====================
        for b=1:nBins

            t0 = pb_start + psth_edges(b);
            t1 = pb_start + psth_edges(b+1);

            idx1 = csp1>=t0 & csp1<t1;
            idx2 = csp2>=t0 & csp2<t1;

            % ---- neuron 1 ----
            if any(idx1)
                t_spk = round(csp1(idx1)*sr_LFP);
                t_spk = max(1,min(t_spk,size(ph1,2)));

                ph = ph1(t_spk);

                phase_n1(p,pb,b) = angle(mean(exp(1i*ph)));
                plv_n1(p,pb,b)   = abs(mean(exp(1i*ph)));
            end

            % ---- neuron 2 ----
            if any(idx2)
                t_spk = round(csp2(idx2)*sr_LFP);
                t_spk = max(1,min(t_spk,size(ph2,2)));

                ph = ph2(t_spk);

                phase_n2(p,pb,b) = angle(mean(exp(1i*ph)));
                plv_n2(p,pb,b)   = abs(mean(exp(1i*ph)));
            end

        end
    end
end

%% =========================
% OUTPUT STRUCT (STACKABLE)
%% =========================

synch_structure.rate_n1  = rate_n1;
synch_structure.rate_n2  = rate_n2;

synch_structure.coin_n1  = coin_n1;
synch_structure.coin_n2  = coin_n2;

synch_structure.phase_n1 = phase_n1;
synch_structure.phase_n2 = phase_n2;

synch_structure.plv_n1   = plv_n1;
synch_structure.plv_n2   = plv_n2;

synch_structure.psth_edges = psth_edges;

synch_structure.cluster_info = cluster_info;
synch_structure.synch_comb   = synch_comb;

synch_structure.bin_size   = bin_size;
synch_structure.hist_range = hist_range;
synch_structure.time_precision = time_precision;



%% =========================
% SYNC DETECTOR
%% =========================


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