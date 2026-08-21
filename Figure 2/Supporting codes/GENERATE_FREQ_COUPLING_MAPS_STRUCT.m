function coupling_struct = GENERATE_FREQ_COUPLING_MAPS_STRUCT(current_dir,Hd_freq1,Hd_freq2 )


% synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
chan_map_folder     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\StarndarChannMap';
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\Area_limits_GoodLooking.xlsx';
% npx_raw_data = 
animal_code         = strsplit(current_dir, '\');
animal_code         = animal_code{end};
animal_code_params  = strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
repeated_animal     = animal_code_params{3};
%% define parameters


%% load synch from synch folder
% load([synch_directory,'\', animal_code, '\synch_model_video2NPX.mat'],'synch_model_video2NPX')




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
%% Create channel map
disp('Loading Channel Map')
areas_by_channel = cell(384,1);
channel_map      = nan(384,2);
hard_coded_x_coords = [8 40;258 290; 508 540; 758 790];
area_limit = readtable(area_limit_table);

% Build animal identifier for area selection
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal,animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName,this_animal),:);

if NPX_Type == 1
    load([chan_map_folder,'\neuropixPhase3A_kilosortChanMap.mat'], 'xcoords','ycoords', 'chanMap' )
    
 



  for ch_n=1:384
      ch = chanMap(ch_n)+1;
      channel_map(ch,1) = xcoords(ch_n);
      channel_map(ch,2) = ycoords(ch_n);
      areas_by_channel{ch} = area_limit.area{ycoords(ch_n)>=area_limit.depth_start &  ycoords(ch_n)<area_limit.depth_end+1 & ismember(area_limit.Probe_Area, 'PAG') };
  end

else
    load([current_dir,'\ChannelMap.mat'], 'xcoords', 'ycoords','chanMap')

  
    
   for ch_n=1:384
       probe_n = find(any(ismember(hard_coded_x_coords,xcoords(ch_n)),2));
      ch = chanMap(ch_n)+1;
      channel_map(ch,1) = xcoords(ch_n);
      channel_map(ch,2) = ycoords(ch_n);
      areas_by_channel{ch} = area_limit.area{ycoords(ch_n)>=area_limit.depth_start &  ycoords(ch_n)<area_limit.depth_end+1 & area_limit.ProbeNum==probe_n};
  end

      

end
%% obtain_psth



PAG_LFP         = double(LFP);
clear LFP
sr_LFP = 2500;
nBins = 72;       % number of phase bins

% --- Step 1: real data binning ---
edges               = linspace(-pi, pi, nBins+1);
meanAmp_real        = zeros(size(PAG_LFP,1),nBins);
coupling_struct     = [];
r_real              = nan(size(PAG_LFP,1),2);
power_hig_freq      = nan(size(PAG_LFP,1),2);
power_low_freq      = nan(size(PAG_LFP,1),2);
phase_diff          = nan(size(PAG_LFP,1),2);
MI_real             = nan(size(PAG_LFP,1),2);

for ch_n=1:size(PAG_LFP,1)
    if mod(ch_n,25)==1
        disp(['Processing ch #', num2str(ch_n)])
    end
    
   

    filtered_signal_freq1 = filtfilt(Hd_freq1.Coefficients, 1, PAG_LFP(ch_n,:));
    hiblert_data1 = hilbert(filtered_signal_freq1);

    filtered_signal_freq2 = filtfilt(Hd_freq2.Coefficients, 1, PAG_LFP(ch_n,:));
    hiblert_data2 = hilbert(filtered_signal_freq2);


        
    phase_data1         = angle(hiblert_data1);

    amplitud_data1      = abs(hiblert_data1);
    amplitud_data2      = abs(hiblert_data2);
    if ch_n==1
        reference_phase =  phase_data1;
    end


    power_hig_freq(ch_n)      = median(amplitud_data2);
    power_low_freq(ch_n)      = median(amplitud_data1);
    phase_diff(ch_n)          = angle(mean(exp(1i*phase_data1) - exp(1i*reference_phase)));

    std_amp1 = std(amplitud_data1);
    [~,max_locs] = findpeaks(amplitud_data1, 'MinPeakProminence',.5*std_amp1, 'MinPeakDistance', sr_LFP/(Hd_freq1.CutoffFrequency2  )) ;

    original_distribution   = phase_data1(max_locs);       
    mean_angle_original     = angle(mean(exp(1i*original_distribution)));
    mean_angle_original     = mod( mean_angle_original+2*pi,2*pi);  
    mean_angle_deg          = mean_angle_original*180/pi;
    phase_data1             = mod(phase_data1  - mean_angle_original + 5*pi , 2*pi) - pi; %% centering step
  

    % Inputs:
    % phase_data_1      -> vector of phases (radians, -pi to pi)
    % amplitud_data_2   -> vector of amplitudes (same length)

   
    [~,~,binIdx] = histcounts(phase_data1, edges);
    for b = 1:nBins
        meanAmp_real(ch_n,b) = mean(amplitud_data2(binIdx == b), 'omitnan');
    end

    % Normalize to get a probability-like distribution
    p_real= meanAmp_real(ch_n,:) / sum(meanAmp_real(ch_n,:), 'omitmissing');

    % Compute Modulation Index (MI) via KL divergence
    MI_real(ch_n) = (log(nBins) + sum(p_real .* log(p_real+eps))) / log(nBins);
    [r_real(ch_n), ~] = corr(amplitud_data1(:), amplitud_data2(:), 'Type', 'Pearson');

    % --- Step 2: permutation test (shuffle phase indices) ---
   
   
    % --- Step 3: significance test ---
  

  


end

        coupling_struct.MI_real             = MI_real;
        coupling_struct.meanAmp_real        = meanAmp_real;
        coupling_struct.r_real              = r_real;
        coupling_struct.Hd_freq1            = Hd_freq1;
        coupling_struct.Hd_freq2            = Hd_freq2;
        coupling_struct.power_hig_freq      = power_hig_freq;
        coupling_struct.power_low_freq      = power_low_freq;
        coupling_struct.edges               = edges;
        coupling_struct.phase_diff      = phase_diff;
        coupling_struct.channel_map         = channel_map;
        coupling_struct.areas_by_channel    = areas_by_channel;
         

      




end
