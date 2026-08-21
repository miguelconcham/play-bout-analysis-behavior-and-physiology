function coupling_struct = GENERATE_FREQ_COUPLING_STRUCT(current_dir,Hd_freq1,Hd_freq2 )


% synch_directory     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Synch data';
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
%% select_mid_pag_channel
disp('Loading Channel Map')
hard_coded_x_coords = [8 40;258 290; 508 540; 758 790];
area_limit = readtable(area_limit_table);
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
this_animal = ['Batch', animal_batch(2), repeated_animal,animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName,this_animal),:);

if NPX_Type == 1

    PAG_channels = area_limit{ismember(area_limit.area, {'LPAG'}), {'ch_start', 'ch_end'}};
    PAG_channels = str2double(PAG_channels);
    channel_Range = [min(PAG_channels(:)) max(PAG_channels(:))];
    mid_PAG_channel = round(mean(channel_Range));
else
    load([current_dir,'\ChannelMap.mat'], 'xcoords', 'ycoords','chanMap')
    Y_Range = area_limit{ismember(area_limit.area, {'LPAG'}), {'ProbeNum','depth_start', 'depth_end'}};
    mid_PAG_channel = nan(size(Y_Range,1),1);
    figure
    plot(xcoords,ycoords, 'k.')
    hold on


    for j=1:size(Y_Range,1)
        this_indexes = ycoords>=Y_Range(j,2) & ycoords<=Y_Range(j,3) & ismember(xcoords,hard_coded_x_coords(Y_Range(j,1),:));

        all_locs = [xcoords(this_indexes) ycoords(this_indexes)];
        plot(all_locs(:,1),all_locs(:,2), 'r.')

        mean_loc = mean(all_locs);
        [~, closest_channel]= min(sum(abs([xcoords ycoords]-repmat(mean_loc,numel(ycoords),1)),2));
        plot(xcoords(closest_channel), ycoords(closest_channel), 'xb')
        mid_PAG_channel(j) = chanMap(closest_channel);

    end

end
%% obtain_psth



PAG_LFP         = double(LFP(mid_PAG_channel,:));
clear LFP
sr_LFP = 2500;


coupling_struct     = [];
for ch_n=1:numel(mid_PAG_channel)
    disp('Filtering Signal')

    filtered_signal_freq1 = filtfilt(Hd_freq1.Coefficients, 1, PAG_LFP(ch_n,:));
    hiblert_data1 = hilbert(filtered_signal_freq1);

    filtered_signal_freq2 = filtfilt(Hd_freq2.Coefficients, 1, PAG_LFP(ch_n,:));
    hiblert_data2 = hilbert(filtered_signal_freq2);


    phase_data1         = angle(hiblert_data1);
    amplitud_data1      = abs(hiblert_data1);
    amplitud_data2      = abs(hiblert_data2);
    std_amp1 = std(amplitud_data1);
    [~,max_locs] = findpeaks(amplitud_data1, 'MinPeakProminence',.5*std_amp1, 'MinPeakDistance', sr_LFP/(Hd_freq1.CutoffFrequency2  )) ;


    original_distribution = phase_data1(max_locs);

    figure
    subplot(1,2,1)
    polarhistogram(original_distribution, linspace(-pi, pi, 180))
    title([animal_code, 'Original'])
    mean_angle_original = angle(mean(exp(1i*original_distribution)));
    mean_angle_original =  mod( mean_angle_original+2*pi,2*pi);

    current_rho = rlim;
    hold on
    polarplot( [mean_angle_original mean_angle_original],current_rho, 'r')
    mean_angle_deg =mean_angle_original*180/pi;

    disp(['Rotating angles in ', num2str(mean_angle_deg),  ' degrees'])
    phase_data1 = mod(phase_data1  - mean_angle_original + 5*pi , 2*pi) - pi; %% centering step



    subplot(1,2,2)
    polarhistogram(phase_data1(max_locs),linspace(-pi, pi, 180))
    current_rho = rlim;
    mean_angle_new = angle(mean(exp(1i*phase_data1(max_locs))));
    hold on
    polarplot( [mean_angle_new mean_angle_new],current_rho, 'r')
    title('Corrected')
    pause(.1)
    
    centered_distribution = phase_data1(max_locs);
    % 
    % min_val = min(phase_data1);
    % phase_data1(phase_data1<0) = pi*phase_data1(phase_data1<0)/abs(min_val); %extending step1: negative values
    % max_val = max(phase_data1);
    % phase_data1(phase_data1>0) = pi*phase_data1(phase_data1>0)/abs(max_val); %extending step2: positive values
    % 
    % 
    % 
    % subplot(2,2,3)
    % polarhistogram(phase_data1(max_locs), linspace(-pi, pi, 180))
    % title([animal_code, 'Extended'])
    % mean_angle_original = angle(mean(exp(1i*phase_data1(max_locs))));
    % mean_angle_original =  mod( mean_angle_original+2*pi,2*pi);
    % 
    % current_rho = rlim;
    % hold on
    % polarplot( [mean_angle_original mean_angle_original],current_rho, 'r')
    % mean_angle_deg =mean_angle_original*180/pi;
    % 
    % disp(['Rotating angles in ', num2str(mean_angle_deg),  ' degrees'])
    % phase_data1 = mod(phase_data1  - mean_angle_original + 5*pi , 2*pi) - pi;
    % 
    % 
    % 
    % subplot(2,2,4)
    % polarhistogram(phase_data1(max_locs),linspace(-pi, pi, 180))
    % current_rho = rlim;
    % mean_angle_new = angle(mean(exp(1i*phase_data1(max_locs))));
    % hold on
    % polarplot( [mean_angle_new mean_angle_new],current_rho, 'r')
    % title('Extended and recentered')
    % 
    % extended_distribution =  phase_data1(max_locs);
    % pause(.1)



    % Inputs:
    % phase_data_1      -> vector of phases (radians, -pi to pi)
    % amplitud_data_2   -> vector of amplitudes (same length)

    nBins = 72;       % number of phase bins
    nPerm = 500;     % number of shuffles for null distribution

    % --- Step 1: real data binning ---
    edges = linspace(-pi, pi, nBins+1);
    [~,~,binIdx] = histcounts(phase_data1, edges);

    meanAmp_real = zeros(1,nBins);
    for b = 1:nBins
        meanAmp_real(b) = mean(amplitud_data2(binIdx == b), 'omitnan');
    end

    % Normalize to get a probability-like distribution
    p_real = meanAmp_real / sum(meanAmp_real, 'omitmissing');

    % Compute Modulation Index (MI) via KL divergence
    MI_real = (log(nBins) + sum(p_real .* log(p_real+eps))) / log(nBins);
    [r_real, ~] = corr(amplitud_data1(:), amplitud_data2(:), 'Type', 'Pearson');

    % --- Step 2: permutation test (shuffle phase indices) ---
    MI_perm = zeros(nPerm,1);
    N = numel(phase_data1);
    
    disp('Estimating Permuted distribution')
    all_mean_perm = nan(nPerm,nBins);
     r_perm = zeros(nPerm,1);
    for p = 1:nPerm
        if mod(p, 100)==0
            disp(['Iteration #', num2str(p)])
        end
        shuffled_idx = randperm(N);
        phase_shuffled = phase_data1(shuffled_idx);

        [~,~,binIdx] = histcounts(phase_shuffled, edges);
        meanAmp_perm = zeros(1,nBins);
        
        for b = 1:nBins
            meanAmp_perm(b) = mean(amplitud_data2(binIdx == b), 'omitnan');
        end
        all_mean_perm(p,:) = meanAmp_perm;
        p_perm = meanAmp_perm / sum(meanAmp_perm, 'omitmissing');
        MI_perm(p) = (log(nBins) + sum(p_perm .* log(p_perm+eps))) / log(nBins);

        shuffled_idx = randperm(N);
        r_perm(p) = corr(amplitud_data1(:), amplitud_data2(shuffled_idx)', 'Type', 'Pearson');

    end
   
    % --- Step 3: significance test ---
    p_val_MI    = mean(MI_perm >= MI_real);
    p_val_corr  = mean(r_perm >= r_real);

    % --- Step 4: obtain mean and std of amplitud to quantify effect ---
    mean_ampl   = mean(amplitud_data2);
    std_ampl    = std(amplitud_data2);

    
    if ch_n==1
        coupling_struct.MI_real             = MI_real;
        coupling_struct.MI_perm             = MI_perm;
        coupling_struct.meanAmp_real        = meanAmp_real;
        coupling_struct.all_mean_perm       = all_mean_perm;
        coupling_struct.r_real              = r_real;
        coupling_struct.r_perm              = r_perm;
        coupling_struct.p_val_MI            = p_val_MI;
        coupling_struct.p_val_corr          = p_val_corr;
        coupling_struct.Hd_freq1            = Hd_freq1;
        coupling_struct.Hd_freq2            = Hd_freq2;
        coupling_struct.mean_ampl           = mean_ampl;
        coupling_struct.std_ampl            = std_ampl;
        coupling_struct.edges               = edges;
        coupling_struct.original_distr      = original_distribution;
        coupling_struct.centered_distr      = centered_distribution;
     
       


    else
        coupling_struct(ch_n).MI_real               = MI_real;
        coupling_struct(ch_n).MI_perm               = MI_perm;
        coupling_struct(ch_n).meanAmp_real          = meanAmp_real;
        coupling_struct(ch_n).all_mean_perm         = all_mean_perm;
        coupling_struct(ch_n).r_real                = r_real;
        coupling_struct(ch_n).r_perm                = r_perm;
        coupling_struct(ch_n).p_val_MI              = p_val_MI;
        coupling_struct(ch_n).p_val_corr            = p_val_corr;
        coupling_struct(ch_n).Hd_freq1              = Hd_freq1;
        coupling_struct(ch_n).Hd_freq2              = Hd_freq2;
        coupling_struct(ch_n).mean_ampl             = mean_ampl;
        coupling_struct(ch_n).std_ampl              = std_ampl;
        coupling_struct(ch_n).edges                 = edges;
        coupling_struct(ch_n).original_distr        = original_distribution;
        coupling_struct(ch_n).centered_distr        = centered_distribution;



    end

end
end
