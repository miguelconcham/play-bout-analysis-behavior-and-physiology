list_of_animals = {'B1D1 1013 Dual','B1S3 1008 Single','B1S3 1009 Single','B2S2 1110 Single2','B2S2 1111 Single2','B3D2 1130 Dual','B4S2 0825 Single'};


saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';
f          = 4:.1:15;      % frequency range for spectrogram
freq_range = [6 12]; 

% f          = .1:.05:6;      % frequency range for spectrogram
% freq_range = [.1 5]; 

f          = 35:1:100;      % frequency range for spectrogram
freq_range = [35 90]; 


wind_params = [.250 .245];
n_strctut = 1;

psth_structure = [];
animal_names = [];
%%
for fn = 1:numel(list_of_animals)

    if fn==1
        psth_structure =GENERATE_THETA_PSTH_CALLS(list_of_animals{fn}, f, freq_range,wind_params);
        n_strctut = n_strctut+numel(psth_structure);
        animal_names = [animal_names;[repmat({list_of_animals{fn}},numel(psth_structure),1), num2cell(1:numel(psth_structure))']]
    else

        transt_psth = GENERATE_THETA_PSTH_CALLS(list_of_animals{fn}, f, freq_range,wind_params);
      
        for sub_j=1:numel(transt_psth)
    
            psth_structure(n_strctut) = transt_psth(sub_j);
            n_strctut = n_strctut+1;
        end
        animal_names = [animal_names;[repmat({list_of_animals{fn}},numel(transt_psth),1) num2cell(1:numel(transt_psth))' ]]

    end


end

%% saving
disp('saving')
% save([saving_folder,'\psth_structure_call_theta.mat'],'psth_structure');
% save([saving_folder,'\animal_names_call_theta.mat'],'animal_names');
% save([saving_folder,'\psth_structure_call_gamma.mat'],'psth_structure');
% save([saving_folder,'\animal_names_call_gamma.mat'],'animal_names');
%%
list_of_animals = {'B1D1 1013 Dual','B1S3 1008 Single','B1S3 1009 Single','B2S2 1110 Single2','B2S2 1111 Single2','B3D2 1130 Dual','B4S2 0825 Single'};


saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';
f          = 4:.1:15;      % frequency range for spectrogram
freq_range = [6 12]; 

f          = .1:.05:6;      % frequency range for spectrogram
freq_range = [.1 5]; 
wind_params = [.250 .245];
n_strctut = 1;
% 
% psth_structure = [];
% animal_names = [];

load([saving_folder,'\psth_structure_call_delta.mat'],'psth_structure');
load([saving_folder,'\animal_names_call_delta.mat'],'animal_names');
%% mergin data

smooth_wind             = 20;
baseline_range          = [-2 0];
animal_label            = {'B1D1','B1S3','B2S2','B3D2', 'B4S2'};
electorde_number        = [1 2];
bin_size                = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
psth_ranges             = psth_structure(1).hist_range;
time                    = psth_ranges(1):bin_size:psth_ranges(2)+bin_size;
baseline_index          = time<baseline_range(2) & time>baseline_range(1);
all_psth_onset          = [];
all_psth_offset         = [];
all_onset_regressors    = [];
all_offset_regressors   = [];
all_Calls               = [];

electrode_index  = [];
animal_index = [];
for j=1:numel(psth_structure)

    if contains(animal_names{j},animal_label)
        
        animal_num      = find(cell2mat(cellfun(@(x) contains(animal_names{j},x), animal_label, 'UniformOutput',false)));  
        electrode_num   = animal_names{j,2};

        this_psth_onset         = psth_structure(j).call_onset;
          animal_index = [animal_index;repmat(animal_num,size(this_psth_onset,1),1)];
          electrode_index = [electrode_index;ones(size(this_psth_onset,1),1)*electrode_num];
        for trial=1:size(this_psth_onset,1)
            this_psth_onset(trial,:) = ( this_psth_onset(trial,:) - mean( this_psth_onset(trial,baseline_index)))/std( this_psth_onset(trial,baseline_index));
            this_psth_onset(trial,:) = movmean(this_psth_onset(trial,:), smooth_wind);
        end
        all_psth_onset      = [all_psth_onset; this_psth_onset];

        this_psth_onset     = psth_structure(j).call_onset;
        this_psth_offset    = psth_structure(j).call_offset;
         for trial=1:size(this_psth_offset,1)
            this_psth_offset(trial,:) = ( this_psth_offset(trial,:) - mean( this_psth_onset(trial,baseline_index)))/std( this_psth_onset(trial,baseline_index));
            this_psth_offset(trial,:) = movmean(this_psth_offset(trial,:), smooth_wind);
         end         
        all_psth_offset = [all_psth_offset; this_psth_offset];

      
        all_onset_regressors = [all_onset_regressors; psth_structure(j).call_onset_regressor];
        all_offset_regressors = [all_offset_regressors; psth_structure(j).call_onset_regressor];
    

        all_Calls = [all_Calls;psth_structure(j).CallStats];
    end
end

call_lengths = all_Calls.CallLengths;


[sorted_call_lengths, order] = sort(call_lengths);

%%
X_lim = [-.5 1];
figure
min_length = .0;
respones_array = [];
for an= 1:numel(animal_label)
    animal_bool = animal_index==an;
    length_bool = call_lengths>min_length;
    electrode_bool = electrode_index==1;
    [sorted_call_lengths, order] = sort(call_lengths(animal_bool & length_bool & electrode_bool,:));
    subplot(5,numel(animal_label),(1:numel(animal_label):2*numel(animal_label)) + an-1)

    array = all_psth_onset(animal_bool & length_bool & electrode_bool,:);
    imagesc(time,1:numel(sorted_call_lengths),array(order,:) )
    xlim(X_lim)
    clim([-2 2])
    axis xy
    hold on
    plot([0 0],[1 numel(sorted_call_lengths)], 'w')
    plot(sorted_call_lengths,1:numel(sorted_call_lengths), 'w')
    title(animal_label{an})

    subplot(5,numel(animal_label),((2*numel(animal_label) + 1):numel(animal_label):5*numel(animal_label)) + an-1)

    [~, ~, ci]  = ttest(array);
    no_nan = ~any(isnan(ci));
    fill([time(no_nan) fliplr(time(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
    hold on
    plot(time,trimmean(array,5), 'k')
    yyaxis right
    plot(time, mean(all_onset_regressors(animal_bool & length_bool,:)), 'r')
    respones_array = [respones_array;trimmean(array,5)];
    xlim(X_lim)
end


%% all session together, no mxied model
figure
X_lim = [-2 2]
min_length = 0.05
 length_bool = call_lengths>min_length;
    length_bool = call_lengths>min_length;
    electrode_bool = electrode_index==1;
    [sorted_call_lengths, order] = sort(call_lengths( length_bool & electrode_bool,:));
    subplot(5,1,1:3)

    array = all_psth_onset( length_bool & electrode_bool,:);
    imagesc(time,1:numel(sorted_call_lengths),array(order,:) )
    xlim(X_lim)
    clim([-2 2])
    axis xy
    hold on
    plot([0 0],[1 numel(sorted_call_lengths)], 'w')
    plot(sorted_call_lengths,1:numel(sorted_call_lengths), 'w')
    title(animal_label{an})

    subplot(5,1,4:5)

    [~, ~, ci]  = ttest(array);
    no_nan = ~any(isnan(ci));
    fill([time(no_nan) fliplr(time(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
    hold on
    plot(time,mean(array,'omitmissing'), 'k')
    yyaxis right
    plot(time, mean(all_onset_regressors(animal_bool & length_bool,:)), 'r')
    xlim(X_lim)

%% estimating mixed model per time

length_limit = .0;

indextinclude = call_lengths>=length_limit;
power = all_psth_onset(indextinclude,:);
% power = all_psth_tw_3points(indextinclude,:);

subject_idx =animal_index(indextinclude);
time_range = [-1 2];
limted_time = time;
valid_bins = find(limted_time >= time_range(1) & limted_time <= time_range(2));
power = power(:,valid_bins);
limted_time = limted_time(valid_bins);

[nTrials, nTime] = size(power);

% Preallocate
est  = nan(nTime,1);
se   = nan(nTime,1);
pvals = nan(nTime,1);
d    = nan(nTime,1);
ci   = nan(nTime,2);

% Reshape to long format
[trial_idx, time_idx] = ndgrid(1:nTrials,1:nTime);
tbl = table;
tbl.Power = power(:);

data2include = abs(tbl.Power )<2.5;
tbl.Subject = categorical(subject_idx(trial_idx(:)));
time_matrix = repmat(limted_time, nTrials, 1);
tbl.Time = time_matrix(:);
tbl = tbl(data2include,:)
% Loop over bins
for i = 1:nTime
    tbl_t = tbl(tbl.Time == limted_time(i),:);
    if sum(~isnan(tbl_t.Power))>10
    lme = fitlme(tbl_t,'Power ~ 1 + (1|Subject)');

    est(i) = lme.Coefficients.Estimate(1);
    se(i)  = lme.Coefficients.SE(1);
    pvals(i) = lme.Coefficients.pValue(1);
    ci_t = coefCI(lme);
    ci(i,:) = ci_t(1,:);

    % Cohen's d-like standardized effect size
    sd_within = std(tbl_t.Power, 'omitmissing');
    d(i) = est(i) / sd_within;
    end
end

% Multiple comparison correction (FDR)
pvals_fdr = mafdr(pvals,'BHFDR',true);

% Package results
results.time = limted_time(:);
results.est = est;
results.se = se;
results.ci = ci;
results.pvals = pvals;
results.pvals_fdr = pvals_fdr;
results.d = d;
%%
save([saving_folder,'\results_call_updated_gamma.mat'],'results','respones_array','time');

%%
load([saving_folder,'\results_call_updated_delta.mat'],'results');

%%


alpha = 0.05;
forced_x_lim = [-1 2];

limited_time = results.time;
est = results.est;
ci = results.ci;
pvals_fdr = results.pvals_fdr;
d = results.d;

figure;
subplot(2,1,1); hold on;
plot(time,respones_array, 'k:')

% Shaded CI
no_nan = ~any(isnan(ci'));
fill([limited_time(no_nan); flipud(limited_time(no_nan))], [ci(no_nan,1); flipud(ci(no_nan,2))], ...
    [0.8 0.8 1], 'EdgeColor','none','FaceAlpha',0.4);
plot(limited_time, est, 'b','LineWidth',2);

% Mark significant bins
sig_idx = pvals_fdr < alpha;
plot(limited_time(sig_idx), est(sig_idx), 'r*','MarkerSize',6);

ylabel('Mean Power');
title('Mixed-Effects Power (CI + FDR-corrected sig)');
grid on;
xlim(forced_x_lim)

subplot(2,1,2);
plot(limited_time, d, 'k','LineWidth',2);
ylabel('Effect size (Cohen''s d)');
xlabel('Time (s)');
grid on;
xlim(forced_x_lim)
