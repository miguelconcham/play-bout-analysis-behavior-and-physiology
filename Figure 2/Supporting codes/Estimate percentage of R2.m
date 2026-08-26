%% 1 LOAD DATA
% Regressor PSTH structs (~5 GB) live on DataSets Theta psth; not copied into repo Data.
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\Theta psth';
disp('Loading')
% load([saving_folder,'\psth_structure_all_regressors_beta_gamma.mat'],'psth_structure'); 
% load([saving_folder,'\animal_names_all_regressors_beta_gamma.mat'],'animal_names');

load([saving_folder,'\psth_structure_all_regressors_delta_theta.mat'],'psth_structure'); 
load([saving_folder,'\animal_names_all_regressors_delta.mat'],'animal_names');
disp('Loading done')
%% 2 stacking data


% for DELTA AND THETA
freq_range = {[1 5],[6 12]}; 
wind_length = {1, .250};
spect_bin_size = 0.005;

% % for BETA AND GAMMA
% freq_range = {[15 30],[40 90]}; 
% wind_length = {.150, .08};
% spect_bin_size = 0.005;

regressors_var_names = fieldnames(psth_structure(1).regressors);
% List of array variables to stack
varNames = { ...
    'play_bout_onset_low_freq', 'play_bout_onset_all_low_freq', ...
    'play_bout_onset_high_freq', 'play_bout_onset_all_high_freq', ...
    'play_bouts_table'};
vars2zscore = {'play_bout_onset_low_freq','play_bout_onset_all_low_freq',...
    'animal_distance_onset_regressor', 'animal_distance_offset_regressor', ...
    'self_acc_onset_regressor','self_acc_offset_regressor'};

regressors2stack = {'animal_angle_onset','animal_angle_speed_onset','animal_angle_acc_onset', ...
                    'animal_speed_onset','animal_accel_onset','animal_speed_kalman_onset',...
                	'animal_accel_kalman_onset','partner_angle_onset',...
                    'partner_angle_speed_onset','partner_angle_acc_onset','partner_speed_onset',...
                	'partner_accel_onset','partner_speed_kalman_onset','partner_accel_kalman_onset',...
                	'relative_distance_onset',	'relative_angle_onset',...
                    'relative_angle_speed_onset','relative_angle_acc_onset','relative_speed_onset',...
                    'relative_acc_onset','relative_acc_kalman_onset','relative_speed_kalman_onset',...
                	'play_bout_onset_regressor','call_onset_regressor','self_onset_regressor',...
                	'other_onset_regressor'};
partner_regressor = {'PartnerNumber'};

regressors2zscore = {'animal_angle_onset','animal_angle_speed_onset','animal_angle_acc_onset', ...
                    'animal_speed_onset','animal_accel_onset','animal_speed_kalman_onset',...
                	'animal_accel_kalman_onset','partner_pos_onset'	'partner_angle_onset',...
                    'partner_angle_speed_onset','partner_angle_acc_onset','partner_speed_onset',...
                	'partner_accel_onset','partner_speed_kalman_onset','partner_accel_kalman_onset',...
                	'relative_distance_onset',	'relative_angle_onset',...
                    'relative_angle_speed_onset','relative_angle_acc_onset','relative_speed_onset',...
                    'relative_acc_onset','relative_acc_kalman_onset','relative_speed_kalman_onset'};

% Initialize stacking container
stacked_data = struct();
for v = 1:numel(varNames)
    stacked_data.(varNames{v}) = [];
end
stacked_data.meta = {};  % {animal_name, partner_id, channel, play_length}
regressors = [];
for v = 1:numel(regressors2stack)
    regressors.(regressors2stack{v}) = [];
end
 regressors.(partner_regressor{1}) = [];
% Loop through each session
for i = 1:numel(psth_structure)
    S = psth_structure(i);

    % Get session info from parallel cell array
    animal_name = animal_names{i,1};
    partner_id  = animal_names{i,2};
    ch          = animal_names{i,3};

    % Compute play lengths
    play_lengths = S.play_bouts_table(:,2) - S.play_bouts_table(:,1);
    nRows = size(S.play_bout_onset_low_freq,1);

    % Stack all variables
    for v = 1:numel(varNames)

        if ismember(varNames{v}, vars2zscore)
            array = S.(varNames{v});
            if ndims(array)==2
            array = (array - repmat(mean(S.(varNames{v})(:), 'omitmissing'), size(array,1), size(array,2)))./repmat(std(S.(varNames{v})(:), 'omitmissing'), size(array,1), size(array,2));
            stacked_data.(varNames{v}) = [stacked_data.(varNames{v}); array];
            else
                array = log10(array);
                array = (array - repmat(mean(array(:), 'omitmissing'), size(array,1), size(array,2),size(array,3)))./repmat(std(array(:), 'omitmissing'), size(array,1), size(array,2),size(array,3));
                stacked_data.(varNames{v}) = [stacked_data.(varNames{v}); array];
            end

        else
            stacked_data.(varNames{v}) = [stacked_data.(varNames{v}); S.(varNames{v})];
        end
    end
    for v = 1:numel(regressors2stack)
        if ismember(regressors2stack{v}, regressors2zscore)
            array = S.regressors.(regressors2stack{v});
            
            array = (array - repmat(mean(S.regressors.(regressors2stack{v})(:), 'omitmissing'), size(array,1), size(array,2)))./repmat(std(S.regressors.(regressors2stack{v})(:), 'omitmissing'), size(array,1), size(array,2));
            regressors.(regressors2stack{v}) = [regressors.(regressors2stack{v}); array];
        else
            regressors.(regressors2stack{v}) = [regressors.(regressors2stack{v}); S.regressors.(regressors2stack{v})];
        end
    end

    regressors.(partner_regressor{1}) =[ regressors.(partner_regressor{1}) ; S.regressors.play_bout_onset_regressor*partner_id];
    % Add metadata rows
    session_meta = [ ...
        repmat({animal_name}, nRows, 1), ...
        num2cell(repmat(partner_id, nRows, 1)), ...
        num2cell(repmat(ch, nRows, 1)), ...
        num2cell(play_lengths) ...
        ];
    stacked_data.meta = [stacked_data.meta; session_meta];
end
%% 3 defined periods in the baseline witoun othe rplaybouts

expanded_regressor = expand_half_intervals(regressors.play_bout_onset_regressor);
%% 4 obtain reciprocity index

bin_size =  spect_bin_size;
time = psth_structure(1).hist_range(1)+bin_size:bin_size:psth_structure(1).hist_range(2);
baseline_index = time<0;

play_bout_lengths = [stacked_data.meta{:,4}]';


recpiprocal_play_val = nan(size(play_bout_lengths,1),1);
detailed_play_self = regressors.self_onset_regressor;


for j=1:size(detailed_play_self,1)

    self_play_events_start =find( diff([0,detailed_play_self(j,:),0])==1);
    self_play_event_end =find( diff([0,detailed_play_self(j,:),0])==-1)-1;

    playbout_events_start =find( diff([0,regressors.play_bout_onset_regressor(j,:),0])==1);
    playbout_event_end =find( diff([0,regressors.play_bout_onset_regressor(j,:),0])==-1)-1;

    for this_pe = 1:numel(self_play_events_start)
        pb_wihin =    find(playbout_events_start<=self_play_events_start(this_pe) & playbout_event_end>=self_play_events_start(this_pe));
        this_pb_index = playbout_events_start(pb_wihin):playbout_event_end(pb_wihin);
        amount_of_self_play = sum(regressors.self_onset_regressor(j,this_pb_index));
        amount_of_other_play = sum(regressors.other_onset_regressor(j,this_pb_index));

        proportion_of_self_play = amount_of_self_play/(amount_of_self_play+amount_of_other_play);

        detailed_play_self(j,self_play_events_start(this_pe):self_play_event_end(this_pe)) = 1 - 2*abs(proportion_of_self_play-.5);
       
    end
    recpiprocal_play_val(j,1) = mean(detailed_play_self(j,time>=0 & time<=play_bout_lengths(j)));

end


detailed_play_other = regressors.other_onset_regressor;

for j=1:size(detailed_play_other,1)

    self_play_events_start =find( diff([0,detailed_play_other(j,:),0])==1);
    self_play_event_end =find( diff([0,detailed_play_other(j,:),0])==-1)-1;

    playbout_events_start =find( diff([0,regressors.play_bout_onset_regressor(j,:),0])==1);
    playbout_event_end =find( diff([0,regressors.play_bout_onset_regressor(j,:),0])==-1)-1;

    for this_pe = 1:numel(self_play_events_start)
        pb_wihin =    find(playbout_events_start<=self_play_events_start(this_pe) & playbout_event_end>=self_play_events_start(this_pe));
        this_pb_index = playbout_events_start(pb_wihin):playbout_event_end(pb_wihin);
        amount_of_self_play = sum(regressors.self_onset_regressor(j,this_pb_index));
        amount_of_other_play = sum(regressors.other_onset_regressor(j,this_pb_index));

        proportion_of_self_play = amount_of_other_play/(amount_of_self_play+amount_of_other_play);

        detailed_play_other(j,self_play_events_start(this_pe):self_play_event_end(this_pe)) = 1 - 2*abs(proportion_of_self_play-.5);
    end
recpiprocal_play_val(j,2) = mean(detailed_play_other(j,time>=0 & time<=play_bout_lengths(j)));

end



detailed_playbout = regressors.play_bout_onset_regressor;




for j=1:size(detailed_playbout,1)

 
    playbout_events_start =find( diff([0,regressors.play_bout_onset_regressor(j,:),0])==1);
    playbout_event_end =find( diff([0,regressors.play_bout_onset_regressor(j,:),0])==-1)-1;

    for this_pb = 1:numel(playbout_events_start)
        
        this_pb_index = playbout_events_start(this_pb):playbout_event_end(this_pb);
        amount_of_self_play = sum(regressors.self_onset_regressor(j,this_pb_index));
        amount_of_other_play = sum(regressors.other_onset_regressor(j,this_pb_index));

        proportion_of_self_play = amount_of_self_play/(amount_of_self_play+amount_of_other_play);

        detailed_playbout(j,this_pb_index) = 1 - 2*abs(proportion_of_self_play-.5);
       
    end

end

%% 6 CREAT REGRESSOR MATRIX AND PLOT SOME VARAIBLES (relevant section)

play_bout_lengths = [stacked_data.meta{:,4}]';


figure
bin_size =  spect_bin_size;
time = psth_structure(1).hist_range(1)+bin_size:bin_size:psth_structure(1).hist_range(2);
baseline_index = time<0;

play_bout_lengths       = [stacked_data.meta{:,4}]';

[sorted_play_bout_length, order] = sort(play_bout_lengths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Baseline correction and smooting %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pow_matrix_high_Freq = stacked_data.play_bout_onset_high_freq;
for j=1:size(pow_matrix_high_Freq,1)
    pow_matrix_high_Freq(j,:) = (pow_matrix_high_Freq(j,:) - mean(pow_matrix_high_Freq(j, baseline_index), 'omitmissing'))/std(pow_matrix_high_Freq(j, baseline_index), 'omitmissing');
    pow_matrix_high_Freq(j,:)  = movmean( pow_matrix_high_Freq(j,:),.125/bin_size, 'omitmissing');
end

pow_matrix_low_Freq = stacked_data.play_bout_onset_low_freq;
for j=1:size(pow_matrix_low_Freq,1)
    pow_matrix_low_Freq(j,:) = (pow_matrix_low_Freq(j,:) - mean(pow_matrix_low_Freq(j, baseline_index), 'omitmissing'))/std(pow_matrix_low_Freq(j, baseline_index), 'omitmissing');
    pow_matrix_low_Freq(j,:)  = movmean( pow_matrix_low_Freq(j,:),1/bin_size, 'omitmissing');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% frequency selection line %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array = pow_matrix_high_Freq;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array(expanded_regressor==0) = NaN;
pow_regressor = array(:);

subplot(5,6,1:6:6*3)
pcolor(time, 1:numel(sorted_play_bout_length), array(order,:))
shading flat
hold on
plot((1:numel(sorted_play_bout_length))*0, 1:numel(sorted_play_bout_length), 'w')
plot(sorted_play_bout_length, 1:numel(sorted_play_bout_length), 'w')
axis xy
xlim([-2 2])
clim([-1 2])
title('Power')
subplot(5,6,(6*3+ 1):6:6*5)
plot(time, mean(array, 'omitmissing'))
xlim([-2 2])




array = regressors.animal_speed_onset;
array(expanded_regressor==0) = NaN;
self_speed_regressor = array(:);
subplot(5,6,(1:6:6*3)  +1)
pcolor(time, 1:numel(sorted_play_bout_length), array(order,:))
shading flat
hold on
plot((1:numel(sorted_play_bout_length))*0, 1:numel(sorted_play_bout_length), 'w')
plot(sorted_play_bout_length, 1:numel(sorted_play_bout_length), 'w')
axis xy
xlim([-2 2])
title('animal speed onset')

subplot(5,6,((6*3+ 1):6:6*5) + 1)
plot(time, mean(array, 'omitmissing'))
xlim([-2 2])


array = regressors.partner_speed_onset;
array(expanded_regressor==0) = NaN;
other_speed_regressro   = array(:);
subplot(5,6,(1:6:6*3)  +2)
pcolor(time, 1:numel(sorted_play_bout_length), array(order,:))
shading flat
hold on
plot((1:numel(sorted_play_bout_length))*0, 1:numel(sorted_play_bout_length), 'w')
plot(sorted_play_bout_length, 1:numel(sorted_play_bout_length), 'w')
axis xy
xlim([-2 2])
title('partner speed onset')

subplot(5,6,((6*3+ 1):6:6*5) + 2)
plot(time, mean(array, 'omitmissing'))
xlim([-2 2])


array = regressors.animal_accel_onset;
array(expanded_regressor==0) = NaN;
slef_acc_regressro      = array(:);
subplot(5,6,(1:6:6*3)  +3)
pcolor(time, 1:numel(sorted_play_bout_length), array(order,:))
shading flat
hold on
plot((1:numel(sorted_play_bout_length))*0, 1:numel(sorted_play_bout_length), 'w')
plot(sorted_play_bout_length, 1:numel(sorted_play_bout_length), 'w')
axis xy
xlim([-2 2])
title('animal acc onset')

subplot(5,6,((6*3+ 1):6:6*5) + 3)
plot(time, mean(array, 'omitmissing'))
xlim([-2 2])




array = regressors.partner_accel_onset;
array(expanded_regressor==0) = NaN;
other_acc_regressro   = array(:);
subplot(5,6,(1:6:6*3)  +4)
pcolor(time, 1:numel(sorted_play_bout_length), array(order,:))
shading flat
hold on
plot((1:numel(sorted_play_bout_length))*0, 1:numel(sorted_play_bout_length), 'w')
plot(sorted_play_bout_length, 1:numel(sorted_play_bout_length), 'w')
axis xy
xlim([-2 2])
title('partner acc onset')

subplot(5,6,((6*3+ 1):6:6*5) + 4)
plot(time, mean(array, 'omitmissing'))
xlim([-2 2])



call_matrix = regressors.call_onset_regressor;

for j=1:size(pow_matrix_low_Freq,1)
    call_matrix(j,:)  = movmean( call_matrix(j,:),.2/bin_size, 'omitmissing');
end

array = call_matrix;
array(~expanded_regressor) = NaN;
call_regressor =   array(:);

subplot(5,6,(1:6:6*3)  +5)
pcolor(time, 1:numel(sorted_play_bout_length), array(order,:))
shading flat
hold on
plot((1:numel(sorted_play_bout_length))*0, 1:numel(sorted_play_bout_length), 'w')
plot(sorted_play_bout_length, 1:numel(sorted_play_bout_length), 'w')
axis xy
xlim([-2 2])
title('carr rate onset')

subplot(5,6,((6*3+ 1):6:6*5) + 5)
plot(time, smoothdata(mean(array, 'omitmissing'),50, 'movmedian'))
xlim([-2 2])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Now generate remaining regressros %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


array = repmat(stacked_data.meta(:,1),1, size(regressors.call_onset_regressor,2));
array(~expanded_regressor) = {'NaN'};
animal_name = categorical(array);
array = repmat(stacked_data.meta(:,3),1, size(regressors.call_onset_regressor,2));
array = cell2mat(array);
array(~expanded_regressor) = NaN;
channel_num = array;

array = regressors.relative_distance_onset;
array(expanded_regressor==0) = NaN;
animal_dist_regressor   = array(:);

array = regressors.play_bout_onset_regressor;
array(expanded_regressor==0) = NaN;
play_bout = array(:);

array = regressors.self_onset_regressor;
array(expanded_regressor==0) = NaN;
self_regressor = array(:);

array = regressors.other_onset_regressor;
array(expanded_regressor==0) = NaN;
other_regressor = array(:);

array = regressors.partner_angle_speed_onset;
array(expanded_regressor==0) = NaN;
other_angle_speed_regressor = array(:);

array = regressors.partner_angle_acc_onset;
array(expanded_regressor==0) = NaN;
other_angle_acc_regressor = array(:);

array = regressors.animal_angle_speed_onset;
array(expanded_regressor==0) = NaN;
self_angle_speed_regressor = array(:);

array = regressors.animal_angle_acc_onset;
array(expanded_regressor==0) = NaN;
self_angle_acc_regressor = array(:);

array = regressors.relative_speed_onset;
array(expanded_regressor==0) = NaN;
relative_speed_regressor = array(:);


array = regressors.PartnerNumber;
array(expanded_regressor==0) = NaN;
PartnerNumber_regressor = array(:);

array = detailed_play_other;
array(expanded_regressor==0) = NaN;
detailed_play_other_regressor = array(:);

array =detailed_play_self;
array(expanded_regressor==0) = NaN;
detailed_play_self_regressor = array(:);


array = detailed_playbout;
array(expanded_regressor==0) = NaN;
detailed_playbout_regressor = array(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Now create table %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



regressor_table = table(call_regressor, other_speed_regressro,self_speed_regressor,other_acc_regressro,slef_acc_regressro, ...
    self_regressor, other_regressor,other_angle_speed_regressor,other_angle_acc_regressor,relative_speed_regressor, ...
    self_angle_speed_regressor, self_angle_acc_regressor,PartnerNumber_regressor, ...
    detailed_play_other_regressor,detailed_play_self_regressor,detailed_playbout_regressor, ...
    play_bout,animal_name(:),channel_num(:), animal_dist_regressor,pow_regressor, ...
    'VariableNames', {'Call', 'OtherSpeed','SelfSpeed','OtherAcc','SelfAcc', ...
    'SelfBool', 'OtherBool','OtherAngleSpeed','OtherAngleAccc','RelativeSpeed', ...
    'SelfAngleSpeed','SelfAngleAccc', 'PartnerNumber',...
     'SelfPlayRecipr', 'OtherPlayRecipr','PlayboutRecipr',...
      'PlayBout','Animal','ChannelNum','AnimalDist','Power'});


predictors = {'Call', 'OtherSpeed','SelfSpeed','OtherAcc','SelfAcc', ...
    'SelfBool', 'OtherBool','OtherAngleSpeed','OtherAngleAccc','RelativeSpeed', ...
    'SelfAngleSpeed','SelfAngleAccc', 'PartnerNumber', ...
      'PlayBout','AnimalDist', 'SelfPlayRecipr', 'OtherPlayRecipr'};

regressor_table(any(isnan(regressor_table{:,predictors}),2),:) = [];

%%

%%
%% ESTIMATE R2 contribution of each variable
% Parameters
k = 30; % number of folds
cv = cvpartition(height(regressor_table), 'KFold', k);

 predictors = {'Call', 'OtherSpeed','SelfSpeed','OtherAcc','SelfAcc', ...
   'SelfBool', 'OtherBool','OtherAngleSpeed','OtherAngleAccc','RelativeSpeed', ...
    'SelfAngleSpeed','SelfAngleAccc', 'PartnerNumber', ...
      'PlayBout','AnimalDist'}; 
% Full model formula
% formula_full = 'Power ~ Call + SelfSpeed + OtherSpeed + OtherAcc + SelfAcc + PlayBout + AnimalDist + (1|Animal)';
% Join predictor names with ' + '
rhs = strjoin(predictors, ' + ');

% Build formula string
formula_full = ['Power ~ ' rhs, ' +  (1|Animal)'];

% Prepare structure to store fold errors
cvResults = struct();

% Compute full model MSEs & R² once
mse_full_allfolds = zeros(k,1);
r2_full_allfolds = zeros(k,1);
bethas_fulls_allfolds = cell(k,1);
for fold = 1:k
    trainIdx = training(cv, fold);
    testIdx = test(cv, fold);
    trainData = regressor_table(trainIdx,:);
    testData = regressor_table(testIdx,:);
    fullModel = fitlme(trainData, formula_full);
   bethas_fulls_allfolds{fold} = fullModel.Coefficients.Estimate;
    % Predictions
    yPred = predict(fullModel, testData);
    yTrue = testData.Power;

    % MSE
    mse_full_allfolds(fold) = mean((yTrue - yPred).^2,'omitmissing');

    % R² (manual)
    ss_res = nansum((yTrue - yPred).^2);
    ss_tot = nansum((yTrue - mean(yTrue,'omitnan')).^2);
    r2_full_allfolds(fold) = 1 - ss_res/ss_tot;
end
avgMSE_full = mean(mse_full_allfolds);
avgR2_full = mean(r2_full_allfolds);

fprintf('Full model average CV MSE: %.6f\n', avgMSE_full);
fprintf('Full model average CV R²: %.6f\n\n', avgR2_full);

% Loop over predictors to remove one at a time
for i = 1:numel(predictors)
    toRemove = predictors{i};
    fprintf('Testing predictor removal: %s\n', toRemove);
    
    % Regex pattern to safely remove predictor
   
    if i==1
        pattern = ['(\s*\+\s*)?' toRemove '(\s*\+\s*)?'];
        formula_reduced = regexprep(formula_full, pattern, '');
        formula_reduced = regexprep(formula_reduced, '^\s*\+\s*', '');    
        formula_reduced = regexprep(formula_reduced, '\s*\+\s*$', '');    
        formula_reduced = regexprep(formula_reduced, '\s*\+\s*\+', ' + '); 
        formula_reduced = strtrim(formula_reduced);
    else
        formula_reduced = regexprep(formula_full, ['\s*\+\s*' toRemove], '');
    end

    

    mse_reduced_allfolds = zeros(k,1);
    r2_reduced_allfolds = zeros(k,1);
    for fold = 1:13
        trainIdx = training(cv, fold);
        testIdx = test(cv, fold);
        trainData = regressor_table(trainIdx,:);
        testData = regressor_table(testIdx,:);

        reducedModel = fitlme(trainData, formula_reduced);
        yPred = predict(reducedModel, testData);
        yTrue = testData.Power;

        % MSE
        mse_reduced_allfolds(fold) = mean((yTrue - yPred).^2,'omitmissing');

        % R² (manual)
        ss_res = nansum((yTrue - yPred).^2);
        ss_tot = nansum((yTrue - mean(yTrue,'omitnan')).^2);
        r2_reduced_allfolds(fold) = 1 - ss_res/ss_tot;
    end

    avgMSE_reduced = mean(mse_reduced_allfolds);
    avgR2_reduced = mean(r2_reduced_allfolds);

    deltaMSE_folds = mse_reduced_allfolds - mse_full_allfolds;
    deltaR2_folds = r2_reduced_allfolds - r2_full_allfolds;

    fprintf('  Avg MSE increase when removing %s: %.6f\n', toRemove, avgMSE_reduced - avgMSE_full);
    fprintf('  Avg R² decrease when removing %s: %.6f\n\n', toRemove, avgR2_reduced - avgR2_full);

    % Store results for visualization
    cvResults.(toRemove).mse_full = mse_full_allfolds;
    cvResults.(toRemove).mse_reduced = mse_reduced_allfolds;
    cvResults.(toRemove).delta_mse = deltaMSE_folds;
    cvResults.(toRemove).avg_delta_mse = avgMSE_reduced - avgMSE_full;

    cvResults.(toRemove).r2_full = r2_full_allfolds;
    cvResults.(toRemove).r2_reduced = r2_reduced_allfolds;
    cvResults.(toRemove).delta_r2 = deltaR2_folds;
    cvResults.(toRemove).avg_delta_r2 = avgR2_reduced - avgR2_full;
end
cvResults.k = k;
cvResults.mse_full_allfolds = mse_full_allfolds;
cvResults.r2_full_allfolds = r2_full_allfolds;
cvResults.predictors = predictors;

cvResults.bethas_fulls_allfolds = bethas_fulls_allfolds;

play_song([],[],[])


%%  save data

% save([saving_folder,'\cvResults_mean_calls_gamma_AlllVar_play_baut.mat'],'cvResults')



