%% 1.1 LOAD PREDICTOR LIST
% Regressor PSTH structs live under Data\Analysis results\Theta psth.
saving_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\Theta psth';
results_folder = fullfile(fileparts(saving_folder), 'psth power by frequency and behavior');

% Use the saved delta model as the canonical predictor set and ordering for
% follow-up theta/gamma models.
delta_reference_file = fullfile(results_folder, 'cvResults_mean_calls_delta_v2_AlllVar.mat');
delta_reference = load(delta_reference_file, 'cvResults');
reference_predictors = delta_reference.cvResults.predictors(:)';
if isstring(reference_predictors)
    reference_predictors = cellstr(reference_predictors);
elseif ~iscellstr(reference_predictors)
    error('Delta reference does not contain a valid cvResults.predictors list: %s', ...
        delta_reference_file);
end
if isempty(reference_predictors)
    error('Delta reference predictor list is empty: %s', delta_reference_file);
end
clear delta_reference

predictors = reference_predictors;
disp('Standard variable list')
disp(predictors)

%% 1.2 SET ANALYSIS OPTIONS
% Keep the generic PlayBout regressor in the main model. SelfBool and
% OtherBool are tested below as (1) an orthogonal actor contrast and
% (2) separate overlap-aware contributions on balanced actor states.
variables2remove = {'SelfBool','OtherBool'};
variables2add = {};
predictors(ismember(predictors, variables2remove)) = [];
predictors = [predictors, variables2add];
predictors = unique(predictors, 'stable');

if ~ismember('PlayBout', predictors)
    error('The actor analyses require PlayBout in the main predictor list.');
end

disp('Updated variable list')
disp(predictors)

random_seed = 1;
rng(random_seed, 'twister');
k_requested = 30;
n_balance_repeats = 10;
low_freq_bool = true;
frequency_band_names = {'delta','theta'}; % Active input contains delta/theta.

%% 1.3 LOAD STRUCTURE DATA
disp('Loading')
% load([saving_folder,'\psth_structure_all_regressors_beta_gamma.mat'],'psth_structure'); 
% load([saving_folder,'\animal_names_all_regressors_beta_gamma.mat'],'animal_names');

load([saving_folder,'\psth_structure_all_regressors_delta_theta.mat'],'psth_structure'); 
load([saving_folder,'\animal_names_all_regressors_delta.mat'],'animal_names');
disp('Loading done')
%% 2 stacking data


% for DELTA AND THETA
freq_range = psth_structure(1).freq_range;
wind_length = psth_structure(1).wind_length;
spect_bin_size = psth_structure(1).spect_bin_size;
if low_freq_bool
    selected_frequency_index = 1;
    power_source = 'play_bout_onset_low_freq';
else
    selected_frequency_index = 2;
    power_source = 'play_bout_onset_high_freq';
end
analysis_band = frequency_band_names{selected_frequency_index};
selected_frequency_range = freq_range{selected_frequency_index};
output_file = fullfile(saving_folder, ...
    sprintf('cvResults_%s_nur_play_bout_redo_261101.mat', analysis_band));

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
stacked_data.bout_group = {};
stacked_data.bout_start = [];
stacked_data.bout_end = [];
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

    % Partner identity applies to the entire peri-bout window. Samples
    % outside the selected analysis interval are masked below by
    % expanded_regressor, rather than encoding non-play samples as partner 0.
    partner_number = repmat(partner_id, size(S.regressors.play_bout_onset_regressor));
    regressors.(partner_regressor{1}) = [regressors.(partner_regressor{1}); partner_number];
    % Add metadata rows
    session_meta = [ ...
        repmat({animal_name}, nRows, 1), ...
        num2cell(repmat(partner_id, nRows, 1)), ...
        num2cell(repmat(ch, nRows, 1)), ...
        num2cell(play_lengths) ...
        ];
    stacked_data.meta = [stacked_data.meta; session_meta];
    % The same physical bout can appear once per channel. Excluding the
    % channel from this key keeps all copies in the same CV fold.
    bout_group = repmat(string(animal_name), nRows, 1) + "|" + ...
        repmat(string(partner_id), nRows, 1) + "|" + ...
        compose('%.9f', S.play_bouts_table(:,1)) + "|" + ...
        compose('%.9f', S.play_bouts_table(:,2));
    stacked_data.bout_group = [stacked_data.bout_group; cellstr(bout_group)];
    stacked_data.bout_start = [stacked_data.bout_start; S.play_bouts_table(:,1)];
    stacked_data.bout_end = [stacked_data.bout_end; S.play_bouts_table(:,2)];
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
end

pow_matrix_low_Freq = stacked_data.play_bout_onset_low_freq;
for j=1:size(pow_matrix_low_Freq,1)
    pow_matrix_low_Freq(j,:) = (pow_matrix_low_Freq(j,:) - mean(pow_matrix_low_Freq(j, baseline_index), 'omitmissing'))/std(pow_matrix_low_Freq(j, baseline_index), 'omitmissing');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% frequency selection line %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nFrequency selection: %s band (%.3g-%.3g Hz)\n', ...
    upper(analysis_band), selected_frequency_range(1), selected_frequency_range(end));
fprintf('Power source: %s; output: %s\n\n', power_source, output_file);
if low_freq_bool
    array = pow_matrix_low_Freq;
else
    array = pow_matrix_high_Freq;
end
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

% Every matrix row is one play-bout-centred observation. The physical-bout
% key groups duplicate channel copies before it is expanded over time bins.
[bout_id_per_row, ~] = findgroups(string(stacked_data.bout_group));
bout_id_matrix = repmat(bout_id_per_row, ...
    1, size(regressors.play_bout_onset_regressor,2));
bout_id = bout_id_matrix(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Now create table %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



regressor_table = table(call_regressor, other_speed_regressro,self_speed_regressor,other_acc_regressro,slef_acc_regressro, ...
    self_regressor, other_regressor,other_angle_speed_regressor,other_angle_acc_regressor,relative_speed_regressor, ...
    self_angle_speed_regressor, self_angle_acc_regressor,PartnerNumber_regressor, ...
    detailed_play_other_regressor,detailed_play_self_regressor,detailed_playbout_regressor, ...
    play_bout,animal_name(:),channel_num(:), animal_dist_regressor,bout_id,pow_regressor, ...
    'VariableNames', {'Call', 'OtherSpeed','SelfSpeed','OtherAcc','SelfAcc', ...
    'SelfBool', 'OtherBool','OtherAngleSpeed','OtherAngleAccc','RelativeSpeed', ...
    'SelfAngleSpeed','SelfAngleAccc', 'PartnerNumber',...
     'OtherPlayRecipr', 'SelfPlayRecipr','PlayboutRecipr',...
      'PlayBout','Animal','ChannelNum','AnimalDist','BoutID','Power'});


missing_predictors = predictors(~ismember(predictors, regressor_table.Properties.VariableNames));
if ~isempty(missing_predictors)
    error('Delta-reference predictors missing from regressor_table: %s', ...
        strjoin(missing_predictors, ', '));
end

required_variables = unique([predictors, {'Power','SelfBool','OtherBool','PlayBout'}], 'stable');
regressor_table(any(ismissing(regressor_table(:,required_variables)),2),:) = [];

%% 6 ESTIMATE CROSS-VALIDATED R2 CONTRIBUTIONS
% Assign folds by bout, not by individual time bin. This prevents adjacent
% and autocorrelated bins from the same bout appearing in train and test.
[fold_assignment, k] = grouped_kfold(regressor_table.BoutID, k_requested);

rhs = strjoin(predictors, ' + ');
formula_full = ['Power ~ ' rhs ' + (1|Animal)'];
formula_actor_contrast = ['Power ~ ' rhs ' + ActorContrast + (1|Animal)'];

cvResults = struct();
predictor_diagnostics = struct();
mse_full_allfolds = nan(k,1);
r2_full_allfolds = nan(k,1);
bethas_fulls_allfolds = cell(k,1);
full_model_coefficients_allfolds = cell(k,1);

actor_contrast = struct();
actor_contrast.mse = nan(k,1);
actor_contrast.r2 = nan(k,1);
actor_contrast.n_self_only_train = nan(k,1);
actor_contrast.n_other_only_train = nan(k,1);
actor_contrast.self_weight = nan(k,1);
actor_contrast.other_weight = nan(k,1);
actor_contrast.playbout_dot_train = nan(k,1);

for fold = 1:k
    trainIdx = fold_assignment ~= fold;
    testIdx = fold_assignment == fold;
    trainData = regressor_table(trainIdx,:);
    testData = regressor_table(testIdx,:);

    fullModel = fitlme(trainData, formula_full);
    bethas_fulls_allfolds{fold} = fullModel.Coefficients.Estimate;
    full_model_coefficients_allfolds{fold} = fullModel.Coefficients;
    [mse_full_allfolds(fold), r2_full_allfolds(fold)] = ...
        score_lme(fullModel, testData);

    % The contrast is estimated from exclusive self-only and other-only
    % bins. Simultaneous and neither-active play bins receive contrast 0
    % and remain in the model.
    [trainData, testData, contrast_info] = add_actor_contrast(trainData, testData);
    actor_contrast.n_self_only_train(fold) = contrast_info.n_self;
    actor_contrast.n_other_only_train(fold) = contrast_info.n_other;
    actor_contrast.self_weight(fold) = contrast_info.self_weight;
    actor_contrast.other_weight(fold) = contrast_info.other_weight;
    if contrast_info.valid
        actor_contrast.playbout_dot_train(fold) = ...
            trainData.PlayBout' * trainData.ActorContrast;
        actorModel = fitlme(trainData, formula_actor_contrast);
        [actor_contrast.mse(fold), actor_contrast.r2(fold)] = ...
            score_lme(actorModel, testData);
    end
end

avgMSE_full = mean(mse_full_allfolds, 'omitnan');
avgR2_full = mean(r2_full_allfolds, 'omitnan');
actor_contrast.delta_mse = mse_full_allfolds - actor_contrast.mse;
actor_contrast.delta_r2 = actor_contrast.r2 - r2_full_allfolds;
actor_contrast.mse_improvement = actor_contrast.delta_mse;
actor_contrast.r2_gain = actor_contrast.delta_r2;
actor_contrast.avg_delta_mse = mean(actor_contrast.delta_mse, 'omitnan');
actor_contrast.avg_delta_r2 = mean(actor_contrast.delta_r2, 'omitnan');
actor_contrast.formula = formula_actor_contrast;
actor_contrast.definition = ...
    'Self-only=nOther/(nSelf+nOther), other-only=-nSelf/(nSelf+nOther), both/neither/non-play=0; weights estimated in each training fold.';

fprintf('Full model average grouped-CV MSE: %.6f\n', avgMSE_full);
fprintf('Full model average grouped-CV R²: %.6f\n', avgR2_full);
fprintf('Actor contrast average added CV R²: %.6f\n\n', actor_contrast.avg_delta_r2);

% Drop one predictor at a time using the exact same grouped folds.
for i = 1:numel(predictors)
    toRemove = predictors{i};
    fprintf('Testing predictor removal: %s\n', toRemove);

    kept = predictors(~strcmp(predictors, toRemove));
    if numel(kept) ~= numel(predictors) - 1
        error('Failed to drop %s from the reduced model.', toRemove);
    end
    formula_reduced = ['Power ~ ' strjoin(kept, ' + ') ' + (1|Animal)'];

    mse_reduced_allfolds = nan(k,1);
    r2_reduced_allfolds = nan(k,1);
    for fold = 1:k
        trainData = regressor_table(fold_assignment ~= fold,:);
        testData = regressor_table(fold_assignment == fold,:);
        reducedModel = fitlme(trainData, formula_reduced);
        [mse_reduced_allfolds(fold), r2_reduced_allfolds(fold)] = ...
            score_lme(reducedModel, testData);
    end

    avgMSE_reduced = mean(mse_reduced_allfolds, 'omitnan');
    avgR2_reduced = mean(r2_reduced_allfolds, 'omitnan');
    deltaMSE_folds = mse_reduced_allfolds - mse_full_allfolds;
    deltaR2_folds = r2_reduced_allfolds - r2_full_allfolds;

    fprintf('  Avg MSE increase when removing %s: %.6f\n', ...
        toRemove, avgMSE_reduced - avgMSE_full);
    fprintf('  Avg R² change when removing %s: %.6f\n\n', ...
        toRemove, avgR2_reduced - avgR2_full);

    cvResults.(toRemove).mse_full = mse_full_allfolds;
    cvResults.(toRemove).mse_reduced = mse_reduced_allfolds;
    cvResults.(toRemove).delta_mse = deltaMSE_folds;
    cvResults.(toRemove).avg_delta_mse = avgMSE_reduced - avgMSE_full;
    predictor_diagnostics.(toRemove).mse_increase_when_removed = deltaMSE_folds;
    predictor_diagnostics.(toRemove).avg_mse_increase_when_removed = ...
        avgMSE_reduced - avgMSE_full;
    cvResults.(toRemove).r2_full = r2_full_allfolds;
    cvResults.(toRemove).r2_reduced = r2_reduced_allfolds;
    cvResults.(toRemove).delta_r2 = deltaR2_folds;
    cvResults.(toRemove).avg_delta_r2 = avgR2_reduced - avgR2_full;
    predictor_diagnostics.(toRemove).r2_loss_when_removed = -deltaR2_folds;
    predictor_diagnostics.(toRemove).avg_r2_loss_when_removed = ...
        avgR2_full - avgR2_reduced;
end

%% 7 BALANCED SELF-VERSUS-OTHER CONTRIBUTIONS, INCLUDING SIMULTANEOUS PLAY
% Balance the observed members of five possible states: non-play, play with
% neither actor active, self-only, other-only, and both active. ActorState
% is categorical, so simultaneous play has its own freely estimated level.
% Reduced models merge levels according to whether SelfBool or OtherBool
% was removed. At most one bin per physical bout and state is sampled, so
% long bouts cannot dominate merely because they contain more bins.
actor_state_names = {'NonPlay','PlayNeither','SelfOnly','OtherOnly','Both'};
actor_state = actor_state_code(regressor_table);
states_present = unique(actor_state)';
regressor_table.ActorState = categorical( ...
    actor_state, states_present, actor_state_names(states_present));

actor_nuisance_predictors = predictors(~strcmp(predictors, 'PlayBout'));
actor_nuisance_rhs = strjoin(actor_nuisance_predictors, ' + ');
formula_actor_full = ['Power ~ ' actor_nuisance_rhs ' + ActorState + (1|Animal)'];
formula_actor_base = formula_full;
formula_without_self = ['Power ~ ' rhs ' + OtherBool + (1|Animal)'];
formula_without_other = ['Power ~ ' rhs ' + SelfBool + (1|Animal)'];

balanced_actor = struct();
balanced_actor.r2_full = nan(k,n_balance_repeats);
balanced_actor.r2_base = nan(k,n_balance_repeats);
balanced_actor.r2_without_self = nan(k,n_balance_repeats);
balanced_actor.r2_without_other = nan(k,n_balance_repeats);
balanced_actor.mse_full = nan(k,n_balance_repeats);
balanced_actor.mse_base = nan(k,n_balance_repeats);
balanced_actor.mse_without_self = nan(k,n_balance_repeats);
balanced_actor.mse_without_other = nan(k,n_balance_repeats);
balanced_actor.n_per_state_train = nan(k,n_balance_repeats);
balanced_actor.n_per_state_test = nan(k,n_balance_repeats);
balanced_actor.n_bouts_train = nan(k,n_balance_repeats);
balanced_actor.n_bouts_test = nan(k,n_balance_repeats);
balanced_actor.train_indices = cell(k,n_balance_repeats);
balanced_actor.test_indices = cell(k,n_balance_repeats);
balanced_actor.n_repeats = n_balance_repeats;

total_balance_iterations = k * n_balance_repeats;
fprintf('\nStarting balanced actor analysis: %d folds x %d repeats = %d iterations.\n', ...
    k, n_balance_repeats, total_balance_iterations);
for fold = 1:k
    trainCandidates = find(fold_assignment ~= fold);
    testCandidates = find(fold_assignment == fold);
    for repeat = 1:n_balance_repeats
        completed_iterations = (fold - 1) * n_balance_repeats + repeat;
        fprintf('Balanced actor analysis: fold %d/%d, repeat %d/%d (%d/%d, %.1f%%)\n', ...
            fold, k, repeat, n_balance_repeats, completed_iterations, ...
            total_balance_iterations, 100 * completed_iterations / total_balance_iterations);
        [trainBalanced, nTrain] = balanced_state_sample( ...
            trainCandidates, actor_state, states_present, regressor_table.BoutID);
        [testBalanced, nTest] = balanced_state_sample( ...
            testCandidates, actor_state, states_present, regressor_table.BoutID);
        balanced_actor.n_per_state_train(fold,repeat) = nTrain;
        balanced_actor.n_per_state_test(fold,repeat) = nTest;
        balanced_actor.train_indices{fold,repeat} = trainBalanced;
        balanced_actor.test_indices{fold,repeat} = testBalanced;

        if isempty(trainBalanced) || isempty(testBalanced)
            continue
        end
        balanced_actor.n_bouts_train(fold,repeat) = ...
            numel(unique(regressor_table.BoutID(trainBalanced)));
        balanced_actor.n_bouts_test(fold,repeat) = ...
            numel(unique(regressor_table.BoutID(testBalanced)));

        trainData = regressor_table(trainBalanced,:);
        testData = regressor_table(testBalanced,:);
        actorFullModel = fitlme(trainData, formula_actor_full);
        actorBaseModel = fitlme(trainData, formula_actor_base);
        actorWithoutSelfModel = fitlme(trainData, formula_without_self);
        actorWithoutOtherModel = fitlme(trainData, formula_without_other);

        [balanced_actor.mse_full(fold,repeat), balanced_actor.r2_full(fold,repeat)] = ...
            score_lme(actorFullModel, testData);
        [balanced_actor.mse_base(fold,repeat), balanced_actor.r2_base(fold,repeat)] = ...
            score_lme(actorBaseModel, testData);
        [balanced_actor.mse_without_self(fold,repeat), ...
            balanced_actor.r2_without_self(fold,repeat)] = ...
            score_lme(actorWithoutSelfModel, testData);
        [balanced_actor.mse_without_other(fold,repeat), ...
            balanced_actor.r2_without_other(fold,repeat)] = ...
            score_lme(actorWithoutOtherModel, testData);
    end
end
fprintf('Balanced actor analysis complete: %d/%d iterations (100%%).\n\n', ...
    total_balance_iterations, total_balance_iterations);

balanced_actor.delta_r2_joint = balanced_actor.r2_full - balanced_actor.r2_base;
balanced_actor.delta_r2_self = balanced_actor.r2_full - balanced_actor.r2_without_self;
balanced_actor.delta_r2_other = balanced_actor.r2_full - balanced_actor.r2_without_other;
balanced_actor.r2_gain_joint = balanced_actor.delta_r2_joint;
balanced_actor.r2_gain_self = balanced_actor.delta_r2_self;
balanced_actor.r2_gain_other = balanced_actor.delta_r2_other;
balanced_actor.delta_mse_joint = balanced_actor.mse_base - balanced_actor.mse_full;
balanced_actor.delta_mse_self = balanced_actor.mse_without_self - balanced_actor.mse_full;
balanced_actor.delta_mse_other = balanced_actor.mse_without_other - balanced_actor.mse_full;
balanced_actor.avg_delta_r2_joint = mean(balanced_actor.delta_r2_joint(:), 'omitnan');
balanced_actor.avg_delta_r2_self = mean(balanced_actor.delta_r2_self(:), 'omitnan');
balanced_actor.avg_delta_r2_other = mean(balanced_actor.delta_r2_other(:), 'omitnan');
balanced_actor.avg_delta_mse_joint = mean(balanced_actor.delta_mse_joint(:), 'omitnan');
balanced_actor.avg_delta_mse_self = mean(balanced_actor.delta_mse_self(:), 'omitnan');
balanced_actor.avg_delta_mse_other = mean(balanced_actor.delta_mse_other(:), 'omitnan');
balanced_actor.formula_full = formula_actor_full;
balanced_actor.formula_base = formula_actor_base;
balanced_actor.formula_without_self = formula_without_self;
balanced_actor.formula_without_other = formula_without_other;
balanced_actor.state_names = actor_state_names(states_present);
balanced_actor.state_codes = states_present;
balanced_actor.state_counts = arrayfun(@(state) nnz(actor_state == state), states_present);
balanced_actor.delta_r2_self_fold_mean = ...
    mean(balanced_actor.delta_r2_self, 2, 'omitnan');
balanced_actor.delta_r2_other_fold_mean = ...
    mean(balanced_actor.delta_r2_other, 2, 'omitnan');
balanced_actor.n_bouts_train_fold_mean = ...
    mean(balanced_actor.n_bouts_train, 2, 'omitnan');
balanced_actor.n_bouts_test_fold_mean = ...
    mean(balanced_actor.n_bouts_test, 2, 'omitnan');
valid_actor_folds = isfinite(balanced_actor.delta_r2_self_fold_mean) & ...
    isfinite(balanced_actor.delta_r2_other_fold_mean);
balanced_actor.n_valid_folds = nnz(valid_actor_folds);
if nnz(valid_actor_folds) >= 5
    contribution_difference = ...
        balanced_actor.delta_r2_self_fold_mean(valid_actor_folds) - ...
        balanced_actor.delta_r2_other_fold_mean(valid_actor_folds);
    test_train_ratio = mean( ...
        balanced_actor.n_bouts_test_fold_mean(valid_actor_folds) ./ ...
        balanced_actor.n_bouts_train_fold_mean(valid_actor_folds), 'omitnan');
    [balanced_actor.p_self_vs_other, balanced_actor.t_self_vs_other, ...
        balanced_actor.mean_self_minus_other] = corrected_cv_ttest( ...
        contribution_difference, test_train_ratio);
else
    balanced_actor.p_self_vs_other = NaN;
    balanced_actor.t_self_vs_other = NaN;
    balanced_actor.mean_self_minus_other = NaN;
    warning('EstimateR2:InsufficientActorFolds', ...
        ['Only %d grouped folds contained enough examples of every actor state. ' ...
         'At least 5 are required for the corrected self-vs-other comparison.'], ...
        nnz(valid_actor_folds));
end
balanced_actor.comparison_test = ...
    ['Nadeau-Bengio corrected paired t-test on fold-wise self-minus-other ' ...
     'delta R2; train/test ratio uses physical-bout counts. Predictive CV evidence, not a causal test.'];

fprintf('Balanced added CV R², self: %.6f\n', balanced_actor.avg_delta_r2_self);
fprintf('Balanced added CV R², other: %.6f\n', balanced_actor.avg_delta_r2_other);
fprintf('Paired self-vs-other p-value across folds: %.6g\n\n', ...
    balanced_actor.p_self_vs_other);

%% 8 FIT AND SAVE THE REPORTABLE FULL-DATA MODELS
full_model = fitlme(regressor_table, formula_full);
[actor_table, ~, final_contrast_info] = add_actor_contrast( ...
    regressor_table, regressor_table);
full_model_with_actor_contrast = fitlme(actor_table, formula_actor_contrast);
full_actor_state_model = fitlme(regressor_table, formula_actor_full);

cvResults.k = k;
cvResults.mse_full_allfolds = mse_full_allfolds;
cvResults.r2_full_allfolds = r2_full_allfolds;
cvResults.predictors = predictors;
cvResults.bethas_fulls_allfolds = bethas_fulls_allfolds;

% Preserve cvResults exactly as expected by Figure 2 and legacy scripts.
% All new diagnostics and fitted model objects live in a separate variable.
extendedResults = struct();
extendedResults.schema_version = 1;
extendedResults.k_requested = k_requested;
extendedResults.random_seed = random_seed;
extendedResults.fold_assignment = fold_assignment;
extendedResults.fold_group = 'BoutID';
extendedResults.avgMSE_full = avgMSE_full;
extendedResults.avgR2_full = avgR2_full;
extendedResults.formula_full = formula_full;
extendedResults.predictor_delta_r2_definition = ...
    'cvResults delta_r2 fields equal reduced R2 minus full R2; predictor_diagnostics r2_loss_when_removed uses full minus reduced.';
extendedResults.predictor_diagnostics = predictor_diagnostics;
extendedResults.full_model_coefficients_allfolds = full_model_coefficients_allfolds;
extendedResults.actor_contrast = actor_contrast;
extendedResults.actor_contrast.final_weights = final_contrast_info;
extendedResults.balanced_actor = balanced_actor;
extendedResults.full_model = full_model;
extendedResults.full_model_coefficients = full_model.Coefficients;
extendedResults.full_model_coefficient_ci = coefCI(full_model);
extendedResults.full_model_criterion = full_model.ModelCriterion;
extendedResults.full_model_with_actor_contrast = full_model_with_actor_contrast;
extendedResults.full_model_with_actor_contrast_coefficients = ...
    full_model_with_actor_contrast.Coefficients;
extendedResults.full_actor_state_model = full_actor_state_model;
extendedResults.full_actor_state_model_coefficients = ...
    full_actor_state_model.Coefficients;
extendedResults.full_actor_state_model_coefficient_ci = ...
    coefCI(full_actor_state_model);
extendedResults.n_observations = height(regressor_table);
extendedResults.n_bouts = numel(unique(regressor_table.BoutID));
extendedResults.observation_index = table( ...
    (1:height(regressor_table))', regressor_table.BoutID, ...
    regressor_table.Animal, regressor_table.ChannelNum, ...
    regressor_table.PlayBout, regressor_table.SelfBool, ...
    regressor_table.OtherBool, regressor_table.ActorState, fold_assignment, ...
    'VariableNames', {'ObservationIndex','BoutID','Animal','ChannelNum','PlayBout', ...
    'SelfBool','OtherBool','ActorState','Fold'});
[lookup_bout_id, lookup_row] = unique(bout_id_per_row, 'stable');
extendedResults.bout_lookup = table( ...
    lookup_bout_id, string(stacked_data.bout_group(lookup_row)), ...
    stacked_data.bout_start(lookup_row), stacked_data.bout_end(lookup_row), ...
    string(stacked_data.meta(lookup_row,1)), ...
    cell2mat(stacked_data.meta(lookup_row,2)), ...
    'VariableNames', {'BoutID','PhysicalBoutKey','Start','End','Animal','Partner'});
extendedResults.analysis_band = analysis_band;
extendedResults.frequency_range_hz = selected_frequency_range;
extendedResults.low_freq_bool = low_freq_bool;
extendedResults.power_source = power_source;
extendedResults.input_psth_file = 'psth_structure_all_regressors_delta_theta.mat';
extendedResults.input_animal_file = 'animal_names_all_regressors_delta.mat';

validate_figure2_compatibility(cvResults)

play_song([],[],[])

%% 9 SAVE DATA
disp(['Saving ', output_file])
save(output_file, 'cvResults', 'extendedResults', '-v7.3')

function [fold_assignment, k] = grouped_kfold(group_ids, k_requested)
groups = unique(group_ids);
k = min(k_requested, numel(groups));
if k < 2
    error('At least two unique bouts are required for cross-validation.');
end
group_cv = cvpartition(numel(groups), 'KFold', k);
fold_assignment = zeros(size(group_ids));
for fold = 1:k
    test_groups = groups(test(group_cv, fold));
    fold_assignment(ismember(group_ids, test_groups)) = fold;
end
end

function validate_figure2_compatibility(results)
required_top_level = [results.predictors, ...
    {'k','mse_full_allfolds','r2_full_allfolds','predictors', ...
     'bethas_fulls_allfolds'}];
unexpected = setdiff(fieldnames(results)', required_top_level);
if ~isempty(unexpected)
    error('EstimateR2:IncompatibleCvResults', ...
        'Unexpected top-level cvResults fields: %s', strjoin(unexpected, ', '));
end
for predictor_index = 1:numel(results.predictors)
    predictor = results.predictors{predictor_index};
    if ~isfield(results, predictor) || ...
            ~isfield(results.(predictor), 'r2_full') || ...
            ~isfield(results.(predictor), 'r2_reduced')
        error('EstimateR2:IncompatibleCvResults', ...
            'Predictor %s lacks Figure 2 R2 fields.', predictor);
    end
    if numel(results.(predictor).r2_full) ~= results.k || ...
            numel(results.(predictor).r2_reduced) ~= results.k
        error('EstimateR2:IncompatibleCvResults', ...
            'Predictor %s does not contain one R2 value per fold.', predictor);
    end
end
end

function [mse, r2] = score_lme(model, test_data)
y_true = test_data.Power;
y_pred = predict(model, test_data);
valid = isfinite(y_true) & isfinite(y_pred);
if nnz(valid) < 2
    mse = NaN;
    r2 = NaN;
    return
end
y_true = y_true(valid);
y_pred = y_pred(valid);
mse = mean((y_true - y_pred).^2);
ss_res = sum((y_true - y_pred).^2);
ss_tot = sum((y_true - mean(y_true)).^2);
if ss_tot == 0
    r2 = NaN;
else
    r2 = 1 - ss_res/ss_tot;
end
end

function [train_data, test_data, info] = add_actor_contrast(train_data, test_data)
self_only = train_data.PlayBout == 1 & ...
    train_data.SelfBool == 1 & train_data.OtherBool == 0;
other_only = train_data.PlayBout == 1 & ...
    train_data.SelfBool == 0 & train_data.OtherBool == 1;
info.n_self = nnz(self_only);
info.n_other = nnz(other_only);
info.valid = info.n_self > 0 && info.n_other > 0;
if info.valid
    denominator = info.n_self + info.n_other;
    info.self_weight = info.n_other / denominator;
    info.other_weight = -info.n_self / denominator;
else
    info.self_weight = NaN;
    info.other_weight = NaN;
end
train_data.ActorContrast = actor_contrast_values( ...
    train_data, info.self_weight, info.other_weight);
test_data.ActorContrast = actor_contrast_values( ...
    test_data, info.self_weight, info.other_weight);
end

function contrast = actor_contrast_values(data, self_weight, other_weight)
contrast = zeros(height(data),1);
if ~isfinite(self_weight) || ~isfinite(other_weight)
    contrast(:) = NaN;
    return
end
self_only = data.PlayBout == 1 & data.SelfBool == 1 & data.OtherBool == 0;
other_only = data.PlayBout == 1 & data.SelfBool == 0 & data.OtherBool == 1;
contrast(self_only) = self_weight;
contrast(other_only) = other_weight;
end

function state = actor_state_code(data)
state = ones(height(data),1); % 1: non-play
in_play = data.PlayBout == 1;
self_active = data.SelfBool == 1;
other_active = data.OtherBool == 1;
state(in_play & ~self_active & ~other_active) = 2; % play, neither active
state(in_play & self_active & ~other_active) = 3;  % self only
state(in_play & ~self_active & other_active) = 4;  % other only
state(in_play & self_active & other_active) = 5;   % simultaneous
end

function [sampled_indices, n_per_state] = balanced_state_sample( ...
    candidate_indices, state, states_required, group_ids)
state_representatives = cell(numel(states_required),1);
for state_index = 1:numel(states_required)
    value = states_required(state_index);
    state_pool = candidate_indices(state(candidate_indices) == value);
    state_groups = unique(group_ids(state_pool));
    representatives = nan(numel(state_groups),1);
    for group_index = 1:numel(state_groups)
        group_pool = state_pool(group_ids(state_pool) == state_groups(group_index));
        representatives(group_index) = group_pool(randi(numel(group_pool)));
    end
    state_representatives{state_index} = representatives;
end
n_per_state = min(cellfun(@numel, state_representatives));
if isempty(n_per_state) || n_per_state < 2
    sampled_indices = [];
    n_per_state = 0;
    return
end
sampled_indices = [];
for state_index = 1:numel(states_required)
    pool = state_representatives{state_index};
    selected = pool(randperm(numel(pool), n_per_state));
    sampled_indices = [sampled_indices; selected(:)]; %#ok<AGROW>
end
sampled_indices = sampled_indices(randperm(numel(sampled_indices)));
end

function [p_value, t_stat, mean_difference] = corrected_cv_ttest( ...
    differences, test_train_ratio)
differences = differences(isfinite(differences));
n = numel(differences);
mean_difference = mean(differences);
if n < 2 || ~isfinite(test_train_ratio)
    p_value = NaN;
    t_stat = NaN;
    return
end
corrected_variance = (1/n + test_train_ratio) * var(differences, 0);
if corrected_variance == 0
    if mean_difference == 0
        t_stat = 0;
        p_value = 1;
        return
    end
    t_stat = sign(mean_difference) * Inf;
else
    t_stat = mean_difference / sqrt(corrected_variance);
end
p_value = 2 * tcdf(-abs(t_stat), n - 1);
end

