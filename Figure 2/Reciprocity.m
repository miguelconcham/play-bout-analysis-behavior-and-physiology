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


%% 5 create power and reiprocity arrays

play_bout_lengths = [stacked_data.meta{:,4}]';


figure
bin_size =  spect_bin_size;
time = psth_structure(1).hist_range(1)+bin_size:bin_size:psth_structure(1).hist_range(2);
baseline_index = time<0;

play_bout_lengths       = [stacked_data.meta{:,4}]';

[sorted_play_bout_length, order] = sort(play_bout_lengths);


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

array = detailed_play_other;
array(expanded_regressor==0) = NaN;
detailed_play_other_regressor = array(:);

array =detailed_play_self;
array(expanded_regressor==0) = NaN;
detailed_play_self_regressor = array(:);


array = detailed_playbout;
array(expanded_regressor==0) = NaN;
detailed_playbout_regressor = array(:);

array = regressors.play_bout_onset_regressor;
array(expanded_regressor==0) = NaN;
play_bout = array(:);


array = repmat(stacked_data.meta(:,1),1, size(regressors.call_onset_regressor,2));
array(~expanded_regressor) = {'NaN'};
animal_name = categorical(array);

%% 6 (not shwon) plot playboput colored according to reciprocity



[sorted_play_bout_length, order] = sort(play_bout_lengths);

 


figure
colormap(1-gray)
array = detailed_play_self;
array(expanded_regressor==0) = NaN;
subplot(5,2,[1 3 5] )
pcolor(time, 1:numel(sorted_play_bout_length), array(order,:))
shading flat
hold on
plot((1:numel(sorted_play_bout_length))*0, 1:numel(sorted_play_bout_length), 'w')
plot(sorted_play_bout_length, 1:numel(sorted_play_bout_length), 'w')
axis xy
xlim([-2 5])
subplot(5,2,[7 9])
plot(time, mean(array, 'omitmissing'))
xlim([-2 5])


colormap(1-gray)
array = detailed_play_other;
array(expanded_regressor==0) = NaN;
subplot(5,2,[1 3 5] +1)
pcolor(time, 1:numel(sorted_play_bout_length), array(order,:))
shading flat
hold on
plot((1:numel(sorted_play_bout_length))*0, 1:numel(sorted_play_bout_length), 'w')
plot(sorted_play_bout_length, 1:numel(sorted_play_bout_length), 'w')
axis xy
xlim([-2 5])
subplot(5,2,[7 9]+1)
plot(time, mean(array, 'omitmissing'))
xlim([-2 5])

figure
subplot(1,2,1)
y = play_bout_lengths(:);
x = recpiprocal_play_val(:,1);

% Sort for plotting
[xs, idx] = sort(x);
ys = y(idx);

% Smooth estimate with LOESS kernel
fhat = smooth(xs, ys, 0.5, 'loess');   % 0.2 = smoothing span (adjust)

% Plot
scatter(x, y, 20, 'filled'); hold on;
plot(xs, fhat, 'LineWidth', 2);
xlabel('x'); ylabel('y'); title('Kernel-smoothed estimate of y = f(x)');


X = table(x, x.^2, y, 'VariableNames', {'x','x2','y'});

mdl = fitlm(X,'y ~ x + x2');

subplot(1,2,2)
x = recpiprocal_play_val(:,2);

% Sort for plotting
[xs, idx] = sort(x);
ys = y(idx);

% Smooth estimate with LOESS kernel
fhat = smooth(xs, ys, 0.5, 'loess');   % 0.2 = smoothing span (adjust)

% Plot
scatter(x, y, 20, 'filled'); hold on;
plot(xs, fhat, 'LineWidth', 2);
xlabel('x'); ylabel('y'); title('Kernel-smoothed estimate of y = f(x)');


%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
% for this code  you need to run first "estimate psth with all regressors.mat"
% from section 1 to 6.
% DO NOT RUN SESSION BEFORE< SINCE YOU MAY OVERWRITE DATA
% later section also should be avoided.

%% 7 select needed arrays (regressors)
array = detailed_playbout;
array(expanded_regressor==0) = NaN;
detailed_playbout_regressor_matrix = array;

array = regressors.play_bout_onset_regressor;
array(expanded_regressor==0) = NaN;
play_bout_matrix  = array;


array = detailed_play_other;
array(expanded_regressor==0) = NaN;
detailed_play_other_matrix  = array;


array = detailed_play_self;
array(expanded_regressor==0) = NaN;
detailed_play_self_matrix  = array;



array = pow_matrix_low_Freq;
array(expanded_regressor==0) = NaN;
pow_regressor_matrix  = array;


call_matrix = regressors.call_onset_regressor;
for j=1:size(pow_matrix_low_Freq,1)
    call_matrix(j,:)  = smooth( call_matrix(j,:),.5/bin_size);
end
array = call_matrix;
array(~expanded_regressor) = NaN;
call_regressor_matrix =   array;

matrixforcorr = nan(size(detailed_playbout_regressor_matrix,1),5);

for j=1:size(detailed_playbout_regressor_matrix,1)

    playbout_indexes    = play_bout_matrix(j,:)==1;
    matrixforcorr(j,1)  = mean(detailed_playbout_regressor_matrix(j,playbout_indexes));
     matrixforcorr(j,2) = mean(pow_regressor_matrix(j,playbout_indexes));
     matrixforcorr(j,3) = mean(call_regressor_matrix(j,playbout_indexes));
     matrixforcorr(j,4) = max(pow_regressor_matrix(j,playbout_indexes));
     matrixforcorr(j,5) = median(pow_regressor_matrix(j,playbout_indexes));
end


array = repmat(stacked_data.meta(:,1),1, size(regressors.call_onset_regressor,2));
array(~expanded_regressor) = {'NaN'};
animal_name = categorical(array);

animal_list = unique(animal_name);
animal_list(ismember(animal_list, 'NaN')) = [];

%% 8 (Supp Fig 3a) % of reciprocical and nonreciprocal play
index = play_bout_matrix==1 & detailed_playbout_regressor_matrix>=0;
x = detailed_playbout_regressor_matrix(index);


figure
subplot(1,2,2)
histogram(x(x>0), [0,0.01,.025:0.025:1],'FaceColor','k','Normalization','pdf', 'FaceAlpha',.5, 'EdgeColor','none')

[f,xf] = kde(x(x>0));
hold on
plot(xf, f, 'k')
subplot(1,2,1)
pie([sum(x==0) sum(x>0)], {['Nonreciprocal ', num2str(100*mean(x==0))],['Reciprocal ', num2str(100*mean(x>0))]})
%% 9 select here if you want to estimae power duirng play botius (firt group of code) or self play  (second group of code lines)
% 
% index = play_bout_matrix==1 & detailed_playbout_regressor_matrix>=0 & abs(pow_regressor_matrix)<3;
% x = detailed_playbout_regressor_matrix(index);
% y = pow_regressor_matrix(index);
% 
% animal_cat = animal_name(index);

index = play_bout_matrix==1 & detailed_play_self_matrix>=0 & abs(pow_regressor_matrix)<3;
x = detailed_play_self_matrix(index);
y = pow_regressor_matrix(index);
call_rate = call_regressor_matrix(index);

animal_cat = animal_name(index);




%% 10 (not shown) percentage of reciprocity ony when implanted animal is playing 

figure
subplot(1,2,2)
histogram(x(x>0), [0,0.01,.025:0.025:1],'FaceColor','k','Normalization','pdf', 'FaceAlpha',.5, 'EdgeColor','none')

[f,xf] = kde(x(x>0));
hold on
plot(xf, f, 'k')
subplot(1,2,1)
pie([sum(x==0) sum(x>0)])

%% 11 (Fig 2h, Supp Fig 3b,c)

figure

subplot(1,5,1)
histogram(y(x <0.1),-4:0.05:4,'Normalization','pdf', 'FaceAlpha',.5, 'EdgeColor','none')
hold on
histogram(y(x >0.7),-4:0.05:4,'Normalization','pdf', 'FaceAlpha',.5, 'EdgeColor','none')


 p = ranksum(y(x <0.1),y(x>0.7 ));
 title(p)
[f_low,xf_low] = kde(y(x <0.1));

[f_high,xf_high] = kde(y(x>0.7));

plot(xf_low, f_low, 'b')
hold on
plot(xf_high, f_high, 'r')

subplot(1,5,2)
[f_low,xf_low] = ecdf(y(x <0.1));

[f_high,xf_high] = ecdf(y(x>0.7));

plot(xf_low, f_low, 'b')
hold on
plot(xf_high, f_high, 'r')


subplot(1,5,3)
hold on
all_differences = nan(numel(animal_list),3);
for an=1:numel(animal_list)

    if any(x <0.1 & animal_cat==animal_list(an))
        p = ranksum(y(x <0.1 & animal_cat==animal_list(an)),y(x>0.7 & animal_cat==animal_list(an)));
        all_differences(an,1) =p;
        all_differences(an,2:3) = [median(y(x <0.1 & animal_cat==animal_list(an))) median(y(x>0.7 & animal_cat==animal_list(an)))];
    end
end

plot([1 2],all_differences(:,2:3), ':k')
 plot([1 2],all_differences(:,2:3), '.k', 'MarkerSize',5)
 xlim([.5 2.5])
 xticks([1 2])
 xticklabels({'Non-reciprocal','Reciprocal'})

[h,p, ~, tstats]=ttest(all_differences(:,2),all_differences(:,3))
title(num2str([p,tstats.tstat, max(all_differences(:,1))]))





%% 12 (Supp Fig 3e,f) now estiamte correlation for recirpocal index (above threshold = .75)

index = play_bout_matrix==1 & detailed_playbout_regressor_matrix>=.7 & abs(pow_regressor_matrix)<3;
x = detailed_playbout_regressor_matrix(index);
y = pow_regressor_matrix(index);
call_rate = call_regressor_matrix(index);
animal_cat = animal_name(index);



x_round = round(x*8,1)/8;
x_round = x;
x_list = unique(x_round);

mena_pow  = nan(size(x_list));
mean_pow_per_animal = nan(numel(animal_list),size(x_list,1));
predicted_mean = nan(numel(animal_list),size(x_list,1));
ci_pow  = nan(size(x_list,1),2);
number_of_data_points = nan(size(x_list));
corr_per_animal = nan(numel(animal_list),2);
animal_lm = cell(numel(animal_list),1);
    for an=1:numel(animal_list)
    animal_lm{an} = fitlm(x_round(animal_cat==animal_list(an)), y(animal_cat==animal_list(an)));

    [c,p]=corr(x_round(animal_cat==animal_list(an)), y(animal_cat==animal_list(an)), 'Type','Spearman');
         corr_per_animal(an,:) = [c,p];
         predicted_mean(an,:) = predict(animal_lm{an},x_list);
    end

for j=1:numel(x_list)
    mena_pow(j) = mean(y(x_round ==x_list(j)));
   [h,p, ci] = ttest(y(x_round ==x_list(j)));
   
   ci_pow(j,:) = ci;
    number_of_data_points(j) = sum(x_round ==x_list(j));

    for an=1:numel(animal_list)

         mean_pow_per_animal(an,j) = mean(y(x_round ==x_list(j) &  animal_cat==animal_list(an)));
    end


end
ci_pow = ci_pow';
ci_pow(1,:)  = movmean(ci_pow(1,:) ,1);
ci_pow(2,:)  = movmean(ci_pow(2,:) ,1);



subplot(1,5,4)

hold on
plot(x_list,mena_pow, '.k')
lin_model_table = array2table([x_list,mena_pow], 'VariableNames',{'Reciprocity','Power'});
linear_model = fitlm(lin_model_table)
plot(linspace(.7,1,100), predict(linear_model,linspace(.7,1,100)'), 'r')
title(num2str([linear_model.coefTest linear_model.Rsquared.Ordinary]))




subplot(1,5,5)
x_rand_loc= .2*(rand(size(corr_per_animal(:,1)))-.5) + 1;
plot(x_rand_loc,corr_per_animal(:,1), '.k')
hold on
plot(x_rand_loc(corr_per_animal(:,2)<0.01),corr_per_animal(corr_per_animal(:,2)<0.01,1), '.r')
bar(mean(corr_per_animal(:,1)), 'FaceAlpha',.1)

%% 13 (not shwon) Distribution of reciprocal play index for each animal play behavior, isntead that of playbouts
x_lim = [.75 1];
figure
condition2incluide_1 = ~isnan(pow_regressor) & ~isnan(detailed_play_self_regressor) & ~isnan(detailed_play_other_regressor) & play_bout==1;
condition2incluide_2 = ~isnan(pow_regressor) & ~isnan(detailed_play_other_regressor) & ~isnan(detailed_play_self_regressor)   & play_bout==1;

reciprocal_index_for_each_animal_play = [detailed_play_self_regressor(condition2incluide_1),detailed_play_other_regressor(condition2incluide_2)];
histogram(reciprocal_index, 'Normalization','percentage')
n_val2plot = 100000;
zscore_limit = Inf;


