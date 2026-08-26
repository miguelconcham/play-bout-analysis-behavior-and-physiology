


npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\Cross_correlogram';
% Per-animal coincidence structs (~3.5 GB) stay on DataSets; ALL-ANIMALS copy is in repo Data.
coincidence_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\Cross_correlogram';
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
figure_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure coincident phases';

animal_file_names =  cellfun(@(x) ['B', x],strsplit([animal_list.name], 'B'), 'UniformOutput',false)';
animal_file_names(1) = [];
% animal2exclude = {'B4D4 0826 Dual'};
animal2exclude = {''};
animal_list(ismember(animal_file_names,animal2exclude)) = [];
animal_names ={};
n_strctut = 1;

freq_range_1    = [1 5];
freq_range_2    = [6 12];
sr              = 2500;
filter_order    = 2000;



% Parameters for delta
Hd_freq = designfilt('bandpassfir', ...
'FilterOrder', filter_order, ...
'CutoffFrequency1', freq_range_1(1), ...
'CutoffFrequency2', freq_range_1(2), ...
'SampleRate', sr, ...
'DesignMethod', 'window', ...
'Window', 'hamming');
bin_size_freq = 0.01;

% Parameters for theta
% Hd_freq = designfilt('bandpassfir', ...
% 'FilterOrder', filter_order, ...
% 'CutoffFrequency1', freq_range_2(1), ...
% 'CutoffFrequency2', freq_range_2(2), ...
% 'SampleRate', sr, ...
% 'DesignMethod', 'window', ...
% 'Window', 'hamming');
% bin_size_freq = 0.001;
areas2analyse = {'DLPAG'	'DR'	'LPAG'	'SupCol'	'VLPAG'};

time_precision = 0.005;

bin_size = 0.002;
hist_range = [-1.5 1];

psth_edges = hist_range(1):bin_size:hist_range(2);
%%
%% load
fn =12;
disp('LOADING')
load([coincidence_folder,'\',animal_list(fn).name,'_FreqRange_',num2str([Hd_freq.CutoffFrequency1 Hd_freq.CutoffFrequency2])...
            ,'_phase_coincidence_structure.mat'],'phase_struct');
load([coincidence_folder,'\',animal_list(fn).name,'_FreqRange_',num2str([Hd_freq.CutoffFrequency1 Hd_freq.CutoffFrequency2])...
            ,'_phase_coincidence_structure_animal_names.mat'],'animal_names');
% phase_struct(6 )=[];
% animal_names(6,:)=[];
disp('LOAD READY')
%% merge data
disp('MERGING SESSIONS')

synch_spikes_phases_before  = [];
synch_spikes_phases_after   = [];
mvl_before_after = [];


all_area_comb = [];
session_id = [];
clusters_id = [];
cluster_fr = [];

for fn = 1:numel(phase_struct)
    for j=1:size(phase_struct(fn).sync_event_data,1)
        row_before = nan(1,4);
        row_after  = nan(1,4);       
        if ~isempty(phase_struct(fn).sync_event_data{j,1})
            before1 = phase_struct(fn).sync_event_data{j,1}(:,4)==1;
            before2 = phase_struct(fn).sync_event_data{j,2}(:,4)==1;
            after1 = phase_struct(fn).sync_event_data{j,1}(:,4)==2;
            after2 = phase_struct(fn).sync_event_data{j,2}(:,4)==2;
            if any(before1) && any(before2)
                row_before = [circ_mean(phase_struct(fn).sync_event_data{j,1}(before1,[2 3]))  circ_mean(phase_struct(fn).sync_event_data{j,2}(before2,[2 3]))];
            elseif any(before1)
                row_before = [circ_mean(phase_struct(fn).sync_event_data{j,1}(before1,[2 3])) [NaN NaN]];
            elseif any(before2)
                row_before = [[NaN NaN] circ_mean(phase_struct(fn).sync_event_data{j,2}(before2,[2 3]))];
            else
                row_before = [NaN NaN NaN NaN];
            end

            if any(row_before==0)
                disp('0 found')

            end
             if any(after1) && any(after2)
                row_after = [circ_mean(phase_struct(fn).sync_event_data{j,1}(after1,[2 3])) circ_mean(phase_struct(fn).sync_event_data{j,2}(after2,[2 3]))];
            elseif any(after1)
                row_after = [circ_mean(phase_struct(fn).sync_event_data{j,1}(after1,[2 3])) [NaN NaN]];
             elseif any(after2)
                 row_after = [[NaN NaN] circ_mean(phase_struct(fn).sync_event_data{j,2}(after2,[2 3]))];
             else
                 row_after = [NaN NaN NaN NaN];
             end

             if sum(before1)>1
                 % mvl_before_1 = circ_r(phase_struct(fn).sync_event_data{j,1}(before1,2));
                 
                   mvl_before_1 = ppc(phase_struct(fn).sync_event_data{j,1}(before1,2));
             else
                 mvl_before_1 = NaN;
             end

             if sum(after1)>1
                 % mvl_after_1 = circ_r(phase_struct(fn).sync_event_data{rj,1}(after1,2));
                  mvl_after_1 = ppc(phase_struct(fn).sync_event_data{j,1}(after1,2));
                 
             else
                 mvl_after_1 = NaN;
             end

             mvl_before_after = [mvl_before_after;[mvl_before_1 mvl_after_1]];


        else
                 mvl_before_after = [mvl_before_after;[NaN NaN]];

        end
        synch_spikes_phases_before          = [synch_spikes_phases_before;row_before];
        synch_spikes_phases_after           = [synch_spikes_phases_after;row_after];
    end

    all_area_comb                   = [all_area_comb;phase_struct(fn).cluster_info.area(phase_struct(fn).synch_comb)];
    session_id                      = [session_id;repmat(animal_names(fn,1),size(phase_struct(fn).synch_comb,1),1)];
    clusters_id                     = [clusters_id;phase_struct(fn).cluster_info.cluster_id(phase_struct(fn).synch_comb)];
    cluster_fr                      = [cluster_fr;phase_struct(fn).cluster_info.fr(phase_struct(fn).synch_comb)];

end
disp('DONE')
%% Obtain area combinations

area_combinations = cell2table(all_area_comb);
area_combinations.Properties.VariableNames = {'Neuron1','Neuron2'};
unique_area_combinations = unique(area_combinations, 'rows');


relevant_areas = {'SupCol','DLPAG','LPAG','VLPAG','DR'};
unique_area_combinations = unique_area_combinations(ismember(unique_area_combinations.Neuron1,relevant_areas) & ismember(unique_area_combinations.Neuron2,relevant_areas),:);


data = string(table2array(unique_area_combinations));
sortedRows = sort(data, 2);

% 2. Get a unique ID for every unique combination
[~, ~, ic] = unique(sortedRows, 'rows');

% 3. Group the original row indexes based on those unique IDs
% 'rowGroups' will be a cell array where each cell contains the indexes for a combination
rowGroups = accumarray(ic, (1:size(data,1))', [], @(x) {sort(x')});

% 4. Filter to keep only the groups with more than 1 index (the mirrored pairs)
mirroredPairs = rowGroups(cellfun(@(x) length(x) > 1, rowGroups));

for j=1:size(mirroredPairs,1)

    unique_area_combinations{mirroredPairs{j},:} = repmat(unique_area_combinations{mirroredPairs{j}(1),:},2,1);
end

unique_area_combinations = unique(unique_area_combinations, 'rows');
%% V load delta and theta phase data
data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
saving_folder = [data_root, '\Analysis results\phase locking data'];

load([saving_folder,'\theta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) =     {'isRt'  };
all_neurons_TD = all_neurons;
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner1'}))         = {'ThetaPartner1'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner2'}))         = {'ThetaPartner2'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Play'}))             = {'ThetaPlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'PrePlay'}))          = {'ThetaPrePlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'EntireSession'}))    = {'ThetaEntireSession'};

load([saving_folder,'\delta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'}))  =     {'isRt'  };
all_neurons_TD.DeltaPartner1                            = all_neurons.Partner1;
all_neurons_TD.DeltaPartner2                            = all_neurons.Partner2;
all_neurons_TD.DeltaEntireSession                       = all_neurons.EntireSession;
all_neurons_TD.DeltaPlay                                = all_neurons.Play;
all_neurons_TD.DeltaPrePlay                             = all_neurons.PrePlay;
all_neurons_TD.Exited                                   = nan(size(all_neurons_TD,1),1);
all_neurons_TD.Inhibited                                = nan(size(all_neurons_TD,1),1);

nRows = size(clusters_id,1);
idx_pairs = zeros(nRows,2);

for i = 1:nRows
    s = session_id{i};
    id_pair = clusters_id(i,:);
    for j = 1:2
        if any(strcmp(all_neurons_TD.session, s) & ...
                all_neurons_TD.cluster_id == id_pair(j), ...
                1)
            idx_pairs(i,j) = find( ...
                strcmp(all_neurons_TD.session, s) & ...
                all_neurons_TD.cluster_id == id_pair(j), ...
                1);  % only one match expected
        else
            disp('Missing neuron')
            idx_pairs(i,j)= NaN;
        end
    end
end

%%
alpha_level = 0.05;

trough          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak            = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained   = all_neurons_TD.DeltaEntireSession.PPCPval>.1;

min_fr  =0;
entrained_group1 = trough;
entrained_group2 = trough;
comparison_group1 = peak;
comparison_group2 = peak;

x_tick_labels               = cell(size(unique_area_combinations,1),1);
r_per_condition_before      = nan(size(unique_area_combinations,1),3);
angle_per_condition_before  = nan(size(unique_area_combinations,1),3);

r_per_condition_after       = nan(size(unique_area_combinations,1),3);
angle_per_condition_after   = nan(size(unique_area_combinations,1),3);
counts_per_condition        = zeros(size(unique_area_combinations,1),3);

p_val_before_after_condition         = nan(size(unique_area_combinations,1),3);

all_angles = [];



for j=1:size(unique_area_combinations,1)

   

   area_combination_indexes      = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
       (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
   trough_trough        = (trough(idx_pairs(:,1)) &  trough(idx_pairs(:,2)) & area_combination_indexes & mean(cluster_fr,2)>min_fr);
   peak_peak          = (peak(idx_pairs(:,1)) &  peak(idx_pairs(:,2)) & area_combination_indexes  & mean(cluster_fr,2)>min_fr);

   non_entrained = ( ( non_entrained(idx_pairs(:,1)) &  non_entrained(idx_pairs(:,2)) | (non_entrained(idx_pairs(:,1)) &  non_entrained(idx_pairs(:,2))) ) & area_combination_indexes  & mean(cluster_fr,2)>min_fr);
    
   if ~isempty(trough_trough)
       counts_per_condition(j,1)       = sum(trough_trough);
       alpha_before =synch_spikes_phases_before(trough_trough & ~isnan(synch_spikes_phases_before(:,1))  & ~isnan(synch_spikes_phases_after(:,1)),1);
       alpha_after =synch_spikes_phases_after(trough_trough & ~isnan(synch_spikes_phases_after(:,1)) & ~isnan(synch_spikes_phases_before(:,1)),1);

       all_angles = [all_angles;[[alpha_before;alpha_after],[alpha_before*0;alpha_after*0+1],[alpha_before;alpha_after]*0+1,[alpha_before;alpha_after]*0+j]];
       r_per_condition_before(j,1)     = circ_r(alpha_before);
       angle_per_condition_before(j,1) = circ_mean(alpha_before);
       r_per_condition_after(j,1)      = circ_r(alpha_after);
       angle_per_condition_after(j,1)  = circ_mean(alpha_after);
       % Calculate circular differences
       diff = circ_dist(alpha_before, alpha_after);
       % Test if the mean difference is 0
       if ~isempty(diff)
           [pval, z] = circ_mtest(diff, 0, alpha_level);
           p_val_before_after_condition(j,1) = pval;
       end
   end

   if ~isempty(peak_peak)
       alpha_before = synch_spikes_phases_before(peak_peak & ~isnan(synch_spikes_phases_before(:,1)) & ~isnan(synch_spikes_phases_after(:,1)),1);
       alpha_after = synch_spikes_phases_after(peak_peak & ~isnan(synch_spikes_phases_after(:,1)) & ~isnan(synch_spikes_phases_before(:,1)),1);
        all_angles = [all_angles;[[alpha_before;alpha_after],[alpha_before*0;alpha_after*0+1],[alpha_before;alpha_after]*0+2,[alpha_before;alpha_after]*0+j]];
       counts_per_condition(j,2)       = sum(peak_peak);
       r_per_condition_before(j,2)     = circ_r(alpha_before);
       angle_per_condition_before(j,2) = circ_mean(alpha_before);
       r_per_condition_after(j,2)      = circ_r(alpha_after);
       angle_per_condition_after(j,2)  = circ_mean(alpha_after);

       % Calculate circular differences
       diff = circ_dist(alpha_before, alpha_after);
       if ~isempty(diff)
           % Test if the mean difference is 0
           [pval, z] = circ_mtest(diff, 0, alpha_level);
           p_val_before_after_condition(j,2) = pval;
       end


   end

   if ~isempty(non_entrained)
       alpha_before = synch_spikes_phases_before(non_entrained & ~isnan(synch_spikes_phases_before(:,1)) & ~isnan(synch_spikes_phases_after(:,1)),1);
       alpha_after = synch_spikes_phases_after(non_entrained & ~isnan(synch_spikes_phases_after(:,1)) & ~isnan(synch_spikes_phases_before(:,1)),1);
       all_angles = [all_angles;[[alpha_before;alpha_after],[alpha_before*0;alpha_after*0+1],[alpha_before;alpha_after]*0+3,[alpha_before;alpha_after]*0+j]];
       counts_per_condition(j,3)       = sum(non_entrained );
       r_per_condition_before(j,3)     = circ_r(alpha_before);
       angle_per_condition_before(j,3) = circ_mean(alpha_before);
       r_per_condition_after(j,3)      = circ_r(alpha_after);
       angle_per_condition_after(j,3)  = circ_mean(alpha_after);

       % Calculate circular differences
       diff = circ_dist(alpha_before, alpha_after);
       if ~isempty(diff)
           % Test if the mean difference is 0
           [pval, z] = circ_mtest(diff, 0, alpha_level);
           p_val_before_after_condition(j,3) = pval;
       end
   end

end

%%  test if angle changes
min_count = 25;
var2compare = 4;
% for first column 0 is beofre and 1 is after
selct_groups_by_count = find(counts_per_condition(:,1)>min_count);
slect_specific_groups = [10 11 12 13 14 15]
index = ismember(all_angles(:,2),1) & ismember(all_angles(:,3),1) & ismember(all_angles(:,4),slect_specific_groups) & ~isnan(all_angles(:,1)) ;
% [p, table] = circ_wwtest(all_angles(index,1), all_angles(index,4));
groups = all_angles(index,var2compare);
angles = all_angles(index,1);
[p, med] = circ_cmtest(angles, groups);
%% Supplementary plot 1

circular_violin_plot(angles, groups, true);
current_labels = xticklabels;

new_labels = current_labels;

for j=1:numel(current_labels)
    j_str = str2double(current_labels{j});
    new_labels{j} = [unique_area_combinations.Neuron1{j_str}, ' ',unique_area_combinations.Neuron2{j_str}];
end
xticklabels(new_labels)

%%

print(gcf,'-vector','-dsvg',[figure_folder,'/violing plots of coincident phases only significnat combs.svg'])

%%

% [pval_area, table] = circ_hktest(all_angles(index,1), all_angles(index,4), all_angles(index,3), 1);

% 1. Setup Groups
group_ids        = unique(groups);
num_groups       = length(group_ids);
p_matrix_cm      = nan(num_groups); % Pre-allocate p-value matrix
p_matrix_ww      = nan(num_groups);
num_comparisons  = (num_groups * (num_groups - 1)) / 2;
alpha_adj        = 0.05 / num_comparisons; % Bonferroni correction

% 2. Nested Loop for Pairwise Comparison
for i = 1:num_groups
    for j = i+1:num_groups
        % Extract data for current pair
        idx1 = (groups == group_ids(i));
        idx2 = (groups == group_ids(j));
        
        subset_angles = [angles(idx1); angles(idx2)];
        subset_groups = [ones(sum(idx1),1); 2*ones(sum(idx2),1)];
        
        % Run the circular median test
        p_val = circ_cmtest(subset_angles, subset_groups);
        
        % Store in symmetric matrix
        p_matrix_cm(i, j) = p_val*num_comparisons;
        p_matrix_cm(j, i) = p_val*num_comparisons;

        p_val = circ_wwtest(subset_angles, subset_groups);
        
        % Store  p-value (with Bonferroni correction)
        p_matrix_ww(i, j) = p_val * num_comparisons;
        p_matrix_ww(j, i) = p_val * num_comparisons;

    end
end
%%
% p_value_matrix = p_matrix_cm;

p_value_matrix = p_matrix_ww;
% 3. Plotting the results
figure('units','normalized','outerposition',[0 0 1 1],'Color', 'w');
% Log-scale makes small p-values (significant) easier to see
h = heatmap(group_ids, group_ids, p_value_matrix, ...
    'Colormap', flipud(hot), ... % Significant values in red/yellow
    'ColorLimits', [0 0.1], ... % Cap at 0.1 to highlight significance
    'MissingDataColor', [0.5 0.5 0.5]);

h.Title = sprintf('Pairwise circ\\_cmtest (\\alpha_{adj} = %.4f)', alpha_adj);
h.XLabel = 'Group ID';
h.YLabel = 'Group ID';
if var2compare==4
h.XDisplayLabels=group_labels(find(counts_per_condition(:,1)>min_count));
h.YDisplayLabels=group_labels(find(counts_per_condition(:,1)>min_count));
elseif var2compare==3
    h.XDisplayLabels={'Trough','Peak','Peak2Trough'};
h.YDisplayLabels={'Trough','Peak','Peak2Trough'};
end
set(h, 'InnerPosition', [0.15 0.15 0.7 0.7]); 

% 3. Alternative: If you want the whole chart object to stay square
% we can set the 'Units' to 'pixels' and make width and height equal
h.Units = 'pixels';
current_pos = h.Position;
side_length = min(current_pos(3:4)); % Find the smaller dimension
h.Position = [current_pos(1) current_pos(2) side_length side_length];

%%

% 1. Identify which groups are directional (Rayleigh Test)
is_directional = zeros(num_groups, 1);
means = nan(num_groups,1);
for i = 1:num_groups
    group_data = angles(groups == group_ids(i));
    % Check if group has enough data and is not uniform
    if ~isempty(group_data)
        means(i) = circ_mean(group_data)
        p_rayleigh = circ_rtest(group_data); 
        if p_rayleigh < 0.05
            is_directional(i) = 1; 
        end
    end
end

% 2. Create the "Filtered" variables for plotting and clustering
valid_idx = find(is_directional == 1);

if isempty(valid_idx)
    error('No groups were significantly directional.');
end

% Extract means for the directional groups only
num_filtered = length(valid_idx);
filtered_means = zeros(num_filtered, 1);
for i = 1:num_filtered
    current_actual_id = group_ids(valid_idx(i));
    group_data = angles(groups == current_actual_id);
    filtered_means(i) = circ_mean(group_data);
end

% Extract the corresponding labels from your heatmap/display labels
filtered_labels = h.XDisplayLabels(valid_idx);

% Ensure labels are in a cell array format for the text functions
if ~iscell(filtered_labels)
    filtered_labels = cellstr(filtered_labels);
end

fprintf('Found %d directional groups out of %d total.\n', num_filtered, num_groups);

%%
% 1. Pre-calculate circular SD for each group correctly
group_stds = zeros(num_filtered, 1);
for i = 1:num_filtered
    % Use group_ids(valid_idx(i)) to get the actual ID value
    current_actual_id = group_ids(valid_idx(i));
    group_data = angles(groups == current_actual_id);
    group_stds(i) = circ_std(group_data); 
end

% 2. Calculate the "Effect Size" Distance Matrix
dist_matrix_effect = zeros(num_filtered);
for i = 1:num_filtered
    for j = i+1:num_filtered
        ang_diff = abs(circ_dist(filtered_means(i), filtered_means(j)));
        pooled_std = sqrt((group_stds(i)^2 + group_stds(j)^2) / 2);
        
        % Avoid division by zero if data is perfectly aligned
        if pooled_std == 0, d_circ = 0; else, d_circ = ang_diff / pooled_std; end
        
        dist_matrix_effect(i,j) = d_circ;
        dist_matrix_effect(j,i) = d_circ;
    end
end

% 3. Cluster using the effect size matrix
Z_effect = linkage(squareform(dist_matrix_effect), 'ward');

% 4. Plotting (Make sure to use Z_effect here!)
figure('Color', 'w');
% Ensure filtered_labels is a cell array of strings
if isstring(filtered_labels) || ischar(filtered_labels)
    filtered_labels = cellstr(filtered_labels);
end

[H, T, outperm] = dendrogram(Z_effect, 'Labels', filtered_labels);
xtickangle(45);
title('Clustering Areas by Circular Effect Size (Directional Groups Only)');
ylabel('Standardized Distance (d_c)');
grid on;

%%

% 1. Identify which groups to keep (from previous step)
valid_idx = find(is_directional == 1);

% p_matrix =p_matrix_cm;
p_matrix= p_matrix_ww;
filtered_p_matrix = p_matrix(valid_idx, valid_idx);
filtered_group_labels = h.XDisplayLabels(valid_idx); 

% 2. Get the reordering indices from the dendrogram
% (Make sure you ran the dendrogram code from the previous response first)
% [H, T, outperm] = dendrogram(Z, 'Labels', filtered_labels);

% 3. Reorder the p-value matrix and labels based on the dendrogram
reordered_p_matrix = filtered_p_matrix(outperm, outperm);
reordered_labels   = filtered_group_labels(outperm);

% 4. Plot the Rearranged Heatmap
figure('Color', 'w');
h_reorder = heatmap(reordered_labels, reordered_labels, reordered_p_matrix, ...
    'Colormap', flipud(hot), ... 
    'ColorLimits', [0 0.1], ... 
    'MissingDataColor', [0.5 0.5 0.5]);

h_reorder.Title = 'Pairwise circ\_cmtest (Ordered by Angular Similarity)';
h_reorder.XLabel = 'Group (Clustered)';
h_reorder.YLabel = 'Group (Clustered)';

% Note: We use the reordered_labels directly as the axis labels 
% so the Heatmap and Dendrogram match perfectly.


%% plot means withing circle with grouping colors
% 1. Manual Entry: Set your desired level/number of clusters based on the dendrogram
num_clusters = 2; 

% 2. Assign each area to a cluster
% This "cuts" the linkage tree (Z_effect) into the number of groups you chose
group_assignments = cluster(Z_effect, 'MaxClust', num_clusters);

% 3. Generate Colors
% 'lines' provides high-contrast colors; 'colorcube' or 'jet' are alternatives
cluster_colors = lines(num_clusters); 

% 4. Print results to Command Window to verify
fprintf('\n--- Cluster Assignments (Level: %d) ---\n', num_clusters);
for c = 1:num_clusters
    members = filtered_labels(group_assignments == c);
    fprintf('Cluster %d (%d areas): %s\n', c, length(members), strjoin(members, ', '));
end

figure('Color', 'w', 'Name', 'Clustered Mean Angles');
polar_ax = polaraxes;
hold on;

% --- NEW ADDITION: Plot non-directional groups as faint background lines ---
% We loop through all groups and only plot those NOT in valid_idx
for i = 1:num_groups
    if ~ismember(i, valid_idx)
        current_angle = means(i);
        % Plot as a very thin, light gray dotted line
        polarplot([current_angle current_angle], [0 1], ...
            'Color', [0.85 0.85 0.85], 'LineWidth', 1, 'LineStyle', ':');
    end
end

% --- YOUR EXISTING WORKING CODE: Plot each directional group ---
for i = 1:num_filtered
    current_angle = filtered_means(i);
    current_cluster = group_assignments(i);
    current_color = cluster_colors(current_cluster, :);
    current_label = filtered_labels{i};
    
    % Draw the needle
    polarplot([current_angle current_angle], [0 1], ...
        'Color', current_color, 'LineWidth', 2.5);
    
    % Draw the tip marker
    polarplot(current_angle, 1, 'o', ...
        'MarkerFaceColor', current_color, 'MarkerEdgeColor', 'w', 'MarkerSize', 6);
        
    % --- RADIAL TEXT LOGIC ---
    angle_deg = rad2deg(current_angle);
    text_rot = angle_deg; 
    if angle_deg > 90 && angle_deg < 270
        text_rot = angle_deg + 180; 
    end
    
    r_pos = 1.1 + (mod(i,3) * 0.1); 

    text(current_angle, r_pos, current_label, ...
        'Color', current_color, 'FontSize', 8, 'FontWeight', 'bold', ...
        'Rotation', -text_rot + 90, ... 
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% 3. Plot 

% 4. Formatting
title(sprintf('Mean Angles: %d Clusters', num_clusters), 'FontSize', 12);
polar_ax.ThetaZeroLocation = 'top'; 
polar_ax.ThetaDir = 'clockwise';
rlim([0 1.6]); % Increased limit to accommodate radial text
polar_ax.RTickLabels = {}; % Clean up the plot by removing radial numbers
grid on;


%% add signrfficicant differences
% p_matrix =p_matrix_cm;
p_matrix= p_matrix_ww;
% --- Add this after your needle plotting loop ---

% 1. Settings for significance arcs
sig_threshold = 0.05;
arc_r_start = 1.6; % Start drawing arcs outside the labels
arc_step = 0.15;   % Distance between stacked arcs
comparison_count = 0;

% 2. Loop through the p_matrix to find significant pairs
% Note: Using p_matrix(valid_idx, valid_idx) to match the directional groups
p_sub = p_matrix(valid_idx, valid_idx);

for i = 1:num_filtered
    for j = i+1:num_filtered
        p_val = p_sub(i,j);
        
        if p_val < sig_threshold
            comparison_count = comparison_count + 1;
            
            % Determine arc radial position
            r_arc = arc_r_start + (comparison_count * arc_step);
            
            % Get the two angles
            ang1 = filtered_means(i);
            ang2 = filtered_means(j);
            
            % Calculate the shortest path for the arc
            arc_angles = linspace(ang1, ang1 + circ_dist(ang2, ang1), 50);
            
            % Draw the arc
            polarplot(arc_angles, repmat(r_arc, size(arc_angles)), ...
                '-k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            
            % Add asterisks at the midpoint of the arc
            mid_angle = ang1 + circ_dist(ang2, ang1)/2;
            
            if p_val < 0.001
                txt = '***';
            elseif p_val < 0.01
                txt = '**';
            else
                txt = '*';
            end
            
            text(mid_angle, r_arc + 0.05, txt, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', 10, 'FontWeight', 'bold');
        end
    end
end

% Adjust radial limit to fit the arcs
rlim([0, arc_r_start + (comparison_count + 1) * arc_step]);
%% same plot but now in tiome

%% Linear Sinusoidal Phase Plot with Significance Brackets
%% Supplementary plot 2: Sinusoidal Plot (Corrected Alignment & Staggering) 

% 1. Setup Phase Space
phase_x = linspace(-pi, pi, 500); 
sine_y  = cos(phase_x);           

figure('Color', 'w', 'Position', [100 100 1000 700]);
hold on;

% 2. Plot Reference Sine Wave
plot(phase_x, sine_y, 'k', 'LineWidth', 2, 'Color', [0.8 0.8 0.8]); 

% --- NEW ADDITION: Background Non-Directional Lines ---
for i = 1:num_groups
    if is_directional(i) == 0
        current_angle = circ_dist(means(i), 0); % Align to -pi to pi
        line([current_angle current_angle], [-1.1 1.1], ...
            'Color', [0.85 0.85 0.85], 'LineWidth', 0.8, 'LineStyle', ':');
    end
end

% --- UPDATED FOREGROUND: Directional Clusters only ---
for i = 1:num_filtered
    current_angle = circ_dist(filtered_means(i), 0); % Align to -pi to pi
    current_cluster = group_assignments(i);
    current_color = cluster_colors(current_cluster, :);
    current_label = filtered_labels{i};
    
    % 3. Draw Vertical Needle
    line([current_angle current_angle], [-1.1 1.1], ...
        'Color', current_color, 'LineWidth', 2.5);
    
    % 4. Intersection Marker
    plot(current_angle, cos(current_angle), 'o', ...
        'MarkerFaceColor', current_color, 'MarkerEdgeColor', 'w', 'MarkerSize', 8);
    
    % 5. Vertical Staggered Labels
    % We use 'i' from the num_filtered loop to stagger correctly
    y_pos = 1.2 + (mod(i,3) * 0.2); 
    text(current_angle, y_pos, current_label, ...
        'Color', current_color, 'FontSize', 8, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'Rotation', 90);
end

% 6. Significance Brackets (Linear)
p_sub = p_matrix(valid_idx, valid_idx);
bracket_y_start = 2.0; % Increased height to clear labels
bracket_step = 0.15;
comparison_count = 0;

y_level = bracket_y_start; % Initialize for the ylim later

for i = 1:num_filtered
    for j = i+1:num_filtered
        p_val = p_sub(i,j);
        if p_val < 0.05
            comparison_count = comparison_count + 1;
            y_level = bracket_y_start + (comparison_count * bracket_step);
            
            % Use filtered_means directly (short path)
            ang1 = circ_dist(filtered_means(i), 0);
            ang2 = circ_dist(filtered_means(j), 0);
            
            % Linear Bracket
            plot([ang1, ang2], [y_level, y_level], '-k', 'LineWidth', 1.2);
            plot([ang1, ang1], [y_level, y_level-0.05], '-k');
            plot([ang2, ang2], [y_level, y_level-0.05], '-k');
            
            % Stars
            if p_val < 0.001, txt = '***'; elseif p_val < 0.01, txt = '**'; else, txt = '*'; end
            text(mean([ang1, ang2]), y_level + 0.05, txt, ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
        end
    end
end

% 7. Formatting
xlabel('Phase (radians)');
ylabel('Amplitude');
set(gca, 'XTick', [-pi, -pi/2, 0 pi/2 2*pi/3 pi], 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '2\pi/3','\pi'});
xlim([-pi*1.1, pi*1.1]);
ylim([-1.2, y_level + 0.5]); 
grid on;
title(sprintf('Phase Relationship (%d Clusters)', num_clusters));
%%
print(gcf,'-vector','-dsvg',[figure_folder,'/sine plot with coincident phases only significnat combs zoom in.svg'])