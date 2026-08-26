%% Figure 3 — example cross-correlograms
% Standalone plots from ANALYZE cross correlogram all aniamls.m
% (former sections 15–16): browse trough–trough CCGs for one area pair,
% then recompute a chosen pair at higher time precision.
%
% Needs Supporting codes/GENERATE_CROSS_CORR_EXAMPLE.m on the MATLAB path.

addpath('\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 3\Supporting codes');
data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';

%% 1 Load cross-correlograms
% Per-session structs from Estimate cross correlogram all aniamls
% (play, non-play, and percentile/shuffle CCGs). Local copies under Data\.

npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
cc_folder = [data_root, '\Analysis results\Cross_correlogram'];
neuron_folder = [data_root, '\Analysis results\phase locking data'];

bin_size_cc = 0.01;
hist_range_cc = [-3 3];
psth_edges_cc = hist_range_cc(1):bin_size_cc:hist_range_cc(2);
time_centers = .5*(psth_edges_cc(1:end-1)+psth_edges_cc(2:end));

load([cc_folder,'\cross_corr_struct_structure_updated_with_non_playbouts_3sec.mat'],'cross_corr_struct');
load([cc_folder,'\cross_corr_struct_animal_names_updated_with_non_playbouts_3sec.mat'],'animal_names');

%% 2 Merge pairs across sessions
% Concatenate CCGs, shuffle CCGs, area labels, session, and cluster IDs.

play_cc = [];
non_play_cc = [];
all_area_comb_cc = [];
session_id_cc = [];
clusters_id_cc = [];
all_cross_corr_play_pctl = [];
all_cross_corr_non_play_pctl = [];
for fn = 1:numel(cross_corr_struct)
    all_cross_corr_play_pctl     = [all_cross_corr_play_pctl;cross_corr_struct(fn).all_cross_corr_play_pctl];
    all_cross_corr_non_play_pctl = [all_cross_corr_non_play_pctl;cross_corr_struct(fn).all_cross_corr_non_play_pctl];
    play_cc                      = cat(1,play_cc,cross_corr_struct(fn).all_cross_corr_play);
    non_play_cc                  = cat(1,non_play_cc,cross_corr_struct(fn).all_cross_corr_non_play);
    all_area_comb_cc             = [all_area_comb_cc;cross_corr_struct(fn).these_neurons_areas(cross_corr_struct(fn).cross_corr_comb)];
    session_id_cc                = [session_id_cc;repmat(animal_names(fn,1),size(cross_corr_struct(fn).cross_corr_comb,1),1)];
    clusters_id_cc               = [clusters_id_cc;cross_corr_struct(fn).cluster_info.cluster_id(cross_corr_struct(fn).cross_corr_comb)];
end

%% 3 Area–area combinations
% Unique directed pairs among SupCol, DLPAG, LPAG, VLPAG, DR.

area_combinations = cell2table(all_area_comb_cc);
area_combinations.Properties.VariableNames = {'Neuron1','Neuron2'};
unique_area_combinations = unique(area_combinations, 'rows');

relevant_areas = {'SupCol','DLPAG','LPAG','VLPAG','DR'};
unique_area_combinations = unique_area_combinations(ismember(unique_area_combinations.Neuron1,relevant_areas) & ismember(unique_area_combinations.Neuron2,relevant_areas),:);

%% 4 Trough / peak / unlocked labels
% Map each CCG pair onto all_neurons_TD, then label neurons by delta-phase
% locking. Sections 5–7 use trough–trough pairs.

load([neuron_folder,'\theta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) = {'isRt'};
all_neurons_TD = all_neurons;
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner1'}))      = {'ThetaPartner1'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner2'}))      = {'ThetaPartner2'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Play'}))          = {'ThetaPlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'PrePlay'}))       = {'ThetaPrePlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'EntireSession'})) = {'ThetaEntireSession'};

load([neuron_folder,'\delta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) = {'isRt'};
all_neurons_TD.DeltaPartner1      = all_neurons.Partner1;
all_neurons_TD.DeltaPartner2      = all_neurons.Partner2;
all_neurons_TD.DeltaEntireSession = all_neurons.EntireSession;
all_neurons_TD.DeltaPlay          = all_neurons.Play;
all_neurons_TD.DeltaPrePlay       = all_neurons.PrePlay;

nRows = size(clusters_id_cc,1);
idx_pairs_cc = zeros(nRows,2);
for i = 1:nRows
    s = session_id_cc{i};
    id_pair = clusters_id_cc(i,:);
    for j = 1:2
        idx_pairs_cc(i,j) = find( ...
            strcmp(all_neurons_TD.session, s) & ...
            all_neurons_TD.cluster_id == id_pair(j), ...
            1);
    end
end

trough = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);

%% 5 Prepare trough CCGs for one area–area pair
% j indexes unique_area_combinations (e.g. 14 = LPAG–SupCol).
% Play CCGs vs shuffle; switch the commented lines to use non-play.

cc2analysie = play_cc;
cc2analyse_control = all_cross_corr_play_pctl;
% cc2analysie = non_play_cc;
% cc2analyse_control = all_cross_corr_non_play_pctl;

j = 14; 
% j = 20; % SupCol and VLPAG

disp('Estimating cross correlogram for the areas : ' )
disp([unique_area_combinations.Neuron1{j}, ' and ', unique_area_combinations.Neuron2{j}])
smoth_wind = 1;
area_combination_indexes = ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j));
cell_type_comb_trough = find(trough(idx_pairs_cc(:,1)) & trough(idx_pairs_cc(:,2)) & area_combination_indexes);

trough_this_cc_control = cc2analyse_control(cell_type_comb_trough,:);
trough_this_cc         = squeeze(cc2analysie(cell_type_comb_trough,:));
for k = 1:size(trough_this_cc,1)
    trough_this_cc(k,:) = movmean(trough_this_cc(k,:),smoth_wind);
end
for k = 1:size(trough_this_cc_control,1)
    trough_this_cc_control(k,:) = movmean(trough_this_cc_control(k,:),smoth_wind);
end

%% 6 Browse trough-cell pairs
% Grid of trough–trough CCGs (real vs shuffle) with mean CCG > 1, to pick
% a pair for section 7. The subplot ylabel is the plotting index for kk.

coincident_level = 0.005;
x_lim = [-.25 .25];

figure('units','normalized','outerposition',[0 0 .5 1]);
nsp = 1;
list_of_neuron = find(mean(trough_this_cc,2)>1)';
fig_n = 1;

for kk_n = 1:numel(list_of_neuron)
    kk = list_of_neuron(kk_n);
    if nsp > 25
        nsp = 1;
        sgtitle(num2str(fig_n))
        fig_n = fig_n + 1;
        figure('units','normalized','outerposition',[0 0 .5 1]);
    end
    subplot(5,5,nsp)
    index1 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_trough(kk),1) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_trough(kk)}));
    index2 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_trough(kk),2) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_trough(kk)}));
    if ~isempty(index1) & ~isempty(index2)
        plot(time_centers, trough_this_cc(kk,:), 'b')
        hold on
        plot(time_centers, trough_this_cc_control(kk,:), 'k')
        xlim(x_lim)
        title([all_neurons_TD.session{index1}, ' ID', num2str(all_neurons_TD.cluster_id(index1)), ' cc#', num2str(kk_n)])
        ylabel({[' ID ',num2str(all_neurons_TD.cluster_id(index2))],['Ploting index: ' num2str(kk_n)]})
        ylim tight
        y_lim = ylim;
        fill([-1 1 1 -1]*coincident_level,y_lim([1 1 2 2]), 'g', 'FaceAlpha',.5, 'EdgeColor','none')
    end
    nsp = nsp + 1;
end

%% 7 Fig 3i, and Supp Fig 5g Plot a chosen pair (GENERATE_CROSS_CORR_EXAMPLE)
% Uses list_of_neuron(plotting index) from section 6: session + cluster IDs
% → example CCG at higher time precision.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% wrote ploting index as selected among all subplots in previous section %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% some examples that work %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk = list_of_neuron(65); % Example in the paper, you need to obtain the right cross correlogram
% from the panesl before
%%% you can chek the neuron id of these neurons, and check the kilosort
%%% output indexes to identify the neuron and trhen the channel where they
%%% belong to and then the aproxiamte location in the midbrain

% kk = list_of_neuron(30); %%  for j = 20 in previouse section, in case you
% want to plot a different combination
time_precision = 0.005;

bin_size_cc = 0.01;
hist_range_cc = [-3 3];

index1 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_trough(kk),1) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_trough(kk)}));
index2 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_trough(kk),2) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_trough(kk)}));
animal_id = all_neurons_TD.session{index1};
id_1 = all_neurons_TD.cluster_id(index1);
id_2 = all_neurons_TD.cluster_id(index2);

animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];

animal_file_names = cellfun(@(x) ['B', x],strsplit([animal_list.name], 'B'), 'UniformOutput',false)';
animal_file_names(1) = [];

animal2analize = animal_id;
animal_list = animal_list(ismember(animal_file_names,animal2analize));
GENERATE_CROSS_CORR_EXAMPLE([npx_Raw_Data, '\', animal_list.name],bin_size_cc, hist_range_cc,[id_1 id_2],time_precision);
