%% Figure 3i–j — coincidence-spike phases
% Polar histograms of coincident-spike LFP phase for trough–trough, peak–peak,
% and unlocked–unlocked pairs. Only the arrays used by the two plots below
% are built: all_angles (pooled events) and angle_per_condition_after
% (circular mean per area pair, after bout onset).
%
% Coincidence structs and delta locking labels are under Data\Analysis results.

%% 1 Load coincidence events and delta locking labels
data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
coincidence_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\Cross_correlogram';
phase_folder = [data_root, '\Analysis results\phase locking data'];
freq_tag = num2str([1 5]);   % delta band used when GENERATE saved the file

load([coincidence_folder, '\ALL AANIMALS_FreqRange_', freq_tag, '_phase_coincidence_structure.mat'], 'phase_struct');
load([coincidence_folder, '\ALL_ANIMALS_FreqRange_', freq_tag, '_phase_coincidence_structure_animal_names.mat'], 'animal_names');

load([phase_folder, '\delta_all_neurons_v2.mat'], 'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) = {'isRt'};

%% 2 Mean coincidence phase per pair (neuron 1), before vs after bout onset
% sync_event_data columns: [spike_rel, phase_on_ch1, phase_on_ch2, period]
% period==1 / ==2 kept as in ANALYZE_phase_coincidence (called before / after).
synch_spikes_phases_before = [];
synch_spikes_phases_after  = [];
all_area_comb = [];
session_id = [];
clusters_id = [];
cluster_fr = [];

for fn = 1:numel(phase_struct)
    n_pairs = size(phase_struct(fn).sync_event_data, 1);
    for j = 1:n_pairs
        row_before = nan(1, 4);
        row_after  = nan(1, 4);
        ev1 = phase_struct(fn).sync_event_data{j, 1};
        ev2 = phase_struct(fn).sync_event_data{j, 2};
        if ~isempty(ev1)
            before1 = ev1(:, 4) == 1;
            before2 = ev2(:, 4) == 1;
            after1  = ev1(:, 4) == 2;
            after2  = ev2(:, 4) == 2;
            if any(before1) && any(before2)
                row_before = [circ_mean(ev1(before1, [2 3])) circ_mean(ev2(before2, [2 3]))];
            elseif any(before1)
                row_before = [circ_mean(ev1(before1, [2 3])) [NaN NaN]];
            elseif any(before2)
                row_before = [[NaN NaN] circ_mean(ev2(before2, [2 3]))];
            end
            if any(after1) && any(after2)
                row_after = [circ_mean(ev1(after1, [2 3])) circ_mean(ev2(after2, [2 3]))];
            elseif any(after1)
                row_after = [circ_mean(ev1(after1, [2 3])) [NaN NaN]];
            elseif any(after2)
                row_after = [[NaN NaN] circ_mean(ev2(after2, [2 3]))];
            end
        end
        synch_spikes_phases_before = [synch_spikes_phases_before; row_before];
        synch_spikes_phases_after  = [synch_spikes_phases_after;  row_after];
    end
    all_area_comb = [all_area_comb; phase_struct(fn).cluster_info.area(phase_struct(fn).synch_comb)];
    session_id    = [session_id; repmat(animal_names(fn, 1), n_pairs, 1)];
    clusters_id   = [clusters_id; phase_struct(fn).cluster_info.cluster_id(phase_struct(fn).synch_comb)];
    cluster_fr    = [cluster_fr; phase_struct(fn).cluster_info.fr(phase_struct(fn).synch_comb)];
end

%% 3 Unique area combinations (treat A–B and B–A as one)
area_combinations = cell2table(all_area_comb);
area_combinations.Properties.VariableNames = {'Neuron1', 'Neuron2'};
unique_area_combinations = unique(area_combinations, 'rows');

relevant_areas = {'SupCol', 'DLPAG', 'LPAG', 'VLPAG', 'DR'};
unique_area_combinations = unique_area_combinations( ...
    ismember(unique_area_combinations.Neuron1, relevant_areas) & ...
    ismember(unique_area_combinations.Neuron2, relevant_areas), :);

data = string(table2array(unique_area_combinations));
sortedRows = sort(data, 2);
[~, ~, ic] = unique(sortedRows, 'rows');
rowGroups = accumarray(ic, (1:size(data, 1))', [], @(x) {sort(x')});
mirroredPairs = rowGroups(cellfun(@(x) length(x) > 1, rowGroups));
for j = 1:size(mirroredPairs, 1)
    unique_area_combinations{mirroredPairs{j}, :} = repmat(unique_area_combinations{mirroredPairs{j}(1), :}, 2, 1);
end
unique_area_combinations = unique(unique_area_combinations, 'rows');

%% 4 Match each coincidence pair to delta peak / trough / unlocked labels
nRows = size(clusters_id, 1);
idx_pairs = zeros(nRows, 2);
for i = 1:nRows
    s = session_id{i};
    id_pair = clusters_id(i, :);
    for j = 1:2
        hit = strcmp(all_neurons.session, s) & all_neurons.cluster_id == id_pair(j);
        if any(hit)
            idx_pairs(i, j) = find(hit, 1);
        else
            disp('Missing neuron')
            idx_pairs(i, j) = NaN;
        end
    end
end

trough        = all_neurons.EntireSession.PPCPval < 0.01 & (all_neurons.EntireSession.PreferedAngle > pi/2 | all_neurons.EntireSession.PreferedAngle < -pi/2);
peak          = all_neurons.EntireSession.PPCPval < 0.01 & ~(all_neurons.EntireSession.PreferedAngle > pi/2 | all_neurons.EntireSession.PreferedAngle < -pi/2);
unlocked      = all_neurons.EntireSession.PPCPval > .1;
min_fr        = 0;

%% 5 Arrays used by the plots: pooled phases and per-area-pair means after onset
% all_angles columns: [phase, unused, cell-type (1 trough / 2 peak / 3 unlocked), area-pair index]
angle_per_condition_after = nan(size(unique_area_combinations, 1), 3);
all_angles = [];

for j = 1:size(unique_area_combinations, 1)
    area_combination_indexes = ...
        (ismember(area_combinations.Neuron1, unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2, unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2, unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1, unique_area_combinations.Neuron2(j)));

    trough_trough = (trough(idx_pairs(:, 1)) & trough(idx_pairs(:, 2)) & area_combination_indexes & mean(cluster_fr, 2) > min_fr);
    peak_peak     = (peak(idx_pairs(:, 1)) & peak(idx_pairs(:, 2)) & area_combination_indexes & mean(cluster_fr, 2) > min_fr);
    unlocked_unlocked = (unlocked(idx_pairs(:, 1)) & unlocked(idx_pairs(:, 2)) & area_combination_indexes & mean(cluster_fr, 2) > min_fr);

    if ~isempty(trough_trough)
        both = trough_trough & ~isnan(synch_spikes_phases_before(:, 1)) & ~isnan(synch_spikes_phases_after(:, 1));
        alpha_before = synch_spikes_phases_before(both, 1);
        alpha_after  = synch_spikes_phases_after(both, 1);
        all_angles = [all_angles; [[alpha_before; alpha_after], [alpha_before*0; alpha_after*0+1], [alpha_before; alpha_after]*0+1, [alpha_before; alpha_after]*0+j]];
        angle_per_condition_after(j, 1) = circ_mean(alpha_after);
    end

    if ~isempty(peak_peak)
        both = peak_peak & ~isnan(synch_spikes_phases_before(:, 1)) & ~isnan(synch_spikes_phases_after(:, 1));
        alpha_before = synch_spikes_phases_before(both, 1);
        alpha_after  = synch_spikes_phases_after(both, 1);
        all_angles = [all_angles; [[alpha_before; alpha_after], [alpha_before*0; alpha_after*0+1], [alpha_before; alpha_after]*0+2, [alpha_before; alpha_after]*0+j]];
        angle_per_condition_after(j, 2) = circ_mean(alpha_after);
    end

    if ~isempty(unlocked_unlocked)
        both = unlocked_unlocked & ~isnan(synch_spikes_phases_before(:, 1)) & ~isnan(synch_spikes_phases_after(:, 1));
        alpha_before = synch_spikes_phases_before(both, 1);
        alpha_after  = synch_spikes_phases_after(both, 1);
        all_angles = [all_angles; [[alpha_before; alpha_after], [alpha_before*0; alpha_after*0+1], [alpha_before; alpha_after]*0+3, [alpha_before; alpha_after]*0+j]];
        angle_per_condition_after(j, 3) = circ_mean(alpha_after);
    end
end

%% Fig 3j: as on the papaer (or almost)
figure

% Define groups and colors
groups = [2 1 3];
colors = {'r','b','g'};


for i = 1:length(groups)

    % Get angles for this group
    angles = all_angles(all_angles(:,3)==groups(i),1);

    % Rayleigh test
    [p, z] = circ_rtest(angles);

    % Circular mean
    mean_angle = circ_mean(angles);

    % Plot histogram
    if i<3
    polarhistogram(angles, -pi:pi/16:pi, ...
        'FaceColor', colors{i}, ...
        'Normalization', 'pdf');
    else
        polarhistogram(angles, -pi:pi/3:pi, ...
        'FaceColor', colors{i}, ...
        'Normalization', 'pdf');
    end


    hold on

    % Get current radial limit
    rlim_current = rlim;

    % Plot mean-angle line
    polarplot([mean_angle mean_angle], ...
              [0 rlim_current(2)], ...
              'Color', colors{i}, ...
              'LineWidth', 2);

    % Add p-value as text
    % Position slightly outside the histogram
    text(mean_angle, rlim_current(2)*1.05, ...
        sprintf('p = %.3g', p), ...
        'Color', colors{i}, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10);
end

hold off

%% Fig 3 j: alternative, mean of each pair before or after play bout onset
figure

colors = {'b','r','g'};

for type = 1:size(angle_per_condition_after,2)

    % Get angles and remove NaNs
    angles = angle_per_condition_after(:,type);
    % angles = angle_per_condition_before(:,type);

     
    angles = angles(~isnan(angles));

    % Skip if there are no data
    if isempty(angles)
        continue
    end

    % Rayleigh test
    [p, z] = circ_rtest(angles);

    % Circular mean
    mean_angle = circ_mean(angles);

    % Plot histogram
    polarhistogram(angles, -pi:pi/16:pi, ...
        'FaceColor', colors{type}, ...
        'Normalization', 'pdf');
    hold on


    % Get radial limit
    rmax = rlim;
    
    % Plot mean angle
    polarplot([mean_angle mean_angle], ...
              [0 rmax(2)], ...
              'Color', colors{type}, ...
              'LineWidth', 2);

    % Add mean angle and p-value
    text(mean_angle, rmax(2)*1.05, ...
        sprintf('\\mu = %.1f°, p = %.3g', ...
        rad2deg(mean_angle), p), ...
        'Color', colors{type}, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10);
end

hold off
