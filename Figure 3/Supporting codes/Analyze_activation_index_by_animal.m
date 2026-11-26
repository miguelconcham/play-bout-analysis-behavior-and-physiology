%% Animal-balanced activation index for SupCol and LPAG
%
% Descriptive analysis only. Neurons are compared only with neurons from
% the same biological animal and area. Animals are then given equal weight.
% No pooled-neuron or animal-level inferential p-values are calculated.
%
% Required checkpoints, produced by
% Estiamte_activation_index_and_other_variables.m:
%   activation_index_A1_A7.mat
%   activation_index_section6.mat

clearvars
close all
clc

%% USER SETTINGS
% Set an absolute folder here. Leave empty to use the default under
% Data/Analysis results/activation index/.
output_folder = '';

selected_area_names = ["SupCol", "LPAG"];

% One common behavior definition is used for both the behavior heatmaps and
% the play/non-play averages. These are the behavior indices displayed in
% section 4 of Figure 3/Figure 3g-h Activation Index.m.
play_behavior_indices = [4 6 8 14 18 20 24 5 7 9 15 19 21 25];
nonplay_behavior_indices = [10 12 16 22 26 34 11 13 17 23 27 35];

settings = struct;
settings.n_permutations = 10000;
settings.minimum_type_neurons = 3;
settings.minimum_reference_neurons = 3;
settings.minimum_behaviors_per_condition = 3;
settings.response_window_seconds = [0 5];
settings.random_seed = 8613;
settings.resume_existing_calculation = true;
settings.figure_visible = 'on';
settings.jpeg_resolution = 300;
settings.activation_color_limits = [-4 4];
settings.paper_summary_y_limits = [-2 2];

%% Paths and inputs
experiment_root = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology'];
activation_folder = fullfile(experiment_root, 'Data', ...
    'Analysis results', 'activation index');
animal_checkpoint_file = fullfile(activation_folder, ...
    'activation_index_A1_A7.mat');
pooled_checkpoint_file = fullfile(activation_folder, ...
    'activation_index_section6.mat');

if isempty(output_folder)
    output_folder = fullfile(activation_folder, ...
        'animal level SupCol LPAG');
end
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

if ~isfile(animal_checkpoint_file)
    error('Missing animal-level input checkpoint: %s', ...
        animal_checkpoint_file);
end
if ~isfile(pooled_checkpoint_file)
    error('Missing pooled activation-index checkpoint: %s', ...
        pooled_checkpoint_file);
end

animal_data = load(animal_checkpoint_file);
required_animal_variables = {'all_activation_order','area_list', ...
    'psth_list'};
missing_variables = required_animal_variables( ...
    ~isfield(animal_data, required_animal_variables));
if ~isempty(missing_variables)
    error('Animal checkpoint is missing: %s', ...
        strjoin(missing_variables, ', '));
end

all_activation_order = animal_data.all_activation_order;
source_area_names = string(animal_data.area_list(:));
n_source_behaviors = size(all_activation_order,1);
n_time = size(all_activation_order{1,1,4},2);
behavior_names = activation_behavior_names(animal_data, ...
    n_source_behaviors);
cell_type_names = ["Peak-locked"; "Trough-locked"; "Unlocked"];
n_types = numel(cell_type_names);

if isfield(animal_data, 'warped_time') ...
        && numel(animal_data.warped_time) == n_time
    activation_time = animal_data.warped_time(:)';
elseif isfield(animal_data, 'time2use') ...
        && numel(animal_data.time2use) == n_time
    activation_time = animal_data.time2use(:)';
else
    error('No checkpoint time vector matches the stored PSTHs.');
end

behavior_indices = [play_behavior_indices, nonplay_behavior_indices];
if any(behavior_indices < 1 | behavior_indices > n_source_behaviors)
    error('A selected behavior index is outside the checkpoint.');
end
if numel(unique(behavior_indices)) ~= numel(behavior_indices)
    error('Play and non-play behavior definitions overlap.');
end
n_play_behaviors = numel(play_behavior_indices);
n_nonplay_behaviors = numel(nonplay_behavior_indices);
n_selected_behaviors = numel(behavior_indices);
selected_behavior_names = behavior_names(behavior_indices);
play_positions = 1:n_play_behaviors;
nonplay_positions = n_play_behaviors + (1:n_nonplay_behaviors);

selected_area_index = nan(numel(selected_area_names),1);
for area_number = 1:numel(selected_area_names)
    match = find(source_area_names == selected_area_names(area_number));
    if numel(match) ~= 1
        error('Area %s was not found exactly once.', ...
            selected_area_names(area_number));
    end
    selected_area_index(area_number) = match;
end
n_areas = numel(selected_area_names);

animal_file_info = dir(animal_checkpoint_file);
input_signature = struct;
input_signature.file_bytes = animal_file_info.bytes;
input_signature.file_datenum = animal_file_info.datenum;
input_signature.source_area_names = source_area_names;
input_signature.behavior_names = behavior_names;
input_signature.selected_area_names = selected_area_names;
input_signature.selected_behavior_names = selected_behavior_names;
input_signature.activation_time = activation_time;

response_index = activation_time >= ...
    settings.response_window_seconds(1) ...
    & activation_time <= settings.response_window_seconds(2);
if ~any(response_index)
    error('The requested response window has no samples.');
end

%% Biological-animal map
all_sessions = strings(0,1);
for area_number = 1:n_areas
    source_area = selected_area_index(area_number);
    sessions = string(all_activation_order{1,source_area,3});
    all_sessions = [all_sessions; sessions(:)]; %#ok<AGROW>
end
all_sessions = unique(all_sessions(strlength(all_sessions) > 0), ...
    'stable');
session_to_animal = table;
session_to_animal.Session = all_sessions;
session_to_animal.AnimalID = biological_animal_id(all_sessions);
animal_ids = unique(session_to_animal.AnimalID, 'stable');
n_animals = numel(animal_ids);
writetable(session_to_animal, fullfile(output_folder, ...
    'session_to_biological_animal.csv'));

fprintf('Animal-balanced activation index\n');
fprintf('Areas: %s\n', strjoin(selected_area_names, ', '));
fprintf('Biological animals: %d\n', n_animals);
disp(session_to_animal)

%% Calculate animal-local activation indices (resumable)
calculation_file = fullfile(output_folder, ...
    'animal_activation_index_calculation.mat');
expected_size = [n_areas, n_types, n_selected_behaviors, ...
    n_animals, n_time];

loaded_progress = false;
if settings.resume_existing_calculation && isfile(calculation_file)
    progress = load(calculation_file);
    if isfield(progress, 'animal_activation_index') ...
            && isequal(size(progress.animal_activation_index), ...
            expected_size) ...
            && isfield(progress, 'completed_behavior') ...
            && isequal(size(progress.completed_behavior), ...
            [n_areas, n_selected_behaviors]) ...
            && isfield(progress, 'saved_behavior_indices') ...
            && isequal(progress.saved_behavior_indices, behavior_indices) ...
            && isfield(progress, 'saved_animal_ids') ...
            && isequal(progress.saved_animal_ids, animal_ids) ...
            && isfield(progress, 'saved_settings') ...
            && compatible_calculation_settings( ...
            progress.saved_settings, settings) ...
            && isfield(progress, 'saved_input_signature') ...
            && isequaln(progress.saved_input_signature, input_signature)
        animal_activation_index = progress.animal_activation_index;
        type_neuron_count = progress.type_neuron_count;
        area_neuron_count = progress.area_neuron_count;
        reference_neuron_count = progress.reference_neuron_count;
        completed_behavior = progress.completed_behavior;
        loaded_progress = true;
        fprintf('Resuming existing calculation checkpoint.\n');
    else
        warning(['Existing calculation checkpoint does not match the ' ...
            'current settings and will be replaced.']);
    end
end

if ~loaded_progress
    animal_activation_index = nan(expected_size);
    type_neuron_count = zeros(n_areas, n_types, ...
        n_selected_behaviors, n_animals);
    area_neuron_count = zeros(n_areas, n_types, ...
        n_selected_behaviors, n_animals);
    reference_neuron_count = zeros(n_areas, n_types, ...
        n_selected_behaviors, n_animals);
    completed_behavior = false(n_areas, n_selected_behaviors);
end

for area_number = 1:n_areas
    source_area = selected_area_index(area_number);
    for selected_behavior_number = 1:n_selected_behaviors
        if completed_behavior(area_number, selected_behavior_number)
            continue
        end

        source_behavior = behavior_indices(selected_behavior_number);
        fprintf('Area %s; behavior %d/%d (%s)\n', ...
            selected_area_names(area_number), ...
            selected_behavior_number, n_selected_behaviors, ...
            selected_behavior_names(selected_behavior_number));

        [responses, type_index, sessions] = activation_block( ...
            all_activation_order, source_behavior, source_area, n_time);
        animal_for_row = biological_animal_id(sessions);
        % Use one fixed, complete neuron pool across the entire trace.
        % This prevents changing target/reference membership between bins.
        valid_row = all(isfinite(responses),2);

        for type_number = 1:n_types
            for animal_number = 1:n_animals
                animal_mask = animal_for_row == animal_ids(animal_number) ...
                    & valid_row;
                target_mask = animal_mask & type_index == type_number;
                reference_mask = animal_mask & type_index ~= type_number;

                target_response = responses(target_mask,:);
                area_response = responses(animal_mask,:);
                n_target = size(target_response,1);
                n_area = size(area_response,1);
                n_reference = sum(reference_mask);

                type_neuron_count(area_number,type_number, ...
                    selected_behavior_number,animal_number) = n_target;
                area_neuron_count(area_number,type_number, ...
                    selected_behavior_number,animal_number) = n_area;
                reference_neuron_count(area_number,type_number, ...
                    selected_behavior_number,animal_number) = n_reference;

                if n_target < settings.minimum_type_neurons ...
                        || n_reference < ...
                        settings.minimum_reference_neurons
                    continue
                end

                % A deterministic seed makes resumed runs identical.
                local_seed = settings.random_seed ...
                    + 1000000*area_number ...
                    + 1000*source_behavior ...
                    + 100*type_number + animal_number;
                rng(local_seed, 'twister');
                animal_activation_index(area_number,type_number, ...
                    selected_behavior_number,animal_number,:) = ...
                    animal_local_activation_index( ...
                    target_response, area_response, settings);
            end
        end

        completed_behavior(area_number,selected_behavior_number) = true;
        saved_behavior_indices = behavior_indices; %#ok<NASGU>
        saved_animal_ids = animal_ids; %#ok<NASGU>
        saved_settings = settings; %#ok<NASGU>
        saved_input_signature = input_signature; %#ok<NASGU>
        save(calculation_file, 'animal_activation_index', ...
            'type_neuron_count', 'area_neuron_count', ...
            'reference_neuron_count', 'completed_behavior', ...
            'saved_behavior_indices', 'saved_animal_ids', ...
            'saved_settings', 'saved_input_signature', ...
            'selected_area_names', ...
            'selected_behavior_names', 'cell_type_names', ...
            'activation_time', '-v7.3');
    end
end

%% Equal-animal summaries
group_behavior_activation = nan(n_areas,n_types, ...
    n_selected_behaviors,n_time);
n_animals_by_behavior = zeros(n_areas,n_types,n_selected_behaviors);

for area_number = 1:n_areas
    for type_number = 1:n_types
        for behavior_number = 1:n_selected_behaviors
            animal_matrix = reshape( ...
                animal_activation_index(area_number,type_number, ...
                behavior_number,:,:), n_animals, n_time);
            group_behavior_activation(area_number,type_number, ...
                behavior_number,:) = ...
                mean(animal_matrix,1,'omitmissing');
            n_animals_by_behavior(area_number,type_number, ...
                behavior_number) = ...
                sum(any(isfinite(animal_matrix),2));
        end
    end
end

% Reproduce the paper's second-stage summary after animal balancing:
% each row is now one equal-animal behavior activation index. Means and
% t-test confidence intervals are calculated across behavior rows.
paper_behavior_mean = nan(n_areas,n_types,2,n_time);
paper_behavior_ci95 = nan(n_areas,n_types,2,2,n_time);
for area_number = 1:n_areas
    for type_number = 1:n_types
        play_matrix = reshape(group_behavior_activation( ...
            area_number,type_number,play_positions,:), ...
            n_play_behaviors,n_time);
        nonplay_matrix = reshape(group_behavior_activation( ...
            area_number,type_number,nonplay_positions,:), ...
            n_nonplay_behaviors,n_time);

        paper_behavior_mean(area_number,type_number,1,:) = ...
            mean(play_matrix,1,'omitmissing');
        paper_behavior_mean(area_number,type_number,2,:) = ...
            mean(nonplay_matrix,1,'omitmissing');
        [~,~,play_ci] = ttest(play_matrix);
        [~,~,nonplay_ci] = ttest(nonplay_matrix);
        paper_behavior_ci95(area_number,type_number,1,:,:) = ...
            reshape(play_ci,1,1,1,2,n_time);
        paper_behavior_ci95(area_number,type_number,2,:,:) = ...
            reshape(nonplay_ci,1,1,1,2,n_time);
    end
end

% Average behaviors inside each animal first. Animals are averaged only
% after this step, so neither neuron count nor behavior count weights them.
animal_condition_activation = nan(n_areas,n_types,2,n_animals,n_time);
animal_condition_behavior_count = zeros(n_areas,n_types,2,n_animals);
for area_number = 1:n_areas
    for type_number = 1:n_types
        for animal_number = 1:n_animals
            play_matrix = reshape(animal_activation_index( ...
                area_number,type_number,play_positions, ...
                animal_number,:), n_play_behaviors, n_time);
            nonplay_matrix = reshape(animal_activation_index( ...
                area_number,type_number,nonplay_positions, ...
                animal_number,:), n_nonplay_behaviors, n_time);

            n_play_valid = sum(all(isfinite(play_matrix),2));
            n_nonplay_valid = sum(all(isfinite(nonplay_matrix),2));
            animal_condition_behavior_count(area_number,type_number, ...
                1,animal_number) = n_play_valid;
            animal_condition_behavior_count(area_number,type_number, ...
                2,animal_number) = n_nonplay_valid;

            if n_play_valid >= settings.minimum_behaviors_per_condition
                animal_condition_activation(area_number,type_number, ...
                    1,animal_number,:) = ...
                    mean(play_matrix,1,'omitmissing');
            end
            if n_nonplay_valid >= settings.minimum_behaviors_per_condition
                animal_condition_activation(area_number,type_number, ...
                    2,animal_number,:) = ...
                    mean(nonplay_matrix,1,'omitmissing');
            end
        end

        % Play and non-play summaries use the same animals. If either
        % condition is unavailable, remove that animal from both.
        for animal_number = 1:n_animals
            play_curve = reshape(animal_condition_activation( ...
                area_number,type_number,1,animal_number,:),1,n_time);
            nonplay_curve = reshape(animal_condition_activation( ...
                area_number,type_number,2,animal_number,:),1,n_time);
            if ~all(isfinite(play_curve)) || ~all(isfinite(nonplay_curve))
                animal_condition_activation(area_number,type_number, ...
                    :,animal_number,:) = NaN;
            end
        end
    end
end

group_condition_activation = nan(n_areas,n_types,2,n_time);
n_animals_by_condition = zeros(n_areas,n_types,2);
for area_number = 1:n_areas
    for type_number = 1:n_types
        for condition_number = 1:2
            animal_matrix = reshape(animal_condition_activation( ...
                area_number,type_number,condition_number,:,:), ...
                n_animals,n_time);
            group_condition_activation(area_number,type_number, ...
                condition_number,:) = ...
                mean(animal_matrix,1,'omitmissing');
            n_animals_by_condition(area_number,type_number, ...
                condition_number) = ...
                sum(any(isfinite(animal_matrix),2));
        end
    end
end

%% Leave-one-animal-out condition curves
leave_one_out_activation = nan(n_areas,n_types,2,n_animals,n_time);
for area_number = 1:n_areas
    for type_number = 1:n_types
        contributor_matrix = reshape(animal_condition_activation( ...
            area_number,type_number,1,:,:),n_animals,n_time);
        contributors = find(all(isfinite(contributor_matrix),2));
        if numel(contributors) < 2
            continue
        end
        for excluded_animal = contributors(:)'
            retained = setdiff(contributors,excluded_animal);
            for condition_number = 1:2
                retained_matrix = reshape(animal_condition_activation( ...
                    area_number,type_number,condition_number, ...
                    retained,:), numel(retained),n_time);
                leave_one_out_activation(area_number,type_number, ...
                    condition_number,excluded_animal,:) = ...
                    mean(retained_matrix,1,'omitmissing');
            end
        end
    end
end

%% Original pooled-neuron activation index, using the same behaviors
pooled_data = load(pooled_checkpoint_file, 'all_psth_cell', ...
    'area_list', 'psth_list', 'warped_time');
if ~isfield(pooled_data, 'all_psth_cell')
    error('Pooled checkpoint does not contain all_psth_cell.');
end
pooled_area_names = string(pooled_data.area_list(:));
pooled_behavior_names = activation_behavior_names(pooled_data, ...
    numel(pooled_data.psth_list));
if ~isfield(pooled_data, 'warped_time') ...
        || numel(pooled_data.warped_time) ~= n_time ...
        || any(abs(pooled_data.warped_time(:)' - activation_time) > ...
        10*eps(max(1,max(abs(activation_time)))))
    error('Pooled and animal-level checkpoint time vectors differ.');
end
[behavior_found, pooled_behavior_indices] = ismember( ...
    selected_behavior_names, pooled_behavior_names);
if ~all(behavior_found)
    error('Pooled checkpoint is missing selected behavior labels: %s', ...
        strjoin(selected_behavior_names(~behavior_found), ', '));
end
pooled_condition_activation = nan(n_areas,n_types,2,n_time);
for area_number = 1:n_areas
    pooled_area = find(pooled_area_names == ...
        selected_area_names(area_number));
    if numel(pooled_area) ~= 1
        error('Pooled checkpoint area mismatch for %s.', ...
            selected_area_names(area_number));
    end
    for type_number = 1:n_types
        pooled_matrix = pooled_data.all_psth_cell{ ...
            type_number,pooled_area};
        if size(pooled_matrix,1) ~= numel(pooled_behavior_names) ...
                || size(pooled_matrix,2) ~= n_time
            error('Pooled checkpoint behavior/time dimensions mismatch.');
        end
        pooled_condition_activation(area_number,type_number,1,:) = ...
            mean(pooled_matrix( ...
            pooled_behavior_indices(play_positions),:),1, ...
            'omitmissing');
        pooled_condition_activation(area_number,type_number,2,:) = ...
            mean(pooled_matrix( ...
            pooled_behavior_indices(nonplay_positions),:),1, ...
            'omitmissing');
    end
end

%% Descriptive tables
coverage_table = build_coverage_table(type_neuron_count, ...
    area_neuron_count, reference_neuron_count, selected_area_names, ...
    cell_type_names, selected_behavior_names, animal_ids);
writetable(coverage_table, fullfile(output_folder, ...
    'animal_neuron_coverage.csv'));

behavior_table = build_behavior_table(animal_activation_index, ...
    activation_time, response_index, selected_area_names, ...
    cell_type_names, selected_behavior_names, animal_ids);
writetable(behavior_table, fullfile(output_folder, ...
    'animal_activation_index_by_behavior.csv'));

[condition_table, consistency_table] = build_condition_tables( ...
    animal_condition_activation, animal_condition_behavior_count, ...
    activation_time, response_index, selected_area_names, ...
    cell_type_names, animal_ids);
writetable(condition_table, fullfile(output_folder, ...
    'animal_play_nonplay_activation_index.csv'));
writetable(consistency_table, fullfile(output_folder, ...
    'play_nonplay_directional_consistency.csv'));

%% Figure 1: animal-balanced activation index per behavior
figure_handle = figure('Visible',settings.figure_visible, ...
    'Color','w', 'Units','normalized', ...
    'OuterPosition',[0 0 1 1]);
tiledlayout(n_areas,n_types, 'TileSpacing','compact', ...
    'Padding','compact');
for area_number = 1:n_areas
    for type_number = 1:n_types
        nexttile
        matrix = reshape(group_behavior_activation( ...
            area_number,type_number,:,:), ...
            n_selected_behaviors,n_time);
        imagesc(activation_time, 1:n_selected_behaviors, matrix)
        axis xy tight
        clim(settings.activation_color_limits)
        hold on
        xline(0,'k','LineWidth',1);
        xline(5,'k','LineWidth',1);
        yline(n_play_behaviors+0.5,'w','LineWidth',2);
        if type_number == 1
            yticks(1:n_selected_behaviors)
            yticklabels(selected_behavior_names)
            ylabel(selected_area_names(area_number))
        else
            yticks([])
        end
        if area_number == 1
            title(cell_type_names(type_number))
        end
        if area_number == n_areas
            xlabel('Warped time (s)')
        end
    end
end
colormap(parula)
colorbar
sgtitle(['Equal-animal activation index per behavior; ' ...
    'play above white line'])
print_figure_both(figure_handle, output_folder, ...
    'animal_balanced_activation_by_behavior', settings);

%% Figure 2: number of animals contributing to each behavior
figure_handle = figure('Visible',settings.figure_visible, ...
    'Color','w', 'Units','normalized', ...
    'OuterPosition',[0 0 1 1]);
tiledlayout(n_areas,n_types, 'TileSpacing','compact', ...
    'Padding','compact');
for area_number = 1:n_areas
    for type_number = 1:n_types
        nexttile
        counts = reshape(n_animals_by_behavior( ...
            area_number,type_number,:), n_selected_behaviors,1);
        imagesc(counts)
        axis xy tight
        clim([0 max(1,n_animals)])
        text(ones(n_selected_behaviors,1), ...
            (1:n_selected_behaviors)', string(counts), ...
            'HorizontalAlignment','center', 'Color','w', ...
            'FontWeight','bold')
        yline(n_play_behaviors+0.5,'w','LineWidth',2);
        xticks([])
        if type_number == 1
            yticks(1:n_selected_behaviors)
            yticklabels(selected_behavior_names)
            ylabel(selected_area_names(area_number))
        else
            yticks([])
        end
        if area_number == 1
            title(cell_type_names(type_number))
        end
    end
end
colormap(parula)
colorbar
sgtitle('Number of biological animals contributing to each estimate')
print_figure_both(figure_handle, output_folder, ...
    'animal_coverage_by_behavior', settings);

%% Figure 3: play and non-play animal-level curves
condition_names = ["Play", "Non-play"];
condition_colors = [0.85 0.10 0.10; 0.05 0.05 0.05];
individual_colors = [1.00 0.70 0.70; 0.65 0.65 0.65];

figure_handle = figure('Visible',settings.figure_visible, ...
    'Color','w', 'Units','normalized', ...
    'OuterPosition',[0 0 1 1]);
tiledlayout(n_areas,n_types, 'TileSpacing','compact', ...
    'Padding','compact');
for area_number = 1:n_areas
    for type_number = 1:n_types
        nexttile
        hold on
        legend_handles = gobjects(2,1);
        for condition_number = 1:2
            animal_matrix = reshape(animal_condition_activation( ...
                area_number,type_number,condition_number,:,:), ...
                n_animals,n_time);
            [lower_quartile, upper_quartile] = ...
                time_percentiles(animal_matrix, [25 75]);
            valid_time = isfinite(lower_quartile) ...
                & isfinite(upper_quartile);
            if any(valid_time)
                fill([activation_time(valid_time), ...
                    fliplr(activation_time(valid_time))], ...
                    [lower_quartile(valid_time), ...
                    fliplr(upper_quartile(valid_time))], ...
                    condition_colors(condition_number,:), ...
                    'FaceAlpha',0.15, 'EdgeColor','none');
            end
            for animal_number = 1:n_animals
                plot(activation_time, animal_matrix(animal_number,:), ...
                    'Color',individual_colors(condition_number,:), ...
                    'LineWidth',0.75);
            end
            group_curve = mean(animal_matrix,1,'omitmissing');
            legend_handles(condition_number) = ...
                plot(activation_time, group_curve, ...
                'Color',condition_colors(condition_number,:), ...
                'LineWidth',2.5);
        end
        xline(0,'k:');
        xline(5,'k:');
        yline(0,'k:');
        xlim([activation_time(1), activation_time(end)])
        if type_number == 1
            ylabel({selected_area_names(area_number), ...
                'Activation index'})
        end
        if area_number == 1
            title(cell_type_names(type_number))
        end
        if area_number == n_areas
            xlabel('Warped time (s)')
        end
        if area_number == 1 && type_number == 1
            legend(legend_handles, condition_names, ...
                'Location','best')
        end
    end
end
sgtitle(['Play and non-play: individual animals, equal-animal mean, ' ...
    'and descriptive IQR'])
print_figure_both(figure_handle, output_folder, ...
    'animal_play_nonplay_time_courses', settings);

%% Figure 3b: paper-matched mean and 95% CI across behaviors
%
% Each behavior curve is animal-balanced before entering this summary.
% The confidence interval is across behaviors, matching the paper figure;
% it is not an animal-level inferential confidence interval.
figure_handle = figure('Visible',settings.figure_visible, ...
    'Color','w', 'Units','normalized', ...
    'OuterPosition',[0 0 1 1]);
tiledlayout(n_areas,n_types, 'TileSpacing','compact', ...
    'Padding','compact');
for area_number = 1:n_areas
    for type_number = 1:n_types
        nexttile
        hold on
        summary_handles = gobjects(2,1);

        % Draw non-play first, as in the paper.
        for condition_number = [2 1]
            mean_curve = reshape(paper_behavior_mean( ...
                area_number,type_number,condition_number,:),1,n_time);
            lower_ci = reshape(paper_behavior_ci95( ...
                area_number,type_number,condition_number,1,:),1,n_time);
            upper_ci = reshape(paper_behavior_ci95( ...
                area_number,type_number,condition_number,2,:),1,n_time);
            valid_time = isfinite(lower_ci) & isfinite(upper_ci);
            if any(valid_time)
                fill([activation_time(valid_time), ...
                    fliplr(activation_time(valid_time))], ...
                    [lower_ci(valid_time), ...
                    fliplr(upper_ci(valid_time))], ...
                    condition_colors(condition_number,:), ...
                    'FaceAlpha',0.20, 'EdgeColor','none');
            end
            summary_handles(condition_number) = ...
                plot(activation_time,mean_curve, ...
                'Color',condition_colors(condition_number,:), ...
                'LineWidth',2);
        end

        xline(0,'k:');
        xline(5,'k:');
        yline(0,'k:');
        xlim([activation_time(1), activation_time(end)])
        ylim(settings.paper_summary_y_limits)
        if type_number == 1
            ylabel({selected_area_names(area_number), ...
                'Activation index'})
        end
        if area_number == 1
            title(cell_type_names(type_number))
        end
        if area_number == n_areas
            xlabel('Warped time (s)')
        end
        if area_number == 1 && type_number == 1
            legend(summary_handles,condition_names, ...
                'Location','best')
        end
    end
end
sgtitle(['Paper-matched summary: mean and 95% CI across ' ...
    'animal-balanced behavior indices'])
print_figure_both(figure_handle, output_folder, ...
    'paper_matched_animal_balanced_mean_and_ci', settings);

%% Figure 4: animal-level 0-5 s play versus non-play values
figure_handle = figure('Visible',settings.figure_visible, ...
    'Color','w', 'Units','normalized', ...
    'OuterPosition',[0 0 1 1]);
tiledlayout(n_areas,n_types, 'TileSpacing','compact', ...
    'Padding','compact');
for area_number = 1:n_areas
    for type_number = 1:n_types
        nexttile
        hold on
        animal_values = nan(n_animals,2);
        for condition_number = 1:2
            animal_matrix = reshape(animal_condition_activation( ...
                area_number,type_number,condition_number,:,:), ...
                n_animals,n_time);
            animal_values(:,condition_number) = ...
                mean(animal_matrix(:,response_index),2,'omitmissing');
        end
        for animal_number = 1:n_animals
            plot([1 2], animal_values(animal_number,:), '-', ...
                'Color',[0.65 0.65 0.65], 'LineWidth',1);
            plot(1,animal_values(animal_number,1),'o', ...
                'MarkerFaceColor',condition_colors(1,:), ...
                'MarkerEdgeColor','none', 'MarkerSize',6);
            plot(2,animal_values(animal_number,2),'o', ...
                'MarkerFaceColor',condition_colors(2,:), ...
                'MarkerEdgeColor','none', 'MarkerSize',6);
        end
        group_values = mean(animal_values,1,'omitmissing');
        plot([0.85 1.15],group_values(1)*[1 1], ...
            'Color',condition_colors(1,:), 'LineWidth',3);
        plot([1.85 2.15],group_values(2)*[1 1], ...
            'Color',condition_colors(2,:), 'LineWidth',3);
        yline(0,'k:');
        xlim([0.5 2.5])
        xticks([1 2])
        xticklabels(condition_names)
        if type_number == 1
            ylabel({selected_area_names(area_number), ...
                'Mean activation index, 0-5 s'})
        end
        if area_number == 1
            title(cell_type_names(type_number))
        end
    end
end
sgtitle('Animal-level descriptive play versus non-play comparison')
print_figure_both(figure_handle, output_folder, ...
    'animal_play_nonplay_bout_average', settings);

%% Figure 5: leave-one-animal-out robustness
figure_handle = figure('Visible',settings.figure_visible, ...
    'Color','w', 'Units','normalized', ...
    'OuterPosition',[0 0 1 1]);
tiledlayout(n_areas,n_types, 'TileSpacing','compact', ...
    'Padding','compact');
for area_number = 1:n_areas
    for type_number = 1:n_types
        nexttile
        hold on
        for condition_number = 1:2
            for excluded_animal = 1:n_animals
                curve = reshape(leave_one_out_activation( ...
                    area_number,type_number,condition_number, ...
                    excluded_animal,:),1,n_time);
                plot(activation_time, curve, ...
                    'Color',individual_colors(condition_number,:), ...
                    'LineWidth',0.75);
            end
            full_curve = reshape(group_condition_activation( ...
                area_number,type_number,condition_number,:),1,n_time);
            plot(activation_time, full_curve, ...
                'Color',condition_colors(condition_number,:), ...
                'LineWidth',2.5);
        end
        xline(0,'k:');
        xline(5,'k:');
        yline(0,'k:');
        xlim([activation_time(1), activation_time(end)])
        if type_number == 1
            ylabel({selected_area_names(area_number), ...
                'Activation index'})
        end
        if area_number == 1
            title(cell_type_names(type_number))
        end
        if area_number == n_areas
            xlabel('Warped time (s)')
        end
    end
end
sgtitle(['Leave-one-animal-out curves (light) and full ' ...
    'equal-animal summaries (dark)'])
print_figure_both(figure_handle, output_folder, ...
    'leave_one_animal_out_play_nonplay', settings);

%% Figure 6: original neuron-pooled versus animal-balanced summaries
figure_handle = figure('Visible',settings.figure_visible, ...
    'Color','w', 'Units','normalized', ...
    'OuterPosition',[0 0 1 1]);
tiledlayout(n_areas,n_types, 'TileSpacing','compact', ...
    'Padding','compact');
for area_number = 1:n_areas
    for type_number = 1:n_types
        nexttile
        hold on
        comparison_handles = gobjects(4,1);
        for condition_number = 1:2
            animal_curve = reshape(group_condition_activation( ...
                area_number,type_number,condition_number,:),1,n_time);
            pooled_curve = reshape(pooled_condition_activation( ...
                area_number,type_number,condition_number,:),1,n_time);
            pooled_handle = plot(activation_time, pooled_curve, '--', ...
                'Color',condition_colors(condition_number,:), ...
                'LineWidth',1.5);
            animal_handle = plot(activation_time, animal_curve, '-', ...
                'Color',condition_colors(condition_number,:), ...
                'LineWidth',2.5);
            comparison_handles(2*condition_number-1) = pooled_handle;
            comparison_handles(2*condition_number) = animal_handle;
        end
        xline(0,'k:');
        xline(5,'k:');
        yline(0,'k:');
        xlim([activation_time(1), activation_time(end)])
        if type_number == 1
            ylabel({selected_area_names(area_number), ...
                'Activation index'})
        end
        if area_number == 1
            title(cell_type_names(type_number))
        end
        if area_number == n_areas
            xlabel('Warped time (s)')
        end
        if area_number == 1 && type_number == 1
            legend(comparison_handles, ...
                {'Pooled play','Animal-balanced play', ...
                'Pooled non-play','Animal-balanced non-play'}, ...
                'Location','best')
        end
    end
end
sgtitle(['Dashed: original neuron-pooled estimate; solid: ' ...
    'equal-animal estimate'])
print_figure_both(figure_handle, output_folder, ...
    'pooled_vs_animal_balanced_play_nonplay', settings);

%% Save analysis bundle
save(fullfile(output_folder, ...
    'animal_balanced_activation_index_results.mat'), ...
    'settings', 'selected_area_names', 'cell_type_names', ...
    'animal_ids', 'session_to_animal', 'behavior_indices', ...
    'selected_behavior_names', 'play_behavior_indices', ...
    'nonplay_behavior_indices', 'activation_time', ...
    'animal_activation_index', 'type_neuron_count', ...
    'area_neuron_count', 'reference_neuron_count', ...
    'group_behavior_activation', 'n_animals_by_behavior', ...
    'paper_behavior_mean', 'paper_behavior_ci95', ...
    'animal_condition_activation', ...
    'animal_condition_behavior_count', ...
    'group_condition_activation', 'n_animals_by_condition', ...
    'leave_one_out_activation', 'pooled_condition_activation', ...
    'coverage_table', 'behavior_table', 'condition_table', ...
    'consistency_table', '-v7.3');

write_analysis_readme(output_folder, settings, selected_area_names, ...
    play_behavior_indices, nonplay_behavior_indices, ...
    selected_behavior_names, animal_ids);

fprintf('\nAnimal-balanced descriptive analysis complete.\n');
fprintf('Outputs: %s\n', output_folder);

%% Local functions
function names = activation_behavior_names(checkpoint, n_behaviors)
if isfield(checkpoint, 'labels') ...
        && numel(checkpoint.labels) == n_behaviors
    names = string(checkpoint.labels(:));
else
    names = string(checkpoint.psth_list(:));
    names = erase(names, "PsthWarped");
    names(names == "") = "PlayBout";
end
names = matlab.lang.makeUniqueStrings( ...
    matlab.lang.makeValidName(cellstr(names)));
names = string(names(:));
end

function animal_id = biological_animal_id(session_names)
% Sessions/probes from one focal animal receive one biological ID.
session_names = string(session_names(:));
animal_id = strings(size(session_names));
for ii = 1:numel(session_names)
    session_name = strtrim(session_names(ii));
    parts = strsplit(char(session_name));
    if numel(parts) < 3
        animal_id(ii) = session_name;
        continue
    end
    animal_batch = parts{1};
    repeated_animal = parts{3};
    if numel(animal_batch) < 2
        animal_id(ii) = session_name;
    elseif strcmp(repeated_animal, 'Single2')
        animal_id(ii) = string( ...
            ['Batch', animal_batch(2), repeated_animal]);
    elseif numel(animal_batch) >= 4
        animal_id(ii) = string( ...
            ['Batch', animal_batch(2), repeated_animal, ...
            animal_batch(4)]);
    else
        animal_id(ii) = session_name;
    end
end
end

function [responses, type_index, sessions] = activation_block( ...
        all_activation_order, behavior_number, area_number, n_time)
responses = all_activation_order{behavior_number,area_number,4};
type_index = all_activation_order{behavior_number,area_number,1};
sessions = string(all_activation_order{behavior_number,area_number,3});
type_index = type_index(:);
sessions = sessions(:);
if size(responses,2) ~= n_time ...
        || size(responses,1) ~= numel(type_index) ...
        || size(responses,1) ~= numel(sessions)
    error('PSTH, type, session, or time dimensions are not aligned.');
end
end

function compatible = compatible_calculation_settings(saved, current)
fields = {'n_permutations','minimum_type_neurons', ...
    'minimum_reference_neurons','response_window_seconds', ...
    'random_seed'};
compatible = true;
for ii = 1:numel(fields)
    field_name = fields{ii};
    if ~isfield(saved, field_name) ...
            || ~isequal(saved.(field_name), current.(field_name))
        compatible = false;
        return
    end
end
end

function activation_index = animal_local_activation_index( ...
        target_response, area_response, settings)
n_target = size(target_response,1);
n_area = size(area_response,1);
n_time = size(area_response,2);
null_median = nan(settings.n_permutations,n_time);

for permutation_number = 1:settings.n_permutations
    selected_rows = randperm(n_area,n_target);
    null_median(permutation_number,:) = ...
        median(area_response(selected_rows,:),1,'omitmissing');
end
observed_median = median(target_response,1,'omitmissing');

activation_index = nan(1,n_time);
for time_number = 1:n_time
    valid_null = isfinite(null_median(:,time_number));
    n_valid_null = sum(valid_null);
    if n_valid_null == 0 || ~isfinite(observed_median(time_number))
        continue
    end
    % Add one pseudo-count in each tail to keep the logit finite.
    probability = (1 + sum(null_median(valid_null,time_number) ...
        < observed_median(time_number))) / (n_valid_null + 2);
    activation_index(time_number) = ...
        log10(probability/(1-probability));
end
end

function [lower, upper] = time_percentiles(matrix, percentiles)
n_time = size(matrix,2);
lower = nan(1,n_time);
upper = nan(1,n_time);
for time_number = 1:n_time
    values = matrix(:,time_number);
    values = values(isfinite(values));
    if isempty(values)
        continue
    end
    limits = prctile(values,percentiles);
    lower(time_number) = limits(1);
    upper(time_number) = limits(2);
end
end

function table_out = build_coverage_table(type_count, area_count, ...
        reference_count, area_names, type_names, behavior_names, animal_ids)
out_area = strings(0,1);
out_type = strings(0,1);
out_behavior = strings(0,1);
out_animal = strings(0,1);
out_type_n = zeros(0,1);
out_reference_n = zeros(0,1);
out_area_n = zeros(0,1);
for aa = 1:numel(area_names)
    for tt = 1:numel(type_names)
        for bb = 1:numel(behavior_names)
            for nn = 1:numel(animal_ids)
                out_area(end+1,1) = area_names(aa); %#ok<AGROW>
                out_type(end+1,1) = type_names(tt); %#ok<AGROW>
                out_behavior(end+1,1) = behavior_names(bb); %#ok<AGROW>
                out_animal(end+1,1) = animal_ids(nn); %#ok<AGROW>
                out_type_n(end+1,1) = ...
                    type_count(aa,tt,bb,nn); %#ok<AGROW>
                out_reference_n(end+1,1) = ...
                    reference_count(aa,tt,bb,nn); %#ok<AGROW>
                out_area_n(end+1,1) = ...
                    area_count(aa,tt,bb,nn); %#ok<AGROW>
            end
        end
    end
end
table_out = table(out_area,out_type,out_behavior,out_animal, ...
    out_type_n,out_reference_n,out_area_n, ...
    'VariableNames', {'Area','CellType','Behavior','AnimalID', ...
    'NTypeNeurons','NReferenceNeurons','NAreaNeurons'});
end

function table_out = build_behavior_table(animal_index, time, ...
        response_index, area_names, type_names, behavior_names, animal_ids)
out_area = strings(0,1);
out_type = strings(0,1);
out_behavior = strings(0,1);
out_animal = strings(0,1);
out_value = zeros(0,1);
n_time = numel(time);
for aa = 1:numel(area_names)
    for tt = 1:numel(type_names)
        for bb = 1:numel(behavior_names)
            for nn = 1:numel(animal_ids)
                curve = reshape(animal_index(aa,tt,bb,nn,:),1,n_time);
                out_area(end+1,1) = area_names(aa); %#ok<AGROW>
                out_type(end+1,1) = type_names(tt); %#ok<AGROW>
                out_behavior(end+1,1) = behavior_names(bb); %#ok<AGROW>
                out_animal(end+1,1) = animal_ids(nn); %#ok<AGROW>
                out_value(end+1,1) = ...
                    mean(curve(response_index),'omitmissing'); %#ok<AGROW>
            end
        end
    end
end
table_out = table(out_area,out_type,out_behavior,out_animal,out_value, ...
    'VariableNames', {'Area','CellType','Behavior','AnimalID', ...
    'MeanActivationIndex0to5s'});
end

function [condition_table, consistency_table] = build_condition_tables( ...
        condition_index, behavior_count, time, response_index, ...
        area_names, type_names, animal_ids)
condition_names = ["Play"; "Non-play"];
n_time = numel(time);
out_area = strings(0,1);
out_type = strings(0,1);
out_condition = strings(0,1);
out_animal = strings(0,1);
out_n_behaviors = zeros(0,1);
out_value = zeros(0,1);

cons_area = strings(0,1);
cons_type = strings(0,1);
cons_n = zeros(0,1);
cons_play_higher = zeros(0,1);
cons_play_lower = zeros(0,1);
cons_equal = zeros(0,1);
cons_mean_difference = zeros(0,1);

for aa = 1:numel(area_names)
    for tt = 1:numel(type_names)
        animal_values = nan(numel(animal_ids),2);
        for cc = 1:2
            for nn = 1:numel(animal_ids)
                curve = reshape(condition_index(aa,tt,cc,nn,:),1,n_time);
                value = mean(curve(response_index),'omitmissing');
                animal_values(nn,cc) = value;
                out_area(end+1,1) = area_names(aa); %#ok<AGROW>
                out_type(end+1,1) = type_names(tt); %#ok<AGROW>
                out_condition(end+1,1) = condition_names(cc); %#ok<AGROW>
                out_animal(end+1,1) = animal_ids(nn); %#ok<AGROW>
                out_n_behaviors(end+1,1) = ...
                    behavior_count(aa,tt,cc,nn); %#ok<AGROW>
                out_value(end+1,1) = value; %#ok<AGROW>
            end
        end
        difference = animal_values(:,1)-animal_values(:,2);
        valid = isfinite(difference);
        cons_area(end+1,1) = area_names(aa); %#ok<AGROW>
        cons_type(end+1,1) = type_names(tt); %#ok<AGROW>
        cons_n(end+1,1) = sum(valid); %#ok<AGROW>
        cons_play_higher(end+1,1) = sum(difference(valid)>0); %#ok<AGROW>
        cons_play_lower(end+1,1) = sum(difference(valid)<0); %#ok<AGROW>
        cons_equal(end+1,1) = sum(difference(valid)==0); %#ok<AGROW>
        cons_mean_difference(end+1,1) = ...
            mean(difference(valid),'omitmissing'); %#ok<AGROW>
    end
end

condition_table = table(out_area,out_type,out_condition,out_animal, ...
    out_n_behaviors,out_value, 'VariableNames', ...
    {'Area','CellType','Condition','AnimalID','NBehaviors', ...
    'MeanActivationIndex0to5s'});
consistency_table = table(cons_area,cons_type,cons_n, ...
    cons_play_higher,cons_play_lower,cons_equal,cons_mean_difference, ...
    'VariableNames', {'Area','CellType','NAnimals', ...
    'NPlayHigher','NPlayLower','NEqual','MeanPlayMinusNonPlay'});
end

function print_figure_both(figure_handle, output_folder, stem, settings)
drawnow
jpg_file = fullfile(output_folder, [stem, '.jpg']);
svg_file = fullfile(output_folder, [stem, '.svg']);
print(figure_handle, jpg_file, '-djpeg', ...
    ['-r',num2str(settings.jpeg_resolution)]);
print(figure_handle, svg_file, '-dsvg', '-vector');
end

function write_analysis_readme(output_folder, settings, area_names, ...
        play_indices, nonplay_indices, behavior_names, animal_ids)
file_name = fullfile(output_folder, 'README_animal_level_analysis.txt');
file_id = fopen(file_name,'w');
if file_id < 0
    warning('Could not write %s.', file_name);
    return
end
cleanup = onCleanup(@() fclose(file_id));
fprintf(file_id, 'ANIMAL-BALANCED ACTIVATION-INDEX ANALYSIS\n\n');
fprintf(file_id, ['This is a descriptive analysis. It contains no ' ...
    'inferential p-values.\n']);
fprintf(file_id, ['Each activation index uses only neurons from the same ' ...
    'animal and area.\n']);
fprintf(file_id, ['Behaviors are averaged within animal before animals ' ...
    'are given equal weight.\n\n']);
fprintf(file_id, ['The paper-matched mean/CI figure first equal-weights ' ...
    'animals within each behavior, then computes its 95%% CI across ' ...
    'behaviors. That interval is not animal-level inference.\n\n']);
fprintf(file_id, 'Areas: %s\n', strjoin(area_names, ', '));
fprintf(file_id, 'Animals (%d): %s\n', numel(animal_ids), ...
    strjoin(animal_ids, ', '));
fprintf(file_id, 'Permutations per animal-local index: %d\n', ...
    settings.n_permutations);
fprintf(file_id, 'Minimum target neurons: %d\n', ...
    settings.minimum_type_neurons);
fprintf(file_id, 'Minimum non-target reference neurons: %d\n', ...
    settings.minimum_reference_neurons);
fprintf(file_id, 'Minimum behaviors per condition: %d\n\n', ...
    settings.minimum_behaviors_per_condition);
fprintf(file_id, 'Play indices: %s\n', mat2str(play_indices));
fprintf(file_id, 'Non-play indices: %s\n', mat2str(nonplay_indices));
fprintf(file_id, 'Displayed behavior order:\n');
for ii = 1:numel(behavior_names)
    fprintf(file_id, '  %d. %s\n', ii, behavior_names(ii));
end
fprintf(file_id, '\nFigures are saved as JPG and vector SVG with print.\n');
fprintf(file_id, ['The original neuron-pooled curves are included only ' ...
    'as descriptive comparisons.\n']);
end
