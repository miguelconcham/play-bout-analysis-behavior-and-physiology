%% Assess whether activation-index inference can use animals as replicates
%
% This script is a feasibility audit, not the final biological analysis.
% It deliberately does not plot or export observed activation-index means,
% signs, or p-values. It answers:
%   1. How many biological animals and neurons support every comparison?
%   2. How precise is an animal-local activation index at the available n?
%   3. What false-positive rate and power would an animal-level test have
%      under simulated effects of known size?
%   4. Would leave-one- or leave-two-animal-out checks retain enough animals?
%
% Required input:
%   Data/Analysis results/activation index/activation_index_A1_A7.mat
% produced by:
%   Figure 3/Supporting codes/Estiamte_activation_index_and_other_variables.m
%
% Outputs are written to:
%   Data/Analysis results/activation index/hierarchical feasibility/
%
% Important:
%   - Sessions/probes from one focal animal are kept in one AnimalID.
%   - Coverage decisions are made before inspecting biological effect signs.
%   - The thresholds below are design choices. Change them before running,
%     then retain the settings MAT file with the results.

clearvars
close all
clc

%% 1) Paths and prespecified feasibility settings
experiment_root = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology'];
checkpoint_folder = fullfile(experiment_root, 'Data', 'Analysis results', ...
    'activation index');
checkpoint_file = fullfile(checkpoint_folder, ...
    'activation_index_A1_A7.mat');
output_folder = fullfile(checkpoint_folder, ...
    'hierarchical feasibility');

settings = struct;
settings.response_window_seconds = [0 5];
settings.minimum_neurons_per_type = 3;
settings.recommended_neurons_per_type = 8;
settings.minimum_animals = 4;
settings.candidate_neuron_counts = [3 5 8 10];
settings.bootstrap_repetitions = 250;
settings.null_repetitions = 1000;
settings.simulation_repetitions = 500;
settings.simulated_effect_sd = [0 0.25 0.5 1];
settings.alpha = 0.05;
settings.maximum_acceptable_ai95_width = 2;
settings.maximum_null_false_positive_rate = 0.075;
settings.minimum_moderate_effect_power = 0.80;
settings.random_seed = 8613;

if ~isfile(checkpoint_file)
    error(['Missing checkpoint:\n%s\nRun ' ...
        'Estiamte_activation_index_and_other_variables.m through the ' ...
        'A1-A7 checkpoint save first.'], checkpoint_file);
end
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
rng(settings.random_seed, 'twister');

required_variables = {'all_activation_order', 'area_list', 'psth_list'};
checkpoint = load(checkpoint_file);
missing_variables = required_variables( ...
    ~isfield(checkpoint, required_variables));
if ~isempty(missing_variables)
    error('Checkpoint is missing: %s', ...
        strjoin(missing_variables, ', '));
end

all_activation_order = checkpoint.all_activation_order;
area_names = string(checkpoint.area_list(:));
n_behaviors = size(all_activation_order, 1);
n_areas = size(all_activation_order, 2);
cell_type_names = ["Peak-locked"; "Trough-locked"; "Unlocked"];
n_types = numel(cell_type_names);

if n_areas ~= numel(area_names)
    error('Area labels do not match all_activation_order columns.');
end

behavior_names = activation_behavior_names(checkpoint, n_behaviors);
n_time = size(all_activation_order{1,1,4}, 2);
if isfield(checkpoint, 'warped_time') ...
        && numel(checkpoint.warped_time) == n_time
    response_time = checkpoint.warped_time(:)';
elseif isfield(checkpoint, 'time2use') ...
        && numel(checkpoint.time2use) == n_time
    response_time = checkpoint.time2use(:)';
else
    error(['No time vector in the checkpoint matches the stored PSTHs. ' ...
        'Rebuild activation_index_A1_A7.mat.']);
end
response_index = response_time >= settings.response_window_seconds(1) ...
    & response_time <= settings.response_window_seconds(2);
if ~any(response_index)
    error('The requested response window has no samples.');
end

%% 2) Validate session-to-animal grouping
all_sessions = strings(0,1);
for area_number = 1:n_areas
    sessions = string(all_activation_order{1,area_number,3});
    all_sessions = [all_sessions; sessions(:)]; %#ok<AGROW>
end
all_sessions = unique(all_sessions(strlength(all_sessions) > 0), 'stable');
session_to_animal = table;
session_to_animal.Session = all_sessions;
session_to_animal.AnimalID = biological_animal_id(all_sessions);
animal_ids = unique(session_to_animal.AnimalID, 'stable');
n_animals_total = numel(animal_ids);

writetable(session_to_animal, fullfile(output_folder, ...
    'session_to_biological_animal.csv'));

fprintf('Found %d sessions grouped into %d biological animals.\n', ...
    height(session_to_animal), n_animals_total);
disp(session_to_animal)

%% 3) Coverage audit: animal x area x type x behavior
cov_area = strings(0,1);
cov_type = strings(0,1);
cov_behavior = strings(0,1);
cov_animal = strings(0,1);
cov_n_type = zeros(0,1);
cov_n_area = zeros(0,1);
cov_n_sessions = zeros(0,1);
cov_session_ids = strings(0,1);

for area_number = 1:n_areas
    for behavior_number = 1:n_behaviors
        [responses, type_index, sessions] = activation_block( ...
            all_activation_order, behavior_number, area_number, n_time);
        animal_for_row = biological_animal_id(sessions);
        scalar_response = median(responses(:,response_index), 2, ...
            'omitmissing');
        valid_area = isfinite(scalar_response);

        for type_number = 1:n_types
            for animal_number = 1:n_animals_total
                this_animal = animal_ids(animal_number);
                animal_mask = animal_for_row == this_animal;
                type_mask = type_index == type_number;
                valid_type = animal_mask & type_mask & valid_area;
                animal_sessions = unique(sessions(animal_mask), 'stable');

                cov_area(end+1,1) = area_names(area_number); %#ok<SAGROW>
                cov_type(end+1,1) = cell_type_names(type_number); %#ok<SAGROW>
                cov_behavior(end+1,1) = ...
                    behavior_names(behavior_number); %#ok<SAGROW>
                cov_animal(end+1,1) = this_animal; %#ok<SAGROW>
                cov_n_type(end+1,1) = sum(valid_type); %#ok<SAGROW>
                cov_n_area(end+1,1) = sum(animal_mask & valid_area); %#ok<SAGROW>
                cov_n_sessions(end+1,1) = numel(animal_sessions); %#ok<SAGROW>
                cov_session_ids(end+1,1) = ...
                    strjoin(animal_sessions, ';'); %#ok<SAGROW>
            end
        end
    end
end

coverage_by_animal = table(cov_area, cov_type, cov_behavior, cov_animal, ...
    cov_n_type, cov_n_area, cov_n_sessions, cov_session_ids, ...
    'VariableNames', {'Area','CellType','Behavior','AnimalID', ...
    'NTypeValid','NAreaValid','NSessions','SessionIDs'});
coverage_by_animal.TypeFraction = coverage_by_animal.NTypeValid ./ ...
    coverage_by_animal.NAreaValid;
coverage_by_animal.TypeFraction( ...
    coverage_by_animal.NAreaValid == 0) = NaN;
writetable(coverage_by_animal, fullfile(output_folder, ...
    'coverage_by_animal.csv'));

coverage_summary = summarize_coverage(coverage_by_animal, area_names, ...
    cell_type_names, behavior_names, settings);
writetable(coverage_summary, fullfile(output_folder, ...
    'coverage_summary.csv'));

%% 4) Bootstrap precision of animal-local activation indices
%
% The observed activation-index location is intentionally not retained.
% Only SD and interval width are exported. Null subsets are sampled from
% neurons in the same animal and area.
stability_area = strings(0,1);
stability_type = strings(0,1);
stability_behavior = strings(0,1);
stability_animal = strings(0,1);
stability_mode = strings(0,1);
stability_sample_n = zeros(0,1);
stability_available_n = zeros(0,1);
stability_area_n = zeros(0,1);
stability_sd = zeros(0,1);
stability_width = zeros(0,1);
stability_stable = false(0,1);

for area_number = 1:n_areas
    fprintf('Precision audit: area %d/%d (%s)\n', area_number, ...
        n_areas, area_names(area_number));
    for behavior_number = 1:n_behaviors
        [responses, type_index, sessions] = activation_block( ...
            all_activation_order, behavior_number, area_number, n_time);
        scalar_response = median(responses(:,response_index), 2, ...
            'omitmissing');
        animal_for_row = biological_animal_id(sessions);

        for type_number = 1:n_types
            for animal_number = 1:n_animals_total
                this_animal = animal_ids(animal_number);
                animal_mask = animal_for_row == this_animal ...
                    & isfinite(scalar_response);
                target_mask = animal_mask & type_index == type_number;
                target_values = scalar_response(target_mask);
                area_values = scalar_response(animal_mask);
                n_target = numel(target_values);
                n_area_valid = numel(area_values);
                if n_target < settings.minimum_neurons_per_type ...
                        || n_area_valid < n_target
                    continue
                end

                sample_sizes = unique([n_target, ...
                    settings.candidate_neuron_counts( ...
                    settings.candidate_neuron_counts <= n_target)]);
                for sample_number = 1:numel(sample_sizes)
                    sample_n = sample_sizes(sample_number);
                    if sample_n == n_target
                        mode_name = "Observed";
                    else
                        mode_name = "Downsample";
                    end

                    [ai_sd, ai_width] = activation_index_precision( ...
                        target_values, area_values, sample_n, settings);

                    stability_area(end+1,1) = ...
                        area_names(area_number); %#ok<SAGROW>
                    stability_type(end+1,1) = ...
                        cell_type_names(type_number); %#ok<SAGROW>
                    stability_behavior(end+1,1) = ...
                        behavior_names(behavior_number); %#ok<SAGROW>
                    stability_animal(end+1,1) = this_animal; %#ok<SAGROW>
                    stability_mode(end+1,1) = mode_name; %#ok<SAGROW>
                    stability_sample_n(end+1,1) = sample_n; %#ok<SAGROW>
                    stability_available_n(end+1,1) = n_target; %#ok<SAGROW>
                    stability_area_n(end+1,1) = n_area_valid; %#ok<SAGROW>
                    stability_sd(end+1,1) = ai_sd; %#ok<SAGROW>
                    stability_width(end+1,1) = ai_width; %#ok<SAGROW>
                    stability_stable(end+1,1) = ai_width <= ...
                        settings.maximum_acceptable_ai95_width; %#ok<SAGROW>
                end
            end
        end
    end
end

stability_by_animal = table(stability_area, stability_type, ...
    stability_behavior, stability_animal, stability_mode, ...
    stability_sample_n, stability_available_n, stability_area_n, ...
    stability_sd, stability_width, stability_stable, ...
    'VariableNames', {'Area','CellType','Behavior','AnimalID','Mode', ...
    'SampleN','AvailableTypeN','AvailableAreaN','ActivationIndexSD', ...
    'ActivationIndex95Width','WithinPrecisionLimit'});
writetable(stability_by_animal, fullfile(output_folder, ...
    'activation_index_precision_by_animal.csv'));

stability_summary = summarize_stability(stability_by_animal, ...
    area_names, cell_type_names, behavior_names, settings);
writetable(stability_summary, fullfile(output_folder, ...
    'activation_index_precision_summary.csv'));

%% 5) Blinded operating-characteristic simulation
%
% Each animal's target and non-target distributions are separately centered,
% removing the observed location effect. Known standardized shifts are then
% added to target neurons. Target neurons remain inside the jointly resampled
% area pool, matching the original activation-index subset relationship. The
% animal statistic is a median target-minus-area contrast. Exact animal-level
% sign flips provide one- and two-sided tests.
sim_area = strings(0,1);
sim_type = strings(0,1);
sim_behavior = strings(0,1);
sim_effect = zeros(0,1);
sim_n_animals = zeros(0,1);
sim_min_two_sided_p = zeros(0,1);
sim_min_one_sided_p = zeros(0,1);
sim_two_sided_detection = zeros(0,1);
sim_one_sided_detection = zeros(0,1);

for area_number = 1:n_areas
    fprintf('Simulation audit: area %d/%d (%s)\n', area_number, ...
        n_areas, area_names(area_number));
    for behavior_number = 1:n_behaviors
        [responses, type_index, sessions] = activation_block( ...
            all_activation_order, behavior_number, area_number, n_time);
        scalar_response = median(responses(:,response_index), 2, ...
            'omitmissing');
        animal_for_row = biological_animal_id(sessions);

        for type_number = 1:n_types
            distributions = cell(0,2);
            included_animals = strings(0,1);
            for animal_number = 1:n_animals_total
                this_animal = animal_ids(animal_number);
                animal_mask = animal_for_row == this_animal ...
                    & isfinite(scalar_response);
                target_values = scalar_response( ...
                    animal_mask & type_index == type_number);
                reference_values = scalar_response( ...
                    animal_mask & type_index ~= type_number);
                if numel(target_values) < ...
                        settings.minimum_neurons_per_type ...
                        || numel(reference_values) < ...
                        settings.minimum_neurons_per_type
                    continue
                end
                distributions(end+1,:) = ...
                    {target_values, reference_values}; %#ok<SAGROW>
                included_animals(end+1,1) = this_animal; %#ok<SAGROW>
            end

            n_included = numel(included_animals);
            for effect_number = 1:numel(settings.simulated_effect_sd)
                effect_sd = settings.simulated_effect_sd(effect_number);
                if n_included < 2
                    two_sided_rate = NaN;
                    one_sided_rate = NaN;
                else
                    [two_sided_rate, one_sided_rate] = ...
                        simulate_animal_inference(distributions, ...
                        effect_sd, settings);
                end

                sim_area(end+1,1) = area_names(area_number); %#ok<SAGROW>
                sim_type(end+1,1) = ...
                    cell_type_names(type_number); %#ok<SAGROW>
                sim_behavior(end+1,1) = ...
                    behavior_names(behavior_number); %#ok<SAGROW>
                sim_effect(end+1,1) = effect_sd; %#ok<SAGROW>
                sim_n_animals(end+1,1) = n_included; %#ok<SAGROW>
                sim_min_two_sided_p(end+1,1) = ...
                    minimum_exact_p(n_included, "two-sided"); %#ok<SAGROW>
                sim_min_one_sided_p(end+1,1) = ...
                    minimum_exact_p(n_included, "one-sided"); %#ok<SAGROW>
                sim_two_sided_detection(end+1,1) = ...
                    two_sided_rate; %#ok<SAGROW>
                sim_one_sided_detection(end+1,1) = ...
                    one_sided_rate; %#ok<SAGROW>
            end
        end
    end
end

simulation_summary = table(sim_area, sim_type, sim_behavior, sim_effect, ...
    sim_n_animals, sim_min_two_sided_p, sim_min_one_sided_p, ...
    sim_two_sided_detection, sim_one_sided_detection, ...
    'VariableNames', {'Area','CellType','Behavior','InjectedEffectSD', ...
    'NAnimals','MinimumExactTwoSidedP','MinimumExactOneSidedP', ...
    'TwoSidedDetectionRate','OneSidedDetectionRate'});
writetable(simulation_summary, fullfile(output_folder, ...
    'animal_level_simulation_summary.csv'));

%% 6) Final design-only decision table
decision_table = coverage_summary;
decision_table.MedianObservedAI95Width = ...
    nan(height(decision_table),1);
decision_table.StableAnimalFraction = ...
    nan(height(decision_table),1);
decision_table.NullTwoSidedFalsePositiveRate = ...
    nan(height(decision_table),1);
decision_table.ModerateEffectTwoSidedPower = ...
    nan(height(decision_table),1);
decision_table.MinimumExactTwoSidedP = ...
    nan(height(decision_table),1);
decision_table.Decision = repmat("Insufficient coverage", ...
    height(decision_table),1);

for row_number = 1:height(decision_table)
    key = stability_summary.Area == decision_table.Area(row_number) ...
        & stability_summary.CellType == ...
        decision_table.CellType(row_number) ...
        & stability_summary.Behavior == ...
        decision_table.Behavior(row_number);
    if any(key)
        decision_table.MedianObservedAI95Width(row_number) = ...
            stability_summary.MedianObservedAI95Width(key);
        decision_table.StableAnimalFraction(row_number) = ...
            stability_summary.StableAnimalFraction(key);
    end

    sim_key = simulation_summary.Area == ...
        decision_table.Area(row_number) ...
        & simulation_summary.CellType == ...
        decision_table.CellType(row_number) ...
        & simulation_summary.Behavior == ...
        decision_table.Behavior(row_number);
    null_key = sim_key & simulation_summary.InjectedEffectSD == 0;
    moderate_key = sim_key & simulation_summary.InjectedEffectSD == 0.5;
    if any(null_key)
        decision_table.NullTwoSidedFalsePositiveRate(row_number) = ...
            simulation_summary.TwoSidedDetectionRate(null_key);
        decision_table.MinimumExactTwoSidedP(row_number) = ...
            simulation_summary.MinimumExactTwoSidedP(null_key);
    end
    if any(moderate_key)
        decision_table.ModerateEffectTwoSidedPower(row_number) = ...
            simulation_summary.TwoSidedDetectionRate(moderate_key);
    end

    if decision_table.NAnimalsRepresented(row_number) < ...
            settings.minimum_animals
        decision_table.Decision(row_number) = ...
            "Insufficient coverage";
    elseif ~isfinite(decision_table.MinimumExactTwoSidedP(row_number)) ...
            || decision_table.MinimumExactTwoSidedP(row_number) >= ...
            settings.alpha
        decision_table.Decision(row_number) = ...
            "Insufficient animals for exact test";
    elseif decision_table.NullTwoSidedFalsePositiveRate(row_number) > ...
            settings.maximum_null_false_positive_rate
        decision_table.Decision(row_number) = ...
            "Simulation type-I error failure";
    elseif ~isfinite( ...
            decision_table.ModerateEffectTwoSidedPower(row_number)) ...
            || decision_table.ModerateEffectTwoSidedPower(row_number) < ...
            settings.minimum_moderate_effect_power
        decision_table.Decision(row_number) = ...
            "Underpowered at 0.5 SD";
    elseif decision_table.NAnimalsRecommended(row_number) < ...
            settings.minimum_animals
        decision_table.Decision(row_number) = ...
            "Use stratified group permutation";
    else
        decision_table.Decision(row_number) = ...
            "Per-animal index feasible";
    end
end

writetable(decision_table, fullfile(output_folder, ...
    'feasibility_decision_table.csv'));

%% 7) Figures (coverage and precision only; no observed effects)
make_coverage_figures(coverage_summary, area_names, cell_type_names, ...
    behavior_names, n_animals_total, output_folder);
make_animal_count_figures(coverage_by_animal, area_names, ...
    cell_type_names, behavior_names, animal_ids, output_folder);
make_precision_figures(stability_summary, area_names, ...
    cell_type_names, behavior_names, settings, output_folder);
make_simulation_figures(simulation_summary, area_names, ...
    cell_type_names, behavior_names, settings, output_folder);
make_decision_figures(decision_table, area_names, cell_type_names, ...
    behavior_names, output_folder);

%% 8) Save reproducibility bundle and concise report
save(fullfile(output_folder, 'hierarchical_feasibility_results.mat'), ...
    'settings', 'session_to_animal', 'coverage_by_animal', ...
    'coverage_summary', 'stability_by_animal', 'stability_summary', ...
    'simulation_summary', 'decision_table', 'area_names', ...
    'cell_type_names', 'behavior_names', '-v7.3');

report_file = fullfile(output_folder, 'README_feasibility_results.txt');
write_feasibility_report(report_file, checkpoint_file, settings, ...
    session_to_animal, decision_table);

fprintf('\nFeasibility audit complete.\nResults: %s\n', output_folder);
fprintf(['Start review with README_feasibility_results.txt and ' ...
    'feasibility_decision_table.csv.\n']);

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
% Sessions of one focal animal share the Batch/recording-animal ID.
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
responses = all_activation_order{behavior_number, area_number, 4};
type_index = all_activation_order{behavior_number, area_number, 1};
sessions = string( ...
    all_activation_order{behavior_number, area_number, 3});
type_index = type_index(:);
sessions = sessions(:);
if size(responses,2) ~= n_time
    error('PSTH width differs across behavior/area blocks.');
end
if size(responses,1) ~= numel(type_index) ...
        || size(responses,1) ~= numel(sessions)
    error('PSTH, type, and session rows are not aligned.');
end
end

function summary = summarize_coverage(coverage, area_names, ...
        type_names, behavior_names, settings)
out_area = strings(0,1);
out_type = strings(0,1);
out_behavior = strings(0,1);
out_any = zeros(0,1);
out_represented = zeros(0,1);
out_recommended = zeros(0,1);
out_min = zeros(0,1);
out_median = zeros(0,1);
out_max = zeros(0,1);
out_loao = false(0,1);
out_ltao = false(0,1);
out_status = strings(0,1);

for aa = 1:numel(area_names)
    for tt = 1:numel(type_names)
        for bb = 1:numel(behavior_names)
            mask = coverage.Area == area_names(aa) ...
                & coverage.CellType == type_names(tt) ...
                & coverage.Behavior == behavior_names(bb);
            counts = coverage.NTypeValid(mask);
            positive_counts = counts(counts > 0);
            n_any = sum(counts > 0);
            n_represented = sum(counts >= ...
                settings.minimum_neurons_per_type);
            n_recommended = sum(counts >= ...
                settings.recommended_neurons_per_type);
            if isempty(positive_counts)
                min_count = 0;
                median_count = 0;
                max_count = 0;
            else
                min_count = min(positive_counts);
                median_count = median(positive_counts);
                max_count = max(positive_counts);
            end

            if n_recommended >= settings.minimum_animals
                status = "Good";
            elseif n_represented >= settings.minimum_animals
                status = "Marginal";
            else
                status = "Insufficient";
            end

            out_area(end+1,1) = area_names(aa); %#ok<AGROW>
            out_type(end+1,1) = type_names(tt); %#ok<AGROW>
            out_behavior(end+1,1) = behavior_names(bb); %#ok<AGROW>
            out_any(end+1,1) = n_any; %#ok<AGROW>
            out_represented(end+1,1) = n_represented; %#ok<AGROW>
            out_recommended(end+1,1) = n_recommended; %#ok<AGROW>
            out_min(end+1,1) = min_count; %#ok<AGROW>
            out_median(end+1,1) = median_count; %#ok<AGROW>
            out_max(end+1,1) = max_count; %#ok<AGROW>
            out_loao(end+1,1) = n_represented >= ...
                settings.minimum_animals + 1; %#ok<AGROW>
            out_ltao(end+1,1) = n_represented >= ...
                settings.minimum_animals + 2; %#ok<AGROW>
            out_status(end+1,1) = status; %#ok<AGROW>
        end
    end
end

summary = table(out_area, out_type, out_behavior, out_any, ...
    out_represented, out_recommended, out_min, out_median, out_max, ...
    out_loao, out_ltao, out_status, ...
    'VariableNames', {'Area','CellType','Behavior','NAnimalsAny', ...
    'NAnimalsRepresented','NAnimalsRecommended','MinimumPositiveN', ...
    'MedianPositiveN','MaximumN','LeaveOneAnimalOutFeasible', ...
    'LeaveTwoAnimalsOutFeasible','CoverageStatus'});
end

function [ai_sd, ai_width] = activation_index_precision( ...
        target_values, area_values, sample_n, settings)
% Bootstrap target neurons; compare each median with one animal-local null.
target_values = target_values(isfinite(target_values));
area_values = area_values(isfinite(area_values));
n_target = numel(target_values);
n_area = numel(area_values);

null_samples = nan(settings.null_repetitions, sample_n);
for null_number = 1:settings.null_repetitions
    null_index = randperm(n_area, sample_n);
    null_samples(null_number,:) = area_values(null_index);
end
null_medians = median(null_samples, 2, 'omitmissing');
target_samples = target_values(randi(n_target, ...
    settings.bootstrap_repetitions, sample_n));
target_medians = median(target_samples, 2, 'omitmissing');

activation_index = nan(settings.bootstrap_repetitions,1);
for ii = 1:settings.bootstrap_repetitions
    % Add one pseudo-count to avoid infinite logits.
    probability = (1 + sum(null_medians < target_medians(ii))) ...
        / (settings.null_repetitions + 2);
    activation_index(ii) = log10(probability / (1-probability));
end

ai_sd = std(activation_index, 'omitmissing');
limits = prctile(activation_index, [2.5 97.5]);
ai_width = limits(2) - limits(1);
end

function summary = summarize_stability(stability, area_names, ...
        type_names, behavior_names, settings)
out_area = strings(0,1);
out_type = strings(0,1);
out_behavior = strings(0,1);
out_n = zeros(0,1);
out_width = zeros(0,1);
out_fraction = zeros(0,1);

observed = stability(stability.Mode == "Observed", :);
for aa = 1:numel(area_names)
    for tt = 1:numel(type_names)
        for bb = 1:numel(behavior_names)
            mask = observed.Area == area_names(aa) ...
                & observed.CellType == type_names(tt) ...
                & observed.Behavior == behavior_names(bb);
            widths = observed.ActivationIndex95Width(mask);
            stable = widths <= ...
                settings.maximum_acceptable_ai95_width;
            out_area(end+1,1) = area_names(aa); %#ok<AGROW>
            out_type(end+1,1) = type_names(tt); %#ok<AGROW>
            out_behavior(end+1,1) = behavior_names(bb); %#ok<AGROW>
            out_n(end+1,1) = numel(widths); %#ok<AGROW>
            if isempty(widths)
                out_width(end+1,1) = NaN; %#ok<AGROW>
                out_fraction(end+1,1) = NaN; %#ok<AGROW>
            else
                out_width(end+1,1) = ...
                    median(widths, 'omitmissing'); %#ok<AGROW>
                out_fraction(end+1,1) = ...
                    mean(stable, 'omitmissing'); %#ok<AGROW>
            end
        end
    end
end

summary = table(out_area, out_type, out_behavior, out_n, out_width, ...
    out_fraction, 'VariableNames', {'Area','CellType','Behavior', ...
    'NAnimalsTested','MedianObservedAI95Width','StableAnimalFraction'});
end

function [two_sided_rate, one_sided_rate] = ...
        simulate_animal_inference(distributions, effect_sd, settings)
n_animals = size(distributions,1);
two_sided_detected = false(settings.simulation_repetitions,1);
one_sided_detected = false(settings.simulation_repetitions,1);

for rr = 1:settings.simulation_repetitions
    animal_contrast = nan(n_animals,1);
    for aa = 1:n_animals
        target = distributions{aa,1};
        reference = distributions{aa,2};
        target = target(isfinite(target));
        reference = reference(isfinite(reference));

        target_residual = target - median(target, 'omitmissing');
        reference_residual = reference - ...
            median(reference, 'omitmissing');
        centered_area = [target_residual; reference_residual];
        robust_scale = 1.4826 * median(abs(centered_area), ...
            'omitmissing');
        if ~isfinite(robust_scale) || robust_scale <= eps
            robust_scale = std(centered_area, 'omitmissing');
        end
        if ~isfinite(robust_scale) || robust_scale <= eps
            robust_scale = 1;
        end

        simulated_target = target_residual(randi(numel(target_residual), ...
            numel(target_residual), 1)) + effect_sd * robust_scale;
        simulated_reference = reference_residual(randi( ...
            numel(reference_residual), numel(reference_residual), 1));
        % Preserve the producer's subset relationship: target neurons are
        % part of the area-wide comparison pool in the same resample.
        simulated_area = [simulated_target; simulated_reference];
        animal_contrast(aa) = ...
            (median(simulated_target, 'omitmissing') - ...
            median(simulated_area, 'omitmissing')) / robust_scale;
    end

    [p_two, p_one] = exact_sign_flip_p(animal_contrast);
    two_sided_detected(rr) = p_two < settings.alpha;
    one_sided_detected(rr) = p_one < settings.alpha;
end

two_sided_rate = mean(two_sided_detected);
one_sided_rate = mean(one_sided_detected);
end

function [p_two, p_one] = exact_sign_flip_p(values)
values = values(isfinite(values));
n_values = numel(values);
if n_values == 0
    p_two = NaN;
    p_one = NaN;
    return
end

n_assignments = 2^n_values;
signs = nan(n_assignments, n_values);
for assignment = 0:n_assignments-1
    signs(assignment+1,:) = 2*bitget(assignment, 1:n_values)-1;
end
null_means = mean(signs .* values(:)', 2);
observed = mean(values);
tolerance = 10*eps(max(1, abs(observed)));
p_two = sum(abs(null_means) >= abs(observed)-tolerance) ...
    / n_assignments;
p_one = sum(null_means >= observed-tolerance) / n_assignments;
end

function value = minimum_exact_p(n_animals, tail)
if n_animals < 1
    value = NaN;
elseif tail == "two-sided"
    value = min(1, 2^(1-n_animals));
else
    value = 2^(-n_animals);
end
end

function make_coverage_figures(summary, area_names, type_names, ...
        behavior_names, n_animals, output_folder)
for aa = 1:numel(area_names)
    figure_handle = figure('Visible','off', 'Color','w', ...
        'Position',[50 50 1500 750]);
    metrics = {'NAnimalsAny','NAnimalsRepresented', ...
        'NAnimalsRecommended'};
    titles = {'Animals with >=1 neuron', ...
        'Animals meeting minimum neuron count', ...
        'Animals meeting recommended neuron count'};
    tiledlayout(3,1, 'TileSpacing','compact', 'Padding','compact');
    for mm = 1:numel(metrics)
        nexttile
        matrix = summary_matrix(summary, area_names(aa), type_names, ...
            behavior_names, metrics{mm});
        imagesc(matrix)
        axis tight
        clim([0 max(1,n_animals)])
        colorbar
        yticks(1:numel(type_names))
        yticklabels(type_names)
        xticks(1:numel(behavior_names))
        xticklabels(behavior_names)
        xtickangle(60)
        title(titles{mm})
    end
    sgtitle(sprintf('%s: animal coverage (no effect estimates)', ...
        area_names(aa)))
    save_diagnostic_figure(figure_handle, output_folder, ...
        "coverage_" + area_names(aa));
end
end

function make_animal_count_figures(coverage, area_names, type_names, ...
        behavior_names, animal_ids, output_folder)
for aa = 1:numel(area_names)
    figure_handle = figure('Visible','off', 'Color','w', ...
        'Position',[50 50 1500 850]);
    tiledlayout(numel(type_names),1, ...
        'TileSpacing','compact', 'Padding','compact');
    for tt = 1:numel(type_names)
        matrix = nan(numel(animal_ids), numel(behavior_names));
        for animal_number = 1:numel(animal_ids)
            for behavior_number = 1:numel(behavior_names)
                mask = coverage.Area == area_names(aa) ...
                    & coverage.CellType == type_names(tt) ...
                    & coverage.AnimalID == animal_ids(animal_number) ...
                    & coverage.Behavior == ...
                    behavior_names(behavior_number);
                if any(mask)
                    matrix(animal_number,behavior_number) = ...
                        coverage.NTypeValid(mask);
                end
            end
        end
        nexttile
        imagesc(matrix)
        axis tight
        colorbar
        yticks(1:numel(animal_ids))
        yticklabels(animal_ids)
        xticks(1:numel(behavior_names))
        xticklabels(behavior_names)
        xtickangle(60)
        title(type_names(tt) + ": valid neurons per animal")
    end
    sgtitle(area_names(aa) + ": neuron-count audit")
    save_diagnostic_figure(figure_handle, output_folder, ...
        "neuron_counts_" + area_names(aa));
end
end

function make_precision_figures(summary, area_names, type_names, ...
        behavior_names, settings, output_folder)
for aa = 1:numel(area_names)
    figure_handle = figure('Visible','off', 'Color','w', ...
        'Position',[50 50 1500 500]);
    tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');

    nexttile
    width_matrix = summary_matrix(summary, area_names(aa), type_names, ...
        behavior_names, 'MedianObservedAI95Width');
    imagesc(width_matrix)
    axis tight
    colorbar
    clim([0 max(2, settings.maximum_acceptable_ai95_width*2)])
    yticklabels(type_names)
    yticks(1:numel(type_names))
    xticks(1:numel(behavior_names))
    xticklabels(behavior_names)
    xtickangle(60)
    title('Median bootstrap 95% width (log10 odds)')

    nexttile
    fraction_matrix = summary_matrix(summary, area_names(aa), ...
        type_names, behavior_names, 'StableAnimalFraction');
    imagesc(fraction_matrix)
    axis tight
    colorbar
    clim([0 1])
    yticklabels(type_names)
    yticks(1:numel(type_names))
    xticks(1:numel(behavior_names))
    xticklabels(behavior_names)
    xtickangle(60)
    title('Fraction of represented animals within precision limit')

    sgtitle(area_names(aa) + ...
        ": activation-index precision (effect location hidden)")
    save_diagnostic_figure(figure_handle, output_folder, ...
        "precision_" + area_names(aa));
end
end

function make_simulation_figures(simulation, area_names, type_names, ...
        behavior_names, settings, output_folder)
for aa = 1:numel(area_names)
    figure_handle = figure('Visible','off', 'Color','w', ...
        'Position',[50 50 1500 850]);
    tiledlayout(numel(type_names),1, ...
        'TileSpacing','compact', 'Padding','compact');
    for tt = 1:numel(type_names)
        matrix = nan(numel(settings.simulated_effect_sd), ...
            numel(behavior_names));
        for ee = 1:numel(settings.simulated_effect_sd)
            for bb = 1:numel(behavior_names)
                mask = simulation.Area == area_names(aa) ...
                    & simulation.CellType == type_names(tt) ...
                    & simulation.Behavior == behavior_names(bb) ...
                    & simulation.InjectedEffectSD == ...
                    settings.simulated_effect_sd(ee);
                if any(mask)
                    matrix(ee,bb) = ...
                        simulation.TwoSidedDetectionRate(mask);
                end
            end
        end
        nexttile
        imagesc(matrix)
        axis tight
        colorbar
        clim([0 1])
        yticks(1:numel(settings.simulated_effect_sd))
        yticklabels(string(settings.simulated_effect_sd))
        ylabel('Injected SD')
        xticks(1:numel(behavior_names))
        xticklabels(behavior_names)
        xtickangle(60)
        title(type_names(tt) + ": exact two-sided detection rate")
    end
    sgtitle(area_names(aa) + ...
        ": blinded animal-level operating characteristics")
    save_diagnostic_figure(figure_handle, output_folder, ...
        "simulation_" + area_names(aa));
end
end

function make_decision_figures(decision, area_names, type_names, ...
        behavior_names, output_folder)
decision_names = ["Insufficient coverage", ...
    "Insufficient animals for exact test", ...
    "Simulation type-I error failure", ...
    "Underpowered at 0.5 SD", ...
    "Use stratified group permutation", ...
    "Per-animal index feasible"];
for aa = 1:numel(area_names)
    matrix = nan(numel(type_names), numel(behavior_names));
    for tt = 1:numel(type_names)
        for bb = 1:numel(behavior_names)
            mask = decision.Area == area_names(aa) ...
                & decision.CellType == type_names(tt) ...
                & decision.Behavior == behavior_names(bb);
            if any(mask)
                matrix(tt,bb) = find(decision_names == ...
                    decision.Decision(mask), 1);
            end
        end
    end

    figure_handle = figure('Visible','off', 'Color','w', ...
        'Position',[50 50 1500 450]);
    imagesc(matrix)
    axis tight
    colormap([0.75 0.75 0.75; 0.55 0.55 0.55; ...
        0.70 0.20 0.20; 0.90 0.45 0.15; ...
        0.95 0.75 0.20; 0.20 0.65 0.35])
    colorbar('Ticks',1:numel(decision_names), ...
        'TickLabels',decision_names)
    clim([0.5 numel(decision_names)+0.5])
    yticks(1:numel(type_names))
    yticklabels(type_names)
    xticks(1:numel(behavior_names))
    xticklabels(behavior_names)
    xtickangle(60)
    title(area_names(aa) + ...
        ": prespecified feasibility decision (no observed effects)")
    save_diagnostic_figure(figure_handle, output_folder, ...
        "decision_" + area_names(aa));
end
end

function matrix = summary_matrix(summary, area_name, type_names, ...
        behavior_names, variable_name)
matrix = nan(numel(type_names), numel(behavior_names));
for tt = 1:numel(type_names)
    for bb = 1:numel(behavior_names)
        mask = summary.Area == area_name ...
            & summary.CellType == type_names(tt) ...
            & summary.Behavior == behavior_names(bb);
        if any(mask)
            values = summary.(variable_name)(mask);
            matrix(tt,bb) = values(1);
        end
    end
end
end

function save_diagnostic_figure(figure_handle, output_folder, stem)
fig_file = fullfile(output_folder, char(stem + ".fig"));
png_file = fullfile(output_folder, char(stem + ".png"));
savefig(figure_handle, fig_file);
try
    exportgraphics(figure_handle, png_file, 'Resolution', 200);
catch
    print(figure_handle, png_file, '-dpng', '-r200');
end
close(figure_handle)
end

function write_feasibility_report(report_file, checkpoint_file, settings, ...
        session_to_animal, decision)
file_id = fopen(report_file, 'w');
if file_id < 0
    warning('Could not write feasibility report: %s', report_file);
    return
end
cleanup = onCleanup(@() fclose(file_id));

fprintf(file_id, 'ACTIVATION-INDEX HIERARCHICAL FEASIBILITY AUDIT\n\n');
fprintf(file_id, 'Input checkpoint:\n%s\n\n', checkpoint_file);
fprintf(file_id, ['This audit does not export observed activation-index ' ...
    'means, directions, or biological p-values.\n\n']);
fprintf(file_id, 'Biological animals: %d\n', ...
    numel(unique(session_to_animal.AnimalID)));
fprintf(file_id, 'Sessions: %d\n\n', height(session_to_animal));
fprintf(file_id, 'Prespecified thresholds\n');
fprintf(file_id, '  Minimum animals: %d\n', settings.minimum_animals);
fprintf(file_id, '  Minimum neurons/type/animal: %d\n', ...
    settings.minimum_neurons_per_type);
fprintf(file_id, '  Recommended neurons/type/animal: %d\n', ...
    settings.recommended_neurons_per_type);
fprintf(file_id, '  Maximum activation-index 95%% width: %.3g\n', ...
    settings.maximum_acceptable_ai95_width);
fprintf(file_id, '  Maximum simulated null false-positive rate: %.3g\n', ...
    settings.maximum_null_false_positive_rate);
fprintf(file_id, '  Minimum simulated power at 0.5 SD: %.3g\n\n', ...
    settings.minimum_moderate_effect_power);

decision_names = unique(decision.Decision, 'stable');
fprintf(file_id, 'Decision counts (area x type x behavior comparisons)\n');
for ii = 1:numel(decision_names)
    fprintf(file_id, '  %s: %d\n', decision_names(ii), ...
        sum(decision.Decision == decision_names(ii)));
end

n_animals = numel(unique(session_to_animal.AnimalID));
fprintf(file_id, '\nExact sign-flip resolution with all animals\n');
fprintf(file_id, '  Minimum two-sided p: %.6f\n', ...
    minimum_exact_p(n_animals, "two-sided"));
fprintf(file_id, '  Minimum one-sided p: %.6f\n', ...
    minimum_exact_p(n_animals, "one-sided"));
fprintf(file_id, ['  A two-sided exact animal-level test cannot reach ' ...
    'alpha 0.05 unless at least 6 animals contribute.\n\n']);
fprintf(file_id, 'Review order\n');
fprintf(file_id, '  1. session_to_biological_animal.csv\n');
fprintf(file_id, '  2. coverage_*.png and neuron_counts_*.png\n');
fprintf(file_id, '  3. precision_*.png\n');
fprintf(file_id, '  4. simulation_*.png\n');
fprintf(file_id, '  5. decision_*.png and feasibility_decision_table.csv\n\n');
fprintf(file_id, ['Only after approving animal grouping, thresholds, and ' ...
    'coverage should the final biological analysis be run.\n']);
end
