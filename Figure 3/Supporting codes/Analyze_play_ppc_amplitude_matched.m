%% Play PPC after matching one-second epochs by delta amplitude/power
% This is the focused referee control. It does not estimate a PPC
% significance value and does not include non-play behaviors.
%
% Outputs:
%   matched_summary         one row per neuron
%   matched_draw_table      repeated amplitude-matched arrangements
%   null_draw_table         play-free pre-play partitions under no condition
%   ppc_power_spike_model   PPC change ~ power change + spike-count terms
%   session_paired_summary  one mean per session, within each plotted area
%   animal_paired_summary   one mean per animal; sessions of that animal pooled
%   paired_level_result     neuron, session-mean, and animal-mean signed-rank tests
%   neuron_count_relation   within-unit significance and effect size versus neuron count

data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
saving_folder = fullfile(data_root, 'Analysis results', ...
    'phase locking power control');
session_files = dir(fullfile(saving_folder, '*_phase_locking_power_control.mat'));
if isempty(session_files)
    error('No per-session power-control files in %s.', saving_folder);
end
power_control_struct = struct([]);
for ii = 1:numel(session_files)
    loaded_session = load(fullfile(session_files(ii).folder, ...
        session_files(ii).name), 'this_control');
    if isempty(power_control_struct)
        power_control_struct = loaded_session.this_control;
    else
        power_control_struct(end+1) = orderfields( ...
            loaded_session.this_control, power_control_struct); %#ok<SAGROW>
    end
end
phase_label_file = fullfile(data_root, 'Analysis results', ...
    'phase locking data', 'delta_all_neurons_v2.mat');
load(phase_label_file, 'all_neurons');

%% Hardcoded matching and reliability parameters
EPOCH_SECONDS = 1;
POWER_CALIPER_SD = 0.25;
TARGET_PLAY_FRACTION = 0.75;
N_DRAWS = 500;
MIN_SPIKES_PER_CONDITION = 20;
MIN_RELIABLE_MATCHED_DRAWS = 50;
MAX_ABS_POWER_SMD = 0.10;
LOCKING_ALPHA = 0.01;
PAIR_ALPHA = 0.05;
N_ANIMAL_BOOTSTRAP = 10000;
N_SIGN_FLIP = 100000;
RANDOM_SEED = 73491;
rng(RANDOM_SEED, 'twister');

% Spike counts never determine which epochs are matched. If a neuron has
% The initial duration is 75% of available play epochs, capped by matching
% capacity. If a neuron has too few spikes, duration is increased (using
% more independently power-matched epochs) until the criterion is met or
% all matchable epochs are exhausted. The power criterion is never relaxed.

total_neurons = sum(arrayfun(@(x) numel(x.neurons), power_control_struct));
summary_rows = {};
excluded_rows = {};
matched_blocks = cell(total_neurons,1);
null_blocks = cell(total_neurons,1);

for ss = 1:numel(power_control_struct)
    session_data = power_control_struct(ss);
    if ~isfield(session_data, 'pre_epoch_table') || ...
            ~isfield(session_data.neurons, 'pre_epoch_phase')
        error(['Saved structures do not contain one-second epochs. Re-run ' ...
            'Estimate_phase_locking_power_control.m with the updated GENERATE function.']);
    end
    if abs(session_data.options.epoch_duration-EPOCH_SECONDS) > eps
        error('Saved epoch duration is %.3g s, but analysis expects %.3g s.', ...
            session_data.options.epoch_duration, EPOCH_SECONDS);
    end

    animal = string(session_data.animal);
    session = string(session_data.session);
    for cc = 1:numel(session_data.channel_data)
        channel_data = session_data.channel_data(cc);
        channel = channel_data.channel;
        pre_power = channel_data.pre_epoch_power_db(:);
        play_power = channel_data.play_epoch_power_db(:);
        [power_pairs_index, power_caliper_db] = nearest_power_pairs( ...
            pre_power, play_power, POWER_CALIPER_SD);
        matching_capacity = size(power_pairs_index,1);
        valid_pre = find(isfinite(pre_power));
        valid_play = find(isfinite(play_power));

        neuron_indices = find([session_data.neurons.channel] == channel);
        for nn = reshape(neuron_indices,1,[])
            neuron = session_data.neurons(nn);
            identity = {animal, session, neuron.cluster_id, ...
                string(neuron.area), channel};
            label_index = string(all_neurons.session) == session & ...
                all_neurons.cluster_id == neuron.cluster_id;
            if sum(label_index) ~= 1
                excluded_rows(end+1,:) = [identity, ... %#ok<SAGROW>
                    {"Missing or duplicate original Figure 3 label", ...
                    matching_capacity, 0}];
                continue
            end
            if string(all_neurons.area(label_index)) ~= string(neuron.area)
                excluded_rows(end+1,:) = [identity, ... %#ok<SAGROW>
                    {"Area mismatch with original Figure 3 label", ...
                    matching_capacity, 0}];
                continue
            end
            original_ppc_p = all_neurons.EntireSession.PPCPval(label_index);
            original_angle = all_neurons.EntireSession.PreferedAngle(label_index);
            is_trough_locked = original_ppc_p < LOCKING_ALPHA && ...
                (original_angle < -pi/2 || original_angle > pi/2);
            if ~is_trough_locked
                continue
            end
            if matching_capacity < 2
                excluded_rows(end+1,:) = [identity, ... %#ok<SAGROW>
                    {"Fewer than two power-matchable epochs", matching_capacity, 0}];
                continue
            end

            % Both null groups come from pre-play epochs, so each can contain
            % at most half of that pool without reusing the same time.
            nonplay_null_capacity = floor(numel(valid_pre)/2);
            if nonplay_null_capacity < 2
                excluded_rows(end+1,:) = [identity, ... %#ok<SAGROW>
                    {"Fewer than four pre-play epochs for a play-free null", ...
                    matching_capacity, 0}];
                continue
            end
            duration_limit = min(matching_capacity, nonplay_null_capacity);
            initial_target_epochs = max(2, ...
                floor(TARGET_PLAY_FRACTION*numel(valid_play)));
            target_epochs = min(initial_target_epochs, duration_limit);
            expanded_for_spikes = false;
            expanded_for_balance = false;
            while true
                [matched_pre_index, matched_play_index] = ...
                    paired_draws(power_pairs_index, ...
                    target_epochs, N_DRAWS);
                matched_metrics = evaluate_draws(neuron, ...
                    matched_pre_index, matched_play_index, ...
                    pre_power, play_power, ...
                    MAX_ABS_POWER_SMD, true);
                n_reliable = sum(matched_metrics(:,14) == 1);
                reliable_metrics = matched_metrics(matched_metrics(:,14) == 1,:);
                enough_power_draws = n_reliable >= min( ...
                    MIN_RELIABLE_MATCHED_DRAWS, size(matched_metrics,1));
                enough_spikes = ~isempty(reliable_metrics) && ...
                    median(reliable_metrics(:,7), 'omitnan') >= ...
                    MIN_SPIKES_PER_CONDITION && ...
                    median(reliable_metrics(:,8), 'omitnan') >= ...
                    MIN_SPIKES_PER_CONDITION;
                if (enough_power_draws && enough_spikes) || ...
                        target_epochs >= duration_limit
                    break
                end
                expanded_for_spikes = expanded_for_spikes || ~enough_spikes;
                expanded_for_balance = expanded_for_balance || ...
                    ~enough_power_draws;
                target_epochs = min(duration_limit, ...
                    max(target_epochs+1, ceil(1.5*target_epochs)));
            end

            reliable_matched = matched_metrics(:,14) == 1;
            reliable_values = matched_metrics(reliable_matched,:);
            required_draws = min(MIN_RELIABLE_MATCHED_DRAWS, ...
                size(matched_metrics,1));
            enough_spikes = ~isempty(reliable_values) && ...
                median(reliable_values(:,7), 'omitnan') >= ...
                MIN_SPIKES_PER_CONDITION && ...
                median(reliable_values(:,8), 'omitnan') >= ...
                MIN_SPIKES_PER_CONDITION;
            if sum(reliable_matched) < required_draws || ~enough_spikes
                excluded_rows(end+1,:) = [identity, ... %#ok<SAGROW>
                    {"Insufficient spikes or residual power imbalance", ...
                    matching_capacity, sum(reliable_matched)}];
                continue
            end

            [null_pre_index, null_post_index] = merged_null_draws( ...
                valid_pre, target_epochs, N_DRAWS);
            null_metrics = evaluate_null_draws(neuron, ...
                null_pre_index, null_post_index, ...
                pre_power, play_power);

            if size(reliable_values,1) > 1
                ppc_draw_percentiles = prctile( ...
                    reliable_values(:,13), [2.5 97.5]);
            else
                ppc_draw_percentiles = [NaN NaN];
            end
            neuron_row = size(summary_rows,1) + 1;
            neuron_key = session + "_ch" + string(channel) + ...
                "_cl" + string(neuron.cluster_id);
            channel_key = session + "_ch" + string(channel);
            cell_type = "Trough";
            summary_rows(end+1,:) = {animal, session, neuron_key, ... %#ok<SAGROW>
                channel_key, neuron.cluster_id, string(neuron.area), ...
                cell_type, channel, numel(valid_play), ...
                initial_target_epochs, target_epochs, ...
                target_epochs/numel(valid_play), matching_capacity, ...
                nonplay_null_capacity, ...
                power_caliper_db, sum(reliable_matched), ...
                expanded_for_spikes, expanded_for_balance, ...
                median(reliable_values(:,3), 'omitnan'), ...
                median(reliable_values(:,4), 'omitnan'), ...
                median(reliable_values(:,5), 'omitnan'), ...
                median(reliable_values(:,7), 'omitnan'), ...
                median(reliable_values(:,8), 'omitnan'), ...
                median(reliable_values(:,9), 'omitnan'), ...
                median(reliable_values(:,10), 'omitnan'), ...
                median(reliable_values(:,11), 'omitnan'), ...
                median(reliable_values(:,12), 'omitnan'), ...
                median(reliable_values(:,13), 'omitnan'), ...
                ppc_draw_percentiles(1), ppc_draw_percentiles(2), ...
                median(abs(reliable_values(:,6)), 'omitnan'), ...
                NaN, NaN, NaN, NaN};

            matched_blocks{neuron_row} = [ ...
                repmat(neuron_row,size(matched_metrics,1),1), ...
                matched_metrics];
            null_blocks{neuron_row} = [ ...
                repmat(neuron_row,size(null_metrics,1),1), ...
                null_metrics];
        end
    end
end

if isempty(summary_rows)
    error('No neurons passed the amplitude-balance and spike-count criteria.');
end

summary_names = {'Animal','Session','NeuronKey','ChannelKey','ClusterID', ...
    'Area','CellType','Channel','AvailablePlayEpochs', ...
    'InitialTargetEpochs','EpochsPerCondition','AnalyzedPlayFraction', ...
    'MatchingCapacity','NonPlayNullCapacity', ...
    'PowerCaliperDB','ReliableMatchedDraws','ExpandedForSpikes', ...
    'ExpandedForPowerBalance','PrePowerDB','PlayPowerDB', ...
    'PowerDifferenceDB','PreSpikeCount','PlaySpikeCount', ...
    'SpikeCountDifference','MeanSpikeCount','PrePPC','PlayPPC', ...
    'PPCDifference','PPCDifferenceDrawP2_5','PPCDifferenceDrawP97_5', ...
    'MedianAbsolutePowerSMD','ModelExpectedPPCDifference', ...
    'ModelExpectedPPCDifferenceCI95Low', ...
    'ModelExpectedPPCDifferenceCI95High','ModelResidual'};
matched_summary = cell2table(summary_rows, 'VariableNames', summary_names);

matched_draw_names = {'NeuronRow','Draw','EpochCount','PrePowerDB','PlayPowerDB', ...
    'PowerDifferenceDB','PowerSMD','PreSpikeCount','PlaySpikeCount', ...
    'SpikeCountDifference','MeanSpikeCount','PrePPC','PlayPPC', ...
    'PPCDifference','Reliable'};
matched_draw_table = array2table(vertcat(matched_blocks{:}), ...
    'VariableNames', matched_draw_names);
null_draw_names = {'NeuronRow','Draw','EpochCount', ...
    'PseudoPrePowerDB','PseudoPostPowerDB','PowerDifferenceDB','PowerSMD', ...
    'PseudoPreSpikeCount','PseudoPostSpikeCount','SpikeCountDifference', ...
    'MeanSpikeCount','PseudoPrePPC','PseudoPostPPC','PPCDifference','Reliable'};
null_draw_table = array2table(vertcat(null_blocks{:}), ...
    'VariableNames', null_draw_names);

excluded_names = {'Animal','Session','ClusterID','Area','Channel', ...
    'Reason','MatchingCapacity','ReliableDraws'};
if isempty(excluded_rows)
    excluded_neurons = cell2table(cell(0,numel(excluded_names)), ...
        'VariableNames', excluded_names);
else
    excluded_neurons = cell2table(excluded_rows, ...
        'VariableNames', excluded_names);
end

%% Model PPC changes under the play-free no-condition null
% MeanSpikeCount is included because PPC variance depends on the available
% spike count; SpikeCountDifference captures condition imbalance requested
% in the proposed model. EpochCount records adaptive duration. Both
% pseudo-conditions are sampled from pre-play epochs only, so neither group
% contains play and their play proportion cannot explain PPC differences.
% Because draws reuse epochs, use this model for conditional prediction;
% its nominal coefficient p-values are not independent-resample inference.
ppc_power_spike_model = [];
model_rows = null_draw_table.Reliable == 1;
model_table = null_draw_table(model_rows,:);
if height(model_table) >= 100
    row_index = model_table.NeuronRow;
    model_table.Animal = categorical(matched_summary.Animal(row_index));
    model_table.NeuronKey = categorical(matched_summary.NeuronKey(row_index));
    try
        ppc_power_spike_model = fitlme(model_table, ...
            ['PPCDifference ~ PowerDifferenceDB + SpikeCountDifference + ' ...
             'MeanSpikeCount + EpochCount + (1|Animal) + (1|NeuronKey)']);

        prediction_table = table;
        prediction_table.PowerDifferenceDB = ...
            matched_summary.PowerDifferenceDB;
        prediction_table.SpikeCountDifference = ...
            matched_summary.SpikeCountDifference;
        prediction_table.MeanSpikeCount = matched_summary.MeanSpikeCount;
        prediction_table.EpochCount = matched_summary.EpochsPerCondition;
        % predict requires the response column but ignores its values.
        prediction_table.PPCDifference = zeros(height(matched_summary),1);
        prediction_table.Animal = categorical(matched_summary.Animal, ...
            categories(model_table.Animal));
        prediction_table.NeuronKey = categorical(matched_summary.NeuronKey, ...
            categories(model_table.NeuronKey));
        predictable = ~isundefined(prediction_table.Animal) & ...
            ~isundefined(prediction_table.NeuronKey);
        [expected_ppc, expected_ppc_ci] = predict( ...
            ppc_power_spike_model, prediction_table(predictable,:), ...
            'Conditional', false, 'Alpha', 0.05, ...
            'Prediction', 'curve');
        matched_summary.ModelExpectedPPCDifference(predictable) = expected_ppc;
        matched_summary.ModelExpectedPPCDifferenceCI95Low(predictable) = ...
            expected_ppc_ci(:,1);
        matched_summary.ModelExpectedPPCDifferenceCI95High(predictable) = ...
            expected_ppc_ci(:,2);
        matched_summary.ModelResidual = matched_summary.PPCDifference - ...
            matched_summary.ModelExpectedPPCDifference;
        disp(ppc_power_spike_model)
    catch exception
        warning('PPC power/spike-count model failed: %s', exception.message);
    end
end

%% Overall trough-locked-cell result, using animals as replicates
valid_result = isfinite(matched_summary.PPCDifference) & ...
    isfinite(matched_summary.ModelExpectedPPCDifference) & ...
    isfinite(matched_summary.ModelResidual);
if any(valid_result)
    [animal_group, animal_name] = findgroups( ...
        matched_summary.Animal(valid_result));
    n_neurons = splitapply(@numel, ...
        matched_summary.ModelResidual(valid_result), animal_group);
    animal_observed = splitapply(@(x) mean(x,'omitnan'), ...
        matched_summary.PPCDifference(valid_result), animal_group);
    animal_expected = splitapply(@(x) mean(x,'omitnan'), ...
        matched_summary.ModelExpectedPPCDifference(valid_result), animal_group);
    animal_residual = splitapply(@(x) mean(x,'omitnan'), ...
        matched_summary.ModelResidual(valid_result), animal_group);
    animal_summary = table(animal_name, n_neurons, animal_observed, ...
        animal_expected, animal_residual, ...
        'VariableNames', {'Animal','NNeurons','ObservedPPCDifference', ...
        'ExpectedNullPPCDifference','ResidualPPCDifference'});

    bootstrap_mean = nan(N_ANIMAL_BOOTSTRAP,1);
    n_animals = height(animal_summary);
    for bb = 1:N_ANIMAL_BOOTSTRAP
        sampled = randi(n_animals,n_animals,1);
        bootstrap_mean(bb) = mean(animal_residual(sampled),'omitnan');
    end
    residual_ci = prctile(bootstrap_mean,[2.5 97.5]);
    residual_p = one_sided_sign_flip(animal_residual,N_SIGN_FLIP);
    overall_result = table(n_animals, sum(n_neurons), ...
        mean(animal_observed,'omitnan'), ...
        mean(animal_expected,'omitnan'), ...
        mean(animal_residual,'omitnan'), residual_ci(1), residual_ci(2), ...
        residual_p, ...
        'VariableNames', {'NAnimals','NNeurons', ...
        'ObservedPPCDifference','ExpectedNullPPCDifference', ...
        'ResidualPPCDifference','ResidualCI95Low','ResidualCI95High', ...
        'OneSidedSignFlipP'});
else
    animal_summary = table;
    overall_result = table(NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
        'VariableNames', {'NAnimals','NNeurons', ...
        'ObservedPPCDifference','ExpectedNullPPCDifference', ...
        'ResidualPPCDifference','ResidualCI95Low','ResidualCI95High', ...
        'OneSidedSignFlipP'});
end
disp(overall_result)

%% Diagnostics
figure('Name', 'Amplitude-matched play PPC control');
subplot(1,3,1)
histogram(matched_summary.PowerDifferenceDB, 20)
xline(0, 'r-')
xlabel('Matched play - pre-play delta power (dB)')
ylabel('Neurons')

subplot(1,3,2)
plot(matched_summary.PowerDifferenceDB, ...
    matched_summary.PPCDifference, 'k.')
xline(0, 'k:'); yline(0, 'k:');
xlabel('\Delta delta power (dB)')
ylabel('\Delta PPC')

subplot(1,3,3)
reliable_null = null_draw_table(null_draw_table.Reliable == 1,:);
plot(reliable_null.PowerDifferenceDB, ...
    reliable_null.PPCDifference, '.', 'Color', [0.7 0.7 0.7])
hold on
plot(matched_summary.PowerDifferenceDB, ...
    matched_summary.PPCDifference, 'r.')
xlabel('\Delta delta power (dB)')
ylabel('\Delta PPC')
legend({'Merged-pool null arrangements','Amplitude-matched estimate'})

%% Power-balanced trough-cell PPC, in the Figure 3e square-plot layout
% Figure 3 labels superior colliculus as SC. The stored area name is SupCol.
paper_areas = ["SupCol","DLPAG","LPAG","VLPAG","DR"];
paper_labels = ["SC","DLPAG","LPAG","VLPAG","DR"];
finite_ppc = isfinite(matched_summary.PrePPC) & isfinite(matched_summary.PlayPPC);
plot_rows = matched_summary(finite_ppc & ...
    ismember(string(matched_summary.Area), paper_areas), :);
pre_z = atanh(min(max(plot_rows.PrePPC, -0.999), 0.999));
play_z = atanh(min(max(plot_rows.PlayPPC, -0.999), 0.999));
present_areas = string(plot_rows.Area);
fprintf('Areas included in the all-areas panel:\n');
disp(unique(present_areas));
axis_limit = [min([-0.005; pre_z; play_z]), max([0.07; pre_z; play_z])];
figure('Name', 'Power-balanced trough PPC, pre-play vs play');
for an = 1:numel(paper_areas)
    subplot(1, numel(paper_areas)+1, an)
    in_area = present_areas == paper_areas(an);
    plot_paired_ppc(pre_z(in_area), play_z(in_area), axis_limit)
    title(sprintf('%s n=%d p=%s', paper_labels(an), sum(in_area), ...
        paired_signrank_text(pre_z(in_area), play_z(in_area))))
    xlabel('Pre-play PPC (atanh)')
    if an == 1
        ylabel('Play PPC (atanh)')
    end
end
subplot(1, numel(paper_areas)+1, numel(paper_areas)+1)
plot_paired_ppc(pre_z, play_z, axis_limit)
xlabel('Pre-play PPC (atanh)')
title(sprintf('ALL AREAS n=%d neurons, %d animals, p=%s', ...
    height(plot_rows), numel(unique(string(plot_rows.Animal))), ...
    paired_signrank_text(pre_z, play_z)))

%% Session means and animal means of the same paired comparison
% Each neuron pair stays on the atanh scale used above. A session mean
% averages the neurons recorded in that session. An animal mean averages
% every neuron from that animal, including neurons from its other sessions.
% The signed-rank p-value in each title treats those means as the pairs.
% A circle is filled when the neuron pairs inside that session or animal
% themselves reject zero at PAIR_ALPHA.
plot_rows.AnimalID = batch_animal_id(plot_rows.Session);
unparsed_animal = plot_rows.AnimalID == string(plot_rows.Session);
if any(unparsed_animal)
    warning('Analyze_play_ppc_amplitude_matched:AnimalId', ...
        ['%d sessions did not match the Batch animal-name pattern. ' ...
        'Each of those sessions is kept as its own animal.'], ...
        numel(unique(plot_rows.Session(unparsed_animal))));
end
animal_names = unique(plot_rows.AnimalID);
animal_colors = paired_animal_colors(numel(animal_names));
area_names = ["ALL"; paper_areas(:)];
area_labels = ["ALL"; paper_labels(:)];
session_blocks = cell(numel(area_names), 1);
animal_blocks = cell(numel(area_names), 1);
level_blocks = cell(numel(area_names), 1);
relation_blocks = cell(numel(area_names)*2, 1);
for an = 1:numel(area_names)
    if area_names(an) == "ALL"
        in_area = true(height(plot_rows), 1);
    else
        in_area = present_areas == area_names(an);
    end
    session_blocks{an} = summarize_session_pairs(area_labels(an), ...
        plot_rows.Session(in_area), plot_rows.AnimalID(in_area), ...
        pre_z(in_area), play_z(in_area), PAIR_ALPHA);
    animal_blocks{an} = summarize_animal_pairs(area_labels(an), ...
        plot_rows.AnimalID(in_area), plot_rows.Session(in_area), ...
        pre_z(in_area), play_z(in_area), PAIR_ALPHA);
    level_blocks{an} = paired_level_rows(area_labels(an), ...
        pre_z(in_area), play_z(in_area), ...
        session_blocks{an}, animal_blocks{an});
    relation_blocks{2*an-1} = one_neuron_count_relation( ...
        "Session", area_labels(an), session_blocks{an});
    relation_blocks{2*an} = one_neuron_count_relation( ...
        "Animal", area_labels(an), animal_blocks{an});
end
session_paired_summary = vertcat(session_blocks{:});
animal_paired_summary = vertcat(animal_blocks{:});
paired_level_result = vertcat(level_blocks{:});
neuron_count_relation = vertcat(relation_blocks{:});
disp(paired_level_result)
fprintf(['Filled circles are within-session or within-animal signed-rank ' ...
    'tests at p < %.2g. A signed-rank test needs several pairs before it ' ...
    'can cross that line, so significance can track neuron count even when ' ...
    'the paired difference does not. The Spearman row asks whether the ' ...
    'mean difference itself grows with the neuron count.\n'], PAIR_ALPHA);
disp(neuron_count_relation)

plot_paired_unit_figure('Power-balanced PPC by session', ...
    paper_labels, session_paired_summary, animal_names, animal_colors, ...
    axis_limit, 'session-mean', PAIR_ALPHA);
plot_paired_unit_figure('Power-balanced PPC by animal', ...
    paper_labels, animal_paired_summary, animal_names, animal_colors, ...
    axis_limit, 'animal-mean', PAIR_ALPHA);
plot_neuron_count_figure(session_paired_summary, animal_paired_summary, ...
    animal_names, animal_colors, neuron_count_relation, PAIR_ALPHA);

analysis_parameters = struct( ...
    'epoch_seconds', EPOCH_SECONDS, ...
    'power_caliper_sd', POWER_CALIPER_SD, ...
    'target_play_fraction', TARGET_PLAY_FRACTION, ...
    'n_draws', N_DRAWS, ...
    'min_spikes_per_condition', MIN_SPIKES_PER_CONDITION, ...
    'min_reliable_matched_draws', MIN_RELIABLE_MATCHED_DRAWS, ...
    'max_abs_power_smd', MAX_ABS_POWER_SMD, ...
    'locking_alpha_from_original_analysis', LOCKING_ALPHA, ...
    'within_unit_pair_alpha', PAIR_ALPHA, ...
    'n_animal_bootstrap', N_ANIMAL_BOOTSTRAP, ...
    'n_sign_flip', N_SIGN_FLIP, ...
    'random_seed', RANDOM_SEED);
save(fullfile(saving_folder, ...
    'play_ppc_amplitude_matched_analysis.mat'), ...
    'matched_summary', 'matched_draw_table', 'null_draw_table', ...
    'excluded_neurons', 'ppc_power_spike_model', ...
    'animal_summary', 'overall_result', ...
    'session_paired_summary', 'animal_paired_summary', ...
    'paired_level_result', 'neuron_count_relation', ...
    'analysis_parameters', '-v7.3');

%% Local functions
function [pairs, caliper_db] = nearest_power_pairs( ...
    pre_power, play_power, caliper_sd)
valid_pre = find(isfinite(pre_power));
valid_play = find(isfinite(play_power));
combined = [pre_power(valid_pre); play_power(valid_play)];
if isempty(valid_pre) || isempty(valid_play)
    pairs = zeros(0,2);
    caliper_db = NaN;
    return
end
power_scale = std(combined, 0, 'omitnan');
if ~isfinite(power_scale) || power_scale == 0
    power_scale = eps(max(abs(combined)));
end
caliper_db = caliper_sd * power_scale;
cost = abs(pre_power(valid_pre) - play_power(valid_play)');
cost(cost > caliper_db) = Inf;
if isempty(cost) || all(isinf(cost), 'all')
    pairs = zeros(0,2);
    return
end
% Minimum-cost one-to-one pairs; unmatched epochs are allowed.
[local_pairs,~] = matchpairs(cost, 2*caliper_db, 'min');
pairs = [valid_pre(local_pairs(:,1)), valid_play(local_pairs(:,2))];
end

function [pre_index, play_index] = paired_draws(pairs, target, n_draws)
capacity = size(pairs,1);
selected_sets = unique_subsets(capacity, target, n_draws);
pre_index = reshape(pairs(selected_sets,1), size(selected_sets));
play_index = reshape(pairs(selected_sets,2), size(selected_sets));
end

function [pseudo_pre, pseudo_post] = merged_null_draws( ...
    valid_merged, k, n_requested)
n_epochs = numel(valid_merged);
if 2*k > n_epochs
    error('The merged epoch pool cannot supply two groups of %d epochs.', k);
end

% Number of oriented assignments: choose K pseudo-pre epochs, then K
% pseudo-post epochs from those remaining.
log_assignment_count = gammaln(n_epochs+1) - gammaln(k+1) - ...
    gammaln(n_epochs-k+1) + gammaln(n_epochs-k+1) - ...
    gammaln(k+1) - gammaln(n_epochs-2*k+1);

if log_assignment_count <= log(n_requested)
    pre_combinations = nchoosek(1:n_epochs, k);
    local_pre = zeros(0,k);
    local_post = zeros(0,k);
    for rr = 1:size(pre_combinations,1)
        remaining = setdiff(1:n_epochs, pre_combinations(rr,:));
        post_combinations = nchoosek(remaining, k);
        local_pre = [local_pre; repmat(pre_combinations(rr,:), ...
            size(post_combinations,1),1)]; %#ok<AGROW>
        local_post = [local_post; post_combinations]; %#ok<AGROW>
    end
else
    % Draw unique unordered partitions, then add each label-swapped
    % complement. This enforces a symmetric no-condition null.
    n_base = floor(n_requested/2);
    canonical_keys = zeros(0,2*k);
    base_pre = zeros(0,k);
    base_post = zeros(0,k);
    while size(base_pre,1) < n_base
        selected = randperm(n_epochs, 2*k);
        this_pre = sort(selected(1:k));
        this_post = sort(selected(k+1:end));
        key = canonical_partition(this_pre, this_post);
        if isempty(canonical_keys) || ~ismember(key, canonical_keys, 'rows')
            canonical_keys(end+1,:) = key; %#ok<AGROW>
            base_pre(end+1,:) = this_pre; %#ok<AGROW>
            base_post(end+1,:) = this_post; %#ok<AGROW>
        end
    end
    local_pre = [base_pre; base_post];
    local_post = [base_post; base_pre];
end

pseudo_pre = reshape(valid_merged(local_pre), size(local_pre));
pseudo_post = reshape(valid_merged(local_post), size(local_post));
end

function key = canonical_partition(first_group, second_group)
forward = [first_group, second_group];
reverse = [second_group, first_group];
first_difference = find(forward ~= reverse, 1);
if isempty(first_difference) || ...
        forward(first_difference) < reverse(first_difference)
    key = forward;
else
    key = reverse;
end
end

function subsets = unique_subsets(population_size, subset_size, n_requested)
log_count = gammaln(population_size+1) - gammaln(subset_size+1) - ...
    gammaln(population_size-subset_size+1);
if log_count <= log(n_requested)
    subsets = nchoosek(1:population_size, subset_size);
    return
end
subsets = zeros(0,subset_size);
while size(subsets,1) < n_requested
    batch_size = 2*(n_requested-size(subsets,1));
    batch = nan(batch_size,subset_size);
    for rr = 1:batch_size
        batch(rr,:) = sort(randperm(population_size,subset_size));
    end
    subsets = unique([subsets; batch], 'rows', 'stable');
end
subsets = subsets(1:n_requested,:);
end

function metrics = evaluate_draws(neuron, pre_index, play_index, ...
    pre_power, play_power, max_abs_smd, require_balance)
pre_count = neuron.pre_epoch_spike_count(:);
play_count = neuron.play_epoch_spike_count(:);
pre_vector = cellfun(@phase_vector, neuron.pre_epoch_phase);
play_vector = cellfun(@phase_vector, neuron.play_epoch_phase);

pre_epochs = values_by_draw(pre_count, pre_index);
play_epochs = values_by_draw(play_count, play_index);
pre_phase = values_by_draw(pre_vector, pre_index);
play_phase = values_by_draw(play_vector, play_index);
pre_power_by_draw = values_by_draw(pre_power, pre_index);
play_power_by_draw = values_by_draw(play_power, play_index);

n_pre = sum(pre_epochs,2);
n_play = sum(play_epochs,2);
r_pre = sum(pre_phase,2);
r_play = sum(play_phase,2);
mean_pre_power = mean(pre_power_by_draw,2,'omitnan');
mean_play_power = mean(play_power_by_draw,2,'omitnan');
pooled_sd = sqrt((var(pre_power_by_draw,0,2,'omitnan') + ...
    var(play_power_by_draw,0,2,'omitnan'))/2);
n_pre = n_pre(:);
n_play = n_play(:);
r_pre = r_pre(:);
r_play = r_play(:);
mean_pre_power = mean_pre_power(:);
mean_play_power = mean_play_power(:);
pooled_sd = pooled_sd(:);
pre_ppc = ppc_from_sufficient_statistics(r_pre, n_pre);
play_ppc = ppc_from_sufficient_statistics(r_play, n_play);
power_difference = mean_play_power-mean_pre_power;
spike_difference = n_play-n_pre;
mean_spike_count = (n_play+n_pre)/2;
ppc_difference = play_ppc-pre_ppc;
power_smd = power_difference./pooled_sd;
power_smd(pooled_sd == 0 & mean_play_power == mean_pre_power) = 0;
reliable = isfinite(pre_ppc) & isfinite(play_ppc);
if require_balance
    reliable = reliable & abs(power_smd) <= max_abs_smd;
end
reliable = double(reliable);

n_draws = size(pre_index,1);
metrics = [(1:n_draws)', repmat(size(pre_index,2),n_draws,1), ...
    mean_pre_power, mean_play_power, power_difference, power_smd, ...
    n_pre, n_play, spike_difference, mean_spike_count, ...
    pre_ppc, play_ppc, ppc_difference, reliable];
end

function metrics = evaluate_null_draws(neuron, pseudo_pre, pseudo_post, ...
    pre_power, play_power)
merged_power = [pre_power(:); play_power(:)];
merged_count = [neuron.pre_epoch_spike_count(:); ...
    neuron.play_epoch_spike_count(:)];
merged_vector = [cellfun(@phase_vector, neuron.pre_epoch_phase); ...
    cellfun(@phase_vector, neuron.play_epoch_phase)];

pre_epochs = values_by_draw(merged_count, pseudo_pre);
post_epochs = values_by_draw(merged_count, pseudo_post);
pre_phase = values_by_draw(merged_vector, pseudo_pre);
post_phase = values_by_draw(merged_vector, pseudo_post);
pre_power_by_draw = values_by_draw(merged_power, pseudo_pre);
post_power_by_draw = values_by_draw(merged_power, pseudo_post);

n_pre = sum(pre_epochs,2);
n_post = sum(post_epochs,2);
r_pre = sum(pre_phase,2);
r_post = sum(post_phase,2);
mean_pre_power = mean(pre_power_by_draw,2,'omitnan');
mean_post_power = mean(post_power_by_draw,2,'omitnan');
pooled_sd = sqrt((var(pre_power_by_draw,0,2,'omitnan') + ...
    var(post_power_by_draw,0,2,'omitnan'))/2);
n_pre = n_pre(:);
n_post = n_post(:);
r_pre = r_pre(:);
r_post = r_post(:);
mean_pre_power = mean_pre_power(:);
mean_post_power = mean_post_power(:);
pooled_sd = pooled_sd(:);
pre_ppc = ppc_from_sufficient_statistics(r_pre, n_pre);
post_ppc = ppc_from_sufficient_statistics(r_post, n_post);
power_difference = mean_post_power-mean_pre_power;
spike_difference = n_post-n_pre;
mean_spike_count = (n_post+n_pre)/2;
ppc_difference = post_ppc-pre_ppc;
power_smd = power_difference./pooled_sd;
power_smd(pooled_sd == 0 & mean_post_power == mean_pre_power) = 0;
reliable = double(isfinite(pre_ppc) & isfinite(post_ppc));

n_draws = size(pseudo_pre,1);
metrics = [(1:n_draws)', repmat(size(pseudo_pre,2),n_draws,1), ...
    mean_pre_power, mean_post_power, power_difference, power_smd, ...
    n_pre, n_post, spike_difference, mean_spike_count, ...
    pre_ppc, post_ppc, ppc_difference, reliable];
end

function selected = values_by_draw(data, index)
% Indexing a vector with a vector keeps the data orientation. A single
% draw is a 1-by-K row, so indexing a column of epoch values returns a
% K-by-1 column. Reducing along dimension 2 then leaves K rows, which
% cannot be concatenated with the one-row draw index. Restore draws-by-
% epochs before summing, averaging, or taking the variance.
selected = reshape(data(index(:)), size(index));
end

function label = paired_signrank_text(pre_z, play_z)
p_value = paired_signrank_p(pre_z, play_z);
if ~isfinite(p_value)
    label = "NA";
    return
end
label = sprintf('%.3g', p_value);
end

function p_value = paired_signrank_p(pre_z, play_z)
pre_z = pre_z(:);
play_z = play_z(:);
finite_pairs = isfinite(pre_z) & isfinite(play_z);
if nnz(finite_pairs) < 1
    p_value = NaN;
    return
end
try
    p_value = signrank(pre_z(finite_pairs), play_z(finite_pairs));
catch
    p_value = NaN;
end
end

function plot_paired_ppc(pre_z, play_z, axis_limit)
plot(pre_z, play_z, 'k.')
hold on
finite_pairs = isfinite(pre_z) & isfinite(play_z);
if any(finite_pairs)
    pair_median = [median(pre_z(finite_pairs)), median(play_z(finite_pairs))];
    plot(pair_median(1), pair_median(2), 'rx', 'LineWidth', 1.2)
end
plot(axis_limit, axis_limit, 'r-')
axis([axis_limit axis_limit])
axis square
end

function value = phase_vector(phase)
value = sum(exp(1i*phase(:)));
end

function value = ppc_from_sufficient_statistics(resultant, n)
value = (abs(resultant).^2-n)./(n.*(n-1));
value(n < 2) = NaN;
end

function p_value = one_sided_sign_flip(values, n_random)
values = values(isfinite(values));
n = numel(values);
if n == 0
    p_value = NaN;
    return
end
observed_mean = mean(values);
if n <= 15
    n_assignments = 2^n;
    null_mean = nan(n_assignments,1);
    for ii = 0:n_assignments-1
        signs = 2*bitget(ii,1:n)-1;
        null_mean(ii+1) = mean(signs(:).*values(:));
    end
    p_value = sum(null_mean >= observed_mean) / n_assignments;
else
    null_mean = nan(n_random,1);
    for ii = 1:n_random
        signs = 2*(rand(n,1) > 0.5)-1;
        null_mean(ii) = mean(signs.*values(:));
    end
    p_value = (1+sum(null_mean >= observed_mean))/(n_random+1);
end
end

function animal_id = batch_animal_id(session_names)
% Same Batch id as the area-limit lookup: sessions of one animal share it.
session_names = string(session_names);
animal_id = strings(size(session_names));
for ii = 1:numel(session_names)
    animal_id(ii) = one_batch_animal_id(session_names(ii));
end
end

function animal_id = one_batch_animal_id(session_name)
parts = strsplit(strtrim(char(session_name)));
if numel(parts) < 3
    animal_id = string(strtrim(session_name));
    return
end
animal_batch = parts{1};
repeated_animal = parts{3};
if numel(animal_batch) < 2
    animal_id = string(strtrim(session_name));
    return
end
if strcmp(repeated_animal, 'Single2')
    animal_id = string(['Batch', animal_batch(2), repeated_animal]);
    return
end
if numel(animal_batch) < 4
    animal_id = string(strtrim(session_name));
    return
end
animal_id = string(['Batch', animal_batch(2), repeated_animal, animal_batch(4)]);
end

function summary = summarize_session_pairs(area_label, session_id, ...
    animal_id, pre_z, play_z, pair_alpha)
names = {'Area','Session','AnimalID','NNeurons','PreMean','PlayMean', ...
    'MeanDifference','SignrankP','Significant'};
session_id = string(session_id(:));
animal_id = string(animal_id(:));
pre_z = pre_z(:);
play_z = play_z(:);
if isempty(session_id)
    summary = table(strings(0,1), strings(0,1), strings(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
        false(0,1), 'VariableNames', names);
    return
end
[group_index, group_names] = findgroups(session_id);
n_groups = numel(group_names);
animal_out = strings(n_groups, 1);
n_neurons = zeros(n_groups, 1);
pre_mean = nan(n_groups, 1);
play_mean = nan(n_groups, 1);
mean_difference = nan(n_groups, 1);
signrank_p = nan(n_groups, 1);
significant = false(n_groups, 1);
for gg = 1:n_groups
    in_group = group_index == gg;
    stats = paired_group_stats(pre_z(in_group), play_z(in_group), pair_alpha);
    these_animals = unique(animal_id(in_group));
    animal_out(gg) = these_animals(1);
    n_neurons(gg) = stats.n_neurons;
    pre_mean(gg) = stats.pre_mean;
    play_mean(gg) = stats.play_mean;
    mean_difference(gg) = stats.mean_difference;
    signrank_p(gg) = stats.p_value;
    significant(gg) = stats.significant;
end
summary = table(repmat(string(area_label), n_groups, 1), ...
    string(group_names(:)), animal_out, n_neurons, pre_mean, play_mean, ...
    mean_difference, signrank_p, significant, 'VariableNames', names);
end

function summary = summarize_animal_pairs(area_label, animal_id, ...
    session_id, pre_z, play_z, pair_alpha)
names = {'Area','AnimalID','NSessions','NNeurons','PreMean','PlayMean', ...
    'MeanDifference','SignrankP','Significant','Sessions'};
animal_id = string(animal_id(:));
session_id = string(session_id(:));
pre_z = pre_z(:);
play_z = play_z(:);
if isempty(animal_id)
    summary = table(strings(0,1), strings(0,1), zeros(0,1), zeros(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), ...
        strings(0,1), 'VariableNames', names);
    return
end
[group_index, group_names] = findgroups(animal_id);
n_groups = numel(group_names);
n_sessions = zeros(n_groups, 1);
n_neurons = zeros(n_groups, 1);
pre_mean = nan(n_groups, 1);
play_mean = nan(n_groups, 1);
mean_difference = nan(n_groups, 1);
signrank_p = nan(n_groups, 1);
significant = false(n_groups, 1);
sessions = strings(n_groups, 1);
for gg = 1:n_groups
    in_group = group_index == gg;
    stats = paired_group_stats(pre_z(in_group), play_z(in_group), pair_alpha);
    these_sessions = unique(session_id(in_group));
    these_sessions = these_sessions(these_sessions ~= "" & ~ismissing(these_sessions));
    n_sessions(gg) = numel(these_sessions);
    if isempty(these_sessions)
        sessions(gg) = "";
    else
        sessions(gg) = string(strjoin(cellstr(these_sessions), ', '));
    end
    n_neurons(gg) = stats.n_neurons;
    pre_mean(gg) = stats.pre_mean;
    play_mean(gg) = stats.play_mean;
    mean_difference(gg) = stats.mean_difference;
    signrank_p(gg) = stats.p_value;
    significant(gg) = stats.significant;
end
summary = table(repmat(string(area_label), n_groups, 1), ...
    string(group_names(:)), n_sessions, n_neurons, pre_mean, play_mean, ...
    mean_difference, signrank_p, significant, sessions, ...
    'VariableNames', names);
end

function stats = paired_group_stats(pre_values, play_values, pair_alpha)
pre_values = pre_values(:);
play_values = play_values(:);
finite = isfinite(pre_values) & isfinite(play_values);
pre_values = pre_values(finite);
play_values = play_values(finite);
stats.n_neurons = numel(pre_values);
if stats.n_neurons == 0
    stats.pre_mean = NaN;
    stats.play_mean = NaN;
    stats.mean_difference = NaN;
    stats.p_value = NaN;
    stats.significant = false;
    return
end
stats.pre_mean = mean(pre_values);
stats.play_mean = mean(play_values);
stats.mean_difference = mean(play_values - pre_values);
stats.p_value = paired_signrank_p(pre_values, play_values);
stats.significant = isfinite(stats.p_value) && stats.p_value < pair_alpha;
end

function rows = paired_level_rows(area_label, pre_z, play_z, ...
    session_summary, animal_summary)
neuron_row = one_level_row(area_label, "Neuron", pre_z, play_z, NaN);
session_row = one_level_row(area_label, "Session", ...
    session_summary.PreMean, session_summary.PlayMean, ...
    within_unit_significant_count(session_summary));
animal_row = one_level_row(area_label, "Animal", ...
    animal_summary.PreMean, animal_summary.PlayMean, ...
    within_unit_significant_count(animal_summary));
rows = [neuron_row; session_row; animal_row];
end

function count = within_unit_significant_count(summary)
if isempty(summary)
    count = NaN;
    return
end
count = sum(summary.Significant);
end

function row = one_level_row(area_label, level_name, pre_values, ...
    play_values, n_within_significant)
pre_values = pre_values(:);
play_values = play_values(:);
finite = isfinite(pre_values) & isfinite(play_values);
pre_values = pre_values(finite);
play_values = play_values(finite);
difference = play_values - pre_values;
row = table(string(area_label), string(level_name), numel(difference), ...
    sum(difference > 0), sum(difference < 0), ...
    paired_signrank_p(pre_values, play_values), n_within_significant, ...
    'VariableNames', {'Area','Level','N','NPositive','NNegative', ...
    'SignrankP','NWithinUnitSignificant'});
end

function row = one_neuron_count_relation(level_name, area_label, summary)
names = {'Level','Area','NUnits','NSignificant','NNotSignificant', ...
    'MedianNeuronsSignificant','MedianNeuronsNotSignificant', ...
    'RanksumP','SpearmanRho','SpearmanP'};
if isempty(summary)
    row = table(string(level_name), string(area_label), 0, 0, 0, ...
        NaN, NaN, NaN, NaN, NaN, 'VariableNames', names);
    return
end
usable = isfinite(summary.NNeurons) & isfinite(summary.MeanDifference);
n_neurons = summary.NNeurons(usable);
significant = summary.Significant(usable);
effect = summary.MeanDifference(usable);
n_significant = sum(significant);
n_not = sum(~significant);
if n_significant > 0
    median_significant = median(n_neurons(significant));
else
    median_significant = NaN;
end
if n_not > 0
    median_not = median(n_neurons(~significant));
else
    median_not = NaN;
end
ranksum_p = NaN;
if n_significant >= 1 && n_not >= 1
    try
        ranksum_p = ranksum(n_neurons(significant), n_neurons(~significant));
    catch
        ranksum_p = NaN;
    end
end
[rho, rho_p] = spearman_pair(n_neurons, effect);
row = table(string(level_name), string(area_label), numel(n_neurons), ...
    n_significant, n_not, median_significant, median_not, ranksum_p, ...
    rho, rho_p, 'VariableNames', names);
end

function [rho, p_value] = spearman_pair(x, y)
rho = NaN;
p_value = NaN;
x = x(:);
y = y(:);
keep = isfinite(x) & isfinite(y);
if nnz(keep) < 3
    return
end
x = x(keep);
y = y(keep);
if numel(unique(x)) < 2 || numel(unique(y)) < 2
    return
end
try
    [rho, p_value] = corr(x, y, 'Type', 'Spearman');
catch
end
end

function plot_paired_unit_figure(figure_name, paper_labels, summary, ...
    animal_names, animal_colors, axis_limit, level_word, pair_alpha)
figure('Name', figure_name, 'Units', 'normalized', ...
    'OuterPosition', [0.02 0.28 0.96 0.42]);
n_areas = numel(paper_labels);
for an = 1:n_areas+1
    subplot(1, n_areas+2, an)
    if an <= n_areas
        area_label = paper_labels(an);
        title_area = char(paper_labels(an));
    else
        area_label = "ALL";
        title_area = 'ALL AREAS';
    end
    in_area = summary.Area == area_label;
    hold on
    plot(axis_limit, axis_limit, 'r-')
    plot_colored_markers(summary.PreMean(in_area), summary.PlayMean(in_area), ...
        summary.AnimalID(in_area), summary.Significant(in_area), ...
        animal_names, animal_colors)
    axis([axis_limit axis_limit])
    axis square
    n_units = nnz(in_area);
    n_positive = sum(summary.MeanDifference(in_area) > 0);
    n_significant = sum(summary.Significant(in_area));
    title({sprintf('%s %s p=%s', title_area, level_word, ...
        format_p_value(paired_signrank_p( ...
        summary.PreMean(in_area), summary.PlayMean(in_area)))), ...
        sprintf('%d/%d > 0, filled %d/%d', n_positive, n_units, ...
        n_significant, n_units)})
    xlabel('Pre-play PPC (atanh)')
    if an == 1
        ylabel('Play PPC (atanh)')
    end
end
subplot(1, n_areas+2, n_areas+2)
draw_pair_color_key(animal_names, animal_colors, pair_alpha)
end

function plot_neuron_count_figure(session_summary, animal_summary, ...
    animal_names, animal_colors, relation, pair_alpha)
figure('Name', 'PPC significance and neuron count', ...
    'Units', 'normalized', 'OuterPosition', [0.08 0.12 0.84 0.46]);
subplot(1, 3, 1)
plot_count_panel(session_summary, animal_names, animal_colors, ...
    relation, "Session", 'Neurons in the session')
subplot(1, 3, 2)
plot_count_panel(animal_summary, animal_names, animal_colors, ...
    relation, "Animal", 'Neurons in the animal')
subplot(1, 3, 3)
draw_pair_color_key(animal_names, animal_colors, pair_alpha)
end

function plot_count_panel(summary, animal_names, animal_colors, ...
    relation, level_name, x_label)
in_area = summary.Area == "ALL";
hold on
yline(0, 'Color', 'r')
plot_colored_markers(summary.NNeurons(in_area), ...
    summary.MeanDifference(in_area), summary.AnimalID(in_area), ...
    summary.Significant(in_area), animal_names, animal_colors)
xlabel(x_label)
ylabel('Mean play - pre PPC (atanh)')
this_relation = relation(relation.Level == level_name & ...
    relation.Area == "ALL", :);
if height(this_relation) == 1
    title({char(level_name), ...
        sprintf('neuron-count ranksum p=%s', ...
        format_p_value(this_relation.RanksumP)), ...
        sprintf('effect Spearman \\rho=%.2f p=%s', ...
        this_relation.SpearmanRho, format_p_value(this_relation.SpearmanP))})
end
end

function plot_colored_markers(x, y, animal_id, is_significant, ...
    animal_names, colors)
x = x(:);
y = y(:);
animal_id = string(animal_id(:));
is_significant = logical(is_significant(:));
for ii = 1:numel(x)
    if ~isfinite(x(ii)) || ~isfinite(y(ii))
        continue
    end
    color = animal_color(animal_id(ii), animal_names, colors);
    if is_significant(ii)
        face = color;
    else
        face = 'none';
    end
    plot(x(ii), y(ii), 'o', ...
        'MarkerEdgeColor', color, ...
        'MarkerFaceColor', face, ...
        'MarkerSize', 8, ...
        'LineWidth', 1.2)
end
end

function draw_pair_color_key(animal_names, animal_colors, pair_alpha)
animal_names = string(animal_names(:));
hold on
key_handles = gobjects(numel(animal_names)+2, 1);
key_handles(1) = plot(nan, nan, 'o', ...
    'MarkerEdgeColor', [0.15 0.15 0.15], ...
    'MarkerFaceColor', [0.15 0.15 0.15], ...
    'MarkerSize', 8, ...
    'DisplayName', sprintf('filled: p < %.2g', pair_alpha));
key_handles(2) = plot(nan, nan, 'o', ...
    'MarkerEdgeColor', [0.15 0.15 0.15], ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize', 8, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'open: not significant');
for aa = 1:numel(animal_names)
    key_handles(aa+2) = plot(nan, nan, 'o', ...
        'MarkerEdgeColor', animal_colors(aa, :), ...
        'MarkerFaceColor', animal_colors(aa, :), ...
        'MarkerSize', 8, ...
        'DisplayName', char(animal_names(aa)));
end
xlim([0 1])
ylim([0 1])
color_key = legend(key_handles, 'Location', 'best', ...
    'FontSize', 7, 'Interpreter', 'none');
set(gca, 'Visible', 'off')
set(color_key, 'Visible', 'on')
end

function color = animal_color(animal_id, animal_names, colors)
index = find(string(animal_names(:)) == string(animal_id), 1);
if isempty(index)
    color = [0.45 0.45 0.45];
    return
end
color = colors(index, :);
end

function colors = paired_animal_colors(n_animals)
if n_animals < 1
    colors = zeros(0, 3);
    return
end
hue = (0:n_animals-1)' / n_animals;
colors = hsv2rgb([hue, repmat(0.75, n_animals, 1), repmat(0.90, n_animals, 1)]);
end

function label = format_p_value(p_value)
if numel(p_value) ~= 1 || ~isfinite(p_value)
    label = 'NA';
    return
end
label = sprintf('%.3g', p_value);
end
