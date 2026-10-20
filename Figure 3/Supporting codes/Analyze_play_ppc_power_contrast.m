%% Play PPC after matching one-second epochs by delta amplitude/power
% Power-contrast null using three near-extreme non-play tails. This script
% reads and writes only Analysis results/phase locking power contrast.
% Analyze_play_ppc_amplitude_matched.m and its saved output are left unchanged.
% It does not estimate a PPC significance value.
%
% Outputs:
%   matched_summary       one row per neuron
%   matched_draw_table    repeated amplitude-matched arrangements
%   null_draw_table       play-free power-tail arrangements
%   neuron_model_table    reliable null rows used to fit the model
%   ppc_power_spike_model PPC change ~ power change + spike-count terms

data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
saving_folder = fullfile(data_root, 'Analysis results', ...
    'phase locking power contrast');
session_files = dir(fullfile(saving_folder, '*_phase_locking_power_contrast.mat'));
if isempty(session_files)
    error('No per-session power-control files in %s.', saving_folder);
end
power_control_struct = struct([]);
for ii = 1:numel(session_files)
    loaded_session = load(fullfile(session_files(ii).folder, ...
        session_files(ii).name), 'this_control');
    if ~isfield(loaded_session.this_control, 'nonplay_epoch_table')
        error(['%s does not contain the session non-play pool. ' ...
            'Re-run Estimate_phase_locking_power_contrast.m.'], ...
            session_files(ii).name);
    end
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
N_TAIL_WIDTHS = 3;
TAIL_WIDTH_FACTOR = 2;
MIN_SPIKES_PER_CONDITION = 20;
MIN_RELIABLE_MATCHED_DRAWS = 50;
MAX_ABS_POWER_SMD = 0.10;
LOCKING_ALPHA = 0.01;
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
    if ~isfield(session_data, 'nonplay_epoch_table') || ...
            ~isfield(session_data.neurons, 'nonplay_epoch_phase')
        error(['Saved structures do not contain session-wide non-play epochs. ' ...
            'Re-run Estimate_phase_locking_power_contrast.m; older files are rebuilt automatically.']);
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
        nonplay_power = channel_data.nonplay_epoch_power_db(:);
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

            % Both null groups come from session non-play epochs. Half the
            % pool is the largest group size that keeps those groups disjoint.
            valid_nonplay = find(isfinite(nonplay_power));
            nonplay_null_capacity = floor(numel(valid_nonplay)/2);
            if nonplay_null_capacity < 2
                excluded_rows(end+1,:) = [identity, ... %#ok<SAGROW>
                    {"Fewer than four session non-play epochs for a play-free null", ...
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

            [null_pre_index, null_post_index, pool_width] = ...
                power_contrast_draws(valid_nonplay, nonplay_power, ...
                target_epochs, N_DRAWS, N_TAIL_WIDTHS, TAIL_WIDTH_FACTOR);
            null_metrics = evaluate_null_draws(neuron, ...
                null_pre_index, null_post_index, nonplay_power);
            null_metrics(:,end+1) = pool_width;

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
    'MeanSpikeCount','PseudoPrePPC','PseudoPostPPC','PPCDifference', ...
    'Reliable','PowerPoolWidth'};
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

%% Model PPC changes under the play-free power-contrast null
% MeanSpikeCount is included because PPC variance depends on the available
% spike count; SpikeCountDifference captures condition imbalance requested
% in the proposed model. EpochCount records adaptive duration. Both
% pseudo-conditions are sampled from session non-play epochs, with one group
% drawn from a low-power tail and the other from the corresponding high-power
% tail. Three widths stay within twice K, and draws are split equally across
% widths that can supply distinct groups. Because draws reuse epochs, use
% this model for conditional prediction;
% its nominal coefficient p-values are not independent-resample inference.
ppc_power_spike_model = [];
model_rows = null_draw_table.Reliable == 1;
model_table = null_draw_table(model_rows,:);
neuron_model_table = model_table;
if ~isempty(model_table)
    neuron_model_table.Animal = matched_summary.Animal(model_table.NeuronRow);
    neuron_model_table.Session = matched_summary.Session(model_table.NeuronRow);
    neuron_model_table.NeuronKey = string(matched_summary.NeuronKey(model_table.NeuronRow));
    neuron_model_table.ClusterID = matched_summary.ClusterID(model_table.NeuronRow);
    neuron_model_table.Area = matched_summary.Area(model_table.NeuronRow);
    neuron_model_table.Channel = matched_summary.Channel(model_table.NeuronRow);
end
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
legend({'Non-play power-contrast null','Amplitude-matched estimate'})

%% Power-balanced pre-play versus play PPC, Figure 3e areas
% SC is stored as SupCol. This is the same comparison as the power-control plot.
paper_areas = ["SupCol","DLPAG","LPAG","VLPAG","DR"];
paper_labels = ["SC","DLPAG","LPAG","VLPAG","DR"];
finite_level = isfinite(matched_summary.PrePPC) & isfinite(matched_summary.PlayPPC);
level_rows = matched_summary(finite_level & ...
    ismember(string(matched_summary.Area), paper_areas), :);
pre_z = atanh(min(max(level_rows.PrePPC, -0.999), 0.999));
play_z = atanh(min(max(level_rows.PlayPPC, -0.999), 0.999));
level_areas = string(level_rows.Area);
fprintf('Areas included in the power-balanced all-areas panel:\n');
disp(unique(level_areas));
level_limit = [min([-0.005; pre_z; play_z]), max([0.07; pre_z; play_z])];
figure('Name', 'Power-balanced trough PPC, pre-play vs play');
for an = 1:numel(paper_areas)
    subplot(1, numel(paper_areas)+1, an)
    in_area = level_areas == paper_areas(an);
    plot_paired_ppc(pre_z(in_area), play_z(in_area), level_limit)
    title(sprintf('%s n=%d p=%s', paper_labels(an), sum(in_area), ...
        paired_signrank_text(pre_z(in_area), play_z(in_area))))
    xlabel('Pre-play PPC (atanh)')
    if an == 1
        ylabel('Play PPC (atanh)')
    end
end
subplot(1, numel(paper_areas)+1, numel(paper_areas)+1)
plot_paired_ppc(pre_z, play_z, level_limit)
xlabel('Pre-play PPC (atanh)')
title(sprintf('ALL AREAS n=%d neurons, %d animals, p=%s', ...
    height(level_rows), numel(unique(string(level_rows.Animal))), ...
    paired_signrank_text(pre_z, play_z)))

%% Observed versus null-model PPC difference, Figure 3e areas
% The prediction uses each neuron's matched power difference, spike-count
% difference, mean spike count, and epoch count.
finite_difference = isfinite(matched_summary.PPCDifference) & ...
    isfinite(matched_summary.ModelExpectedPPCDifference);
plot_rows = matched_summary(finite_difference & ...
    ismember(string(matched_summary.Area), paper_areas), :);
predicted_difference = plot_rows.ModelExpectedPPCDifference;
observed_difference = plot_rows.PPCDifference;
present_areas = string(plot_rows.Area);
fprintf('Areas included in the PPC-difference all-areas panel:\n');
disp(unique(present_areas));
span = [predicted_difference; observed_difference; 0];
pad = 0.05 * max(max(span) - min(span), eps);
axis_limit = [min(span) - pad, max(span) + pad];
figure('Name', 'Observed versus null-model PPC difference');
for an = 1:numel(paper_areas)
    subplot(1, numel(paper_areas)+1, an)
    in_area = present_areas == paper_areas(an);
    plot_paired_ppc(predicted_difference(in_area), ...
        observed_difference(in_area), axis_limit)
    title(sprintf('%s n=%d p=%s', paper_labels(an), sum(in_area), ...
        paired_signrank_text(predicted_difference(in_area), ...
        observed_difference(in_area))))
    xlabel('Null-model PPC difference')
    if an == 1
        ylabel('Observed PPC difference')
    end
end
subplot(1, numel(paper_areas)+1, numel(paper_areas)+1)
plot_paired_ppc(predicted_difference, observed_difference, axis_limit)
xlabel('Null-model PPC difference')
title(sprintf('ALL AREAS n=%d neurons, %d animals, p=%s', ...
    height(plot_rows), numel(unique(string(plot_rows.Animal))), ...
    paired_signrank_text(predicted_difference, observed_difference)))

%% Example matched powers and one neuron's null-model fit
% Estimate stores every epoch power, not the matched pairs. The pairs below
% are rebuilt from those saved values. LFP and spike phases are not estimated again.
if ~isempty(plot_rows)
    null_per_neuron = accumarray(neuron_model_table.NeuronRow, 1, ...
        [height(matched_summary), 1]);
    candidate_rows = find(finite_difference & ...
        ismember(string(matched_summary.Area), paper_areas));
    [~, richest] = max(null_per_neuron(candidate_rows));
    example_row = candidate_rows(richest);
    example_neuron = matched_summary(example_row,:);
    session_index = find(string({power_control_struct.session}) == ...
        string(example_neuron.Session), 1);
    session_data = power_control_struct(session_index);
    channel_index = find([session_data.channel_data.channel] == ...
        example_neuron.Channel, 1);
    channel_data = session_data.channel_data(channel_index);
    [matched_pairs, example_caliper] = nearest_power_pairs( ...
        channel_data.pre_epoch_power_db(:), ...
        channel_data.play_epoch_power_db(:), POWER_CALIPER_SD);
    if isempty(matched_pairs)
        warning('Example neuron has no power-matched epoch pairs.');
        matched_pre_power = [];
        matched_play_power = [];
    else
        matched_pre_power = channel_data.pre_epoch_power_db(matched_pairs(:,1));
        matched_play_power = channel_data.play_epoch_power_db(matched_pairs(:,2));
        matched_pre_power = matched_pre_power(:);
        matched_play_power = matched_play_power(:);
    end
    fprintf(['Example matched powers: %s, channel %d, cluster %d, ' ...
        '%d pairs, caliper %.3g dB.\n'], string(example_neuron.Session), ...
        example_neuron.Channel, example_neuron.ClusterID, ...
        numel(matched_pre_power), example_caliper);

    if ~isempty(matched_pre_power)
        figure('Name', 'Example matched delta-power distributions');
        subplot(1,2,1)
        histogram(matched_pre_power, 'Normalization', 'probability', ...
            'FaceAlpha', 0.5)
        hold on
        histogram(matched_play_power, 'Normalization', 'probability', ...
            'FaceAlpha', 0.5)
        xlabel('Delta power (dB)')
        ylabel('Probability')
        legend({'Matched pre-play','Matched play'})
        title(sprintf('n=%d pairs', numel(matched_pre_power)))
        subplot(1,2,2)
        plot(matched_pre_power, matched_play_power, 'k.')
        hold on
        pair_limit = [min([matched_pre_power; matched_play_power]), ...
            max([matched_pre_power; matched_play_power])];
        plot(pair_limit, pair_limit, 'r-')
        axis square
        xlabel('Matched pre-play power (dB)')
        ylabel('Matched play power (dB)')
    end

    if ~isempty(ppc_power_spike_model)
        example_null = neuron_model_table( ...
            neuron_model_table.NeuronRow == example_row, :);
        coefficient = fixedEffects(ppc_power_spike_model);
        null_predictor = fixed_effect_predictor(example_null, coefficient);
        observed_features = table(example_neuron.PowerDifferenceDB, ...
            example_neuron.SpikeCountDifference, ...
            example_neuron.MeanSpikeCount, ...
            example_neuron.EpochsPerCondition, 'VariableNames', ...
            {'PowerDifferenceDB','SpikeCountDifference', ...
            'MeanSpikeCount','EpochCount'});
        observed_predictor = fixed_effect_predictor( ...
            observed_features, coefficient);
        predictor_span = [min([null_predictor; observed_predictor]), ...
            max([null_predictor; observed_predictor])];
        figure('Name', 'Example null-model fit');
        plot(null_predictor, example_null.PPCDifference, '.', ...
            'Color', [0.65 0.65 0.65])
        hold on
        plot(predictor_span, predictor_span, 'k-')
        plot(observed_predictor, example_neuron.ModelExpectedPPCDifference, ...
            'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 9)
        plot(observed_predictor, example_neuron.PPCDifference, ...
            'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9)
        xlabel('Fixed-effect linear predictor')
        ylabel('PPC difference')
        legend({'Null draws','Fixed-effect model','Predicted','Observed'}, ...
            'Location', 'best')
        title(sprintf('%s ch%d cl%d, n=%d null draws', ...
            string(example_neuron.NeuronKey), example_neuron.Channel, ...
            example_neuron.ClusterID, height(example_null)))
    end
end

analysis_parameters = struct( ...
    'epoch_seconds', EPOCH_SECONDS, ...
    'power_caliper_sd', POWER_CALIPER_SD, ...
    'target_play_fraction', TARGET_PLAY_FRACTION, ...
    'n_draws', N_DRAWS, ...
    'n_tail_widths', N_TAIL_WIDTHS, ...
    'tail_width_factor', TAIL_WIDTH_FACTOR, ...
    'min_spikes_per_condition', MIN_SPIKES_PER_CONDITION, ...
    'min_reliable_matched_draws', MIN_RELIABLE_MATCHED_DRAWS, ...
    'max_abs_power_smd', MAX_ABS_POWER_SMD, ...
    'locking_alpha_from_original_analysis', LOCKING_ALPHA, ...
    'n_animal_bootstrap', N_ANIMAL_BOOTSTRAP, ...
    'n_sign_flip', N_SIGN_FLIP, ...
    'random_seed', RANDOM_SEED);
save(fullfile(saving_folder, ...
    'play_ppc_power_contrast_analysis.mat'), ...
    'matched_summary', 'matched_draw_table', 'null_draw_table', ...
    'neuron_model_table', 'excluded_neurons', 'ppc_power_spike_model', ...
    'animal_summary', 'overall_result', 'analysis_parameters', '-v7.3');

%% Local functions
function predictor = fixed_effect_predictor(data_table, coefficient)
predictor = coefficient(1) + ...
    coefficient(2) * data_table.PowerDifferenceDB + ...
    coefficient(3) * data_table.SpikeCountDifference + ...
    coefficient(4) * data_table.MeanSpikeCount + ...
    coefficient(5) * data_table.EpochCount;
end

function label = paired_signrank_text(x_value, y_value)
finite_pairs = isfinite(x_value) & isfinite(y_value);
if nnz(finite_pairs) < 1
    label = "NA";
    return
end
try
    p_value = signrank(x_value(finite_pairs), y_value(finite_pairs));
catch
    label = "NA";
    return
end
label = sprintf('%.3g', p_value);
end

function plot_paired_ppc(x_value, y_value, axis_limit)
plot(x_value, y_value, 'k.')
hold on
finite_pairs = isfinite(x_value) & isfinite(y_value);
if any(finite_pairs)
    pair_median = [median(x_value(finite_pairs)), median(y_value(finite_pairs))];
    plot(pair_median(1), pair_median(2), 'rx', 'LineWidth', 1.2)
end
plot(axis_limit, axis_limit, 'r-')
axis([axis_limit axis_limit])
axis square
end

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

function [pseudo_low, pseudo_high, pool_width] = power_contrast_draws( ...
    epoch_index, power_db, k, n_requested, n_widths, width_factor)
finite = isfinite(power_db(epoch_index));
epoch_index = epoch_index(finite);
[~, order] = sort(power_db(epoch_index), 'ascend');
sorted_index = epoch_index(order);
sorted_power = power_db(sorted_index);
n_epochs = numel(sorted_index);
half = floor(n_epochs/2);
if k > half
    error('Cannot form two disjoint groups of %d from %d epochs.', k, n_epochs);
end

% Keep every pool inside the extreme tails. The widest pool is twice K,
% or half the epochs when the session is shorter than that. Widths are
% placed so the low-versus-high power gap moves in equal steps. If that
% gap does not change, the widths are spaced by pool size so the sampled
% epochs, and therefore the spike counts, can still differ. Each width that
% can supply distinct groups receives the same number of draws.
widest = min(half, max(k, round(width_factor * k)));
sorted_power = sorted_power(:);
candidate = (k:widest).';
low_sum = cumsum(sorted_power);
high_sum = cumsum(sorted_power(end:-1:1));
power_gap = high_sum(candidate) ./ candidate - low_sum(candidate) ./ candidate;
gap_targets = linspace(power_gap(1), power_gap(end), n_widths);
widths = zeros(1, n_widths);
for ii = 1:n_widths
    [~, pick] = min(abs(power_gap - gap_targets(ii)));
    widths(ii) = candidate(pick);
end
widths = unique(widths, 'stable');
gap_span = power_gap(1) - power_gap(end);
if ~(gap_span > max(1e-6, eps(power_gap(1))))
    widths = unique(round(linspace(k, widest, n_widths)), 'stable');
end
n_base = max(1, floor(n_requested/2));
quotas = equal_width_quotas(widths, k, n_base);

low_rows = zeros(0,k);
high_rows = zeros(0,k);
width_rows = zeros(0,1);
for ii = 1:numel(widths)
    if quotas(ii) < 1
        continue
    end
    low_pool = sorted_index(1:widths(ii));
    high_pool = sorted_index(n_epochs-widths(ii)+1:n_epochs);
    [these_low, these_high] = sample_two_pools( ...
        low_pool, high_pool, k, quotas(ii));
    low_rows = [low_rows; these_low]; %#ok<AGROW>
    high_rows = [high_rows; these_high]; %#ok<AGROW>
    width_rows = [width_rows; repmat(widths(ii), size(these_low,1), 1)]; %#ok<AGROW>
end
pseudo_low = [low_rows; high_rows];
pseudo_high = [high_rows; low_rows];
pool_width = [width_rows; width_rows];
end

function quotas = equal_width_quotas(widths, k, n_base)
n_widths = numel(widths);
capacity = repmat(n_base, n_widths, 1);
for ii = 1:n_widths
    log_pairs = 2 * (gammaln(widths(ii)+1) - gammaln(k+1) - ...
        gammaln(widths(ii)-k+1));
    if log_pairs <= log(n_base)
        capacity(ii) = max(1, round(exp(log_pairs)));
    end
end
quotas = zeros(n_widths, 1);
open = true(n_widths, 1);
left = n_base;
while left > 0 && any(open)
    open_idx = find(open);
    n_open = numel(open_idx);
    share = floor(left / n_open);
    extra = mod(left, n_open);
    give = zeros(n_widths, 1);
    give(open_idx) = share;
    if extra > 0
        give(open_idx(1:extra)) = give(open_idx(1:extra)) + 1;
    end
    add = min(give, capacity - quotas);
    add(~open) = 0;
    if sum(add) < 1
        break
    end
    quotas = quotas + add;
    left = left - sum(add);
    open = quotas < capacity;
end
end

function [low_sets, high_sets] = sample_two_pools(low_pool, high_pool, k, quota)
log_one = gammaln(numel(low_pool)+1) - gammaln(k+1) - ...
    gammaln(numel(low_pool)-k+1);
if log_one <= log(2000)
    low_combinations = nchoosek(low_pool(:)', k);
    high_combinations = nchoosek(high_pool(:)', k);
    n_pairs = size(low_combinations,1) * size(high_combinations,1);
    if n_pairs <= quota
        [low_row, high_row] = ndgrid(1:size(low_combinations,1), ...
            1:size(high_combinations,1));
        low_sets = low_combinations(low_row(:),:);
        high_sets = high_combinations(high_row(:),:);
        return
    end
    selected = randperm(n_pairs, quota);
    [low_row, high_row] = ind2sub([size(low_combinations,1), ...
        size(high_combinations,1)], selected);
    low_sets = low_combinations(low_row(:),:);
    high_sets = high_combinations(high_row(:),:);
    return
end

low_sets = zeros(0,k);
high_sets = zeros(0,k);
seen = zeros(0,2*k);
guard = 0;
while size(low_sets,1) < quota && guard < max(100, 50*quota)
    guard = guard + 1;
    this_low = sort(low_pool(randperm(numel(low_pool), k)));
    this_high = sort(high_pool(randperm(numel(high_pool), k)));
    key = [this_low(:)', this_high(:)'];
    if isempty(seen) || ~ismember(key, seen, 'rows')
        seen = [seen; key]; %#ok<AGROW>
        low_sets = [low_sets; this_low(:)']; %#ok<AGROW>
        high_sets = [high_sets; this_high(:)']; %#ok<AGROW>
    end
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
    nonplay_power)
merged_power = nonplay_power(:);
merged_count = neuron.nonplay_epoch_spike_count(:);
merged_vector = cellfun(@phase_vector, neuron.nonplay_epoch_phase);

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
