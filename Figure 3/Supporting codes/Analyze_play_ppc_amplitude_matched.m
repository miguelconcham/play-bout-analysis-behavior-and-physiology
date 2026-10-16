%% Play PPC after matching one-second epochs by delta amplitude/power
% This is the focused referee control. It does not estimate a PPC
% significance value and does not include non-play behaviors.
%
% Outputs:
%   matched_summary       one row per neuron
%   matched_draw_table    repeated amplitude-matched arrangements
%   null_draw_table       merged-pool label randomizations under no condition
%   ppc_power_spike_model PPC change ~ power change + spike-count terms

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
saving_folder = fullfile(data_root, 'Analysis results', ...
    'phase locking power control');
load(fullfile(saving_folder, ...
    'phase_locking_power_control_all_sessions.mat'), ...
    'power_control_struct');

%% Hardcoded matching and reliability parameters
EPOCH_SECONDS = 1;
POWER_CALIPER_SD = 0.25;
TARGET_PLAY_FRACTION = 0.75;
N_DRAWS = 500;
MIN_SPIKES_PER_CONDITION = 20;
MIN_RELIABLE_MATCHED_DRAWS = 50;
MAX_ABS_POWER_SMD = 0.10;
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
            if matching_capacity < 2
                excluded_rows(end+1,:) = [identity, ... %#ok<SAGROW>
                    {"Fewer than two power-matchable epochs", matching_capacity, 0}];
                continue
            end

            initial_target_epochs = max(2, ...
                floor(TARGET_PLAY_FRACTION*numel(valid_play)));
            target_epochs = min(initial_target_epochs, matching_capacity);
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
                        target_epochs >= matching_capacity
                    break
                end
                expanded_for_spikes = expanded_for_spikes || ~enough_spikes;
                expanded_for_balance = expanded_for_balance || ...
                    ~enough_power_draws;
                target_epochs = min(matching_capacity, ...
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

            valid_merged = [valid_pre; numel(pre_power)+valid_play];
            [null_pre_index, null_post_index] = merged_null_draws( ...
                valid_merged, target_epochs, N_DRAWS);
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
            cell_type = preferred_phase_type( ...
                neuron.entire_preferred_angle);
            summary_rows(end+1,:) = {animal, session, neuron_key, ... %#ok<SAGROW>
                channel_key, neuron.cluster_id, string(neuron.area), ...
                cell_type, channel, numel(valid_play), ...
                initial_target_epochs, target_epochs, ...
                target_epochs/numel(valid_play), matching_capacity, ...
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
    'MatchingCapacity', ...
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

%% Model PPC changes under the merged-pool no-condition null
% MeanSpikeCount is included because PPC variance depends on the available
% spike count; SpikeCountDifference captures condition imbalance requested
% in the proposed model. EpochCount records adaptive duration. Pre-play and
% play labels are absent from these draws: equal-sized pseudo-conditions are
% sampled from the merged epoch pool. Because draws reuse epochs, use this
% model for conditional prediction; its nominal coefficient p-values are
% not independent-resample inference.
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

analysis_parameters = struct( ...
    'epoch_seconds', EPOCH_SECONDS, ...
    'power_caliper_sd', POWER_CALIPER_SD, ...
    'target_play_fraction', TARGET_PLAY_FRACTION, ...
    'n_draws', N_DRAWS, ...
    'min_spikes_per_condition', MIN_SPIKES_PER_CONDITION, ...
    'min_reliable_matched_draws', MIN_RELIABLE_MATCHED_DRAWS, ...
    'max_abs_power_smd', MAX_ABS_POWER_SMD, ...
    'random_seed', RANDOM_SEED);
save(fullfile(saving_folder, ...
    'play_ppc_amplitude_matched_analysis.mat'), ...
    'matched_summary', 'matched_draw_table', 'null_draw_table', ...
    'excluded_neurons', 'ppc_power_spike_model', ...
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

n_pre = sum(pre_count(pre_index),2);
n_play = sum(play_count(play_index),2);
r_pre = sum(pre_vector(pre_index),2);
r_play = sum(play_vector(play_index),2);
pre_ppc = ppc_from_sufficient_statistics(r_pre, n_pre);
play_ppc = ppc_from_sufficient_statistics(r_play, n_play);
mean_pre_power = mean(pre_power(pre_index),2,'omitnan');
mean_play_power = mean(play_power(play_index),2,'omitnan');
pooled_sd = sqrt((var(pre_power(pre_index),0,2,'omitnan') + ...
    var(play_power(play_index),0,2,'omitnan'))/2);
power_smd = (mean_play_power-mean_pre_power)./pooled_sd;
power_smd(pooled_sd == 0 & mean_play_power == mean_pre_power) = 0;
reliable = isfinite(pre_ppc) & isfinite(play_ppc);
if require_balance
    reliable = reliable & abs(power_smd) <= max_abs_smd;
end

n_draws = size(pre_index,1);
metrics = [(1:n_draws)', repmat(size(pre_index,2),n_draws,1), ...
    mean_pre_power, mean_play_power, mean_play_power-mean_pre_power, ...
    power_smd, n_pre, n_play, n_play-n_pre, (n_play+n_pre)/2, ...
    pre_ppc, play_ppc, play_ppc-pre_ppc, double(reliable)];
end

function metrics = evaluate_null_draws(neuron, pseudo_pre, pseudo_post, ...
    pre_power, play_power)
merged_power = [pre_power(:); play_power(:)];
merged_count = [neuron.pre_epoch_spike_count(:); ...
    neuron.play_epoch_spike_count(:)];
merged_vector = [cellfun(@phase_vector, neuron.pre_epoch_phase); ...
    cellfun(@phase_vector, neuron.play_epoch_phase)];

n_pre = sum(merged_count(pseudo_pre),2);
n_post = sum(merged_count(pseudo_post),2);
r_pre = sum(merged_vector(pseudo_pre),2);
r_post = sum(merged_vector(pseudo_post),2);
pre_ppc = ppc_from_sufficient_statistics(r_pre, n_pre);
post_ppc = ppc_from_sufficient_statistics(r_post, n_post);
mean_pre_power = mean(merged_power(pseudo_pre),2,'omitnan');
mean_post_power = mean(merged_power(pseudo_post),2,'omitnan');
pooled_sd = sqrt((var(merged_power(pseudo_pre),0,2,'omitnan') + ...
    var(merged_power(pseudo_post),0,2,'omitnan'))/2);
power_smd = (mean_post_power-mean_pre_power)./pooled_sd;
power_smd(pooled_sd == 0 & mean_post_power == mean_pre_power) = 0;
reliable = isfinite(pre_ppc) & isfinite(post_ppc);

n_draws = size(pseudo_pre,1);
metrics = [(1:n_draws)', repmat(size(pseudo_pre,2),n_draws,1), ...
    mean_pre_power, mean_post_power, ...
    mean_post_power-mean_pre_power, power_smd, ...
    n_pre, n_post, n_post-n_pre, (n_post+n_pre)/2, ...
    pre_ppc, post_ppc, post_ppc-pre_ppc, double(reliable)];
end

function value = phase_vector(phase)
value = sum(exp(1i*phase(:)));
end

function value = ppc_from_sufficient_statistics(resultant, n)
value = (abs(resultant).^2-n)./(n.*(n-1));
value(n < 2) = NaN;
end

function label = preferred_phase_type(angle_value)
if ~isfinite(angle_value)
    label = "Unknown";
elseif angle_value < -pi/2 || angle_value > pi/2
    label = "Trough";
else
    label = "Peak";
end
end
