%% Analyze PPC changes after controlling for paired delta power
% PPC is the metric called "PLI" in some Figure 3 notes. This script keeps
% the PPC name because that is what compute_ppc_fast actually estimates.
% Legacy bout-caliper analysis. For the focused referee control, use
% Analyze_play_ppc_amplitude_matched.m (1-s power-stratified epochs).

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
saving_folder = fullfile(data_root, 'Analysis results', 'phase locking power control');
load(fullfile(saving_folder, 'phase_locking_power_control_all_sessions.mat'), ...
    'power_control_struct');

%% Hardcoded control parameters (change here for sensitivity analyses)
POWER_MATCH_SD_MULTIPLIER = 1; % requested threshold
MIN_MATCHED_PAIRS = 5;
MIN_SPIKES_PER_CONDITION = 20;
ARTIFACT_ROBUST_Z_LIMIT = 6;
AMPLITUDE_PERCENTILE = 20;
N_EQUAL_COUNT_DRAWS = 200;
N_MONTE_CARLO = 1000;
ENTIRE_RECORDING_ALPHA = 0.01;
LOCKED_CELLS_ONLY = false; % no PPC significance selection in this control
RANDOM_SEED = 73491;
rng(RANDOM_SEED, 'twister');

% Power matching is done in dB. Within each session/channel:
% abs(10*log10(Pplay/Ppre)) <= 1 * SD(10*log10(Pplay)).
% This uses the requested variation in post-onset/play power while keeping
% the difference and threshold in the same units. It is preferable to raw
% power because power is positive, skewed, and electrode-scale dependent.

observed_rows = {};
total_neurons = sum(arrayfun(@(x) numel(x.neurons), power_control_struct));
monte_carlo_blocks = cell(total_neurons,1);

for ss = 1:numel(power_control_struct)
    session_data = power_control_struct(ss);
    animal = string(session_data.animal);
    session = string(session_data.session);
    for nn = 1:numel(session_data.neurons)
        neuron = session_data.neurons(nn);
        if LOCKED_CELLS_ONLY && ...
                ~(isfinite(neuron.entire_ppc_p) && ...
                neuron.entire_ppc_p < ENTIRE_RECORDING_ALPHA)
            continue
        end

        channel_index = find([session_data.channel_data.channel] == ...
            neuron.channel, 1);
        if isempty(channel_index), continue; end
        channel = session_data.channel_data(channel_index);
        pre_power_db = channel.pre_power_db(:);
        play_power_db = channel.play_power_db(:);
        power_difference_db = channel.power_difference_db(:);

        valid_power = isfinite(pre_power_db) & isfinite(play_power_db);
        valid_power = valid_power & ...
            abs(robust_z(pre_power_db)) <= ARTIFACT_ROBUST_Z_LIMIT & ...
            abs(robust_z(play_power_db)) <= ARTIFACT_ROBUST_Z_LIMIT;
        if sum(valid_power) < MIN_MATCHED_PAIRS, continue; end

        post_power_sd = std(play_power_db(valid_power), 0, 'omitnan');
        threshold_db = POWER_MATCH_SD_MULTIPLIER * post_power_sd;
        matched = valid_power & abs(power_difference_db) <= threshold_db;
        matched_index = find(matched);
        if numel(matched_index) < MIN_MATCHED_PAIRS, continue; end

        pre_phase = concatenate_cells(neuron.pre_phase, matched_index);
        play_phase = concatenate_cells(neuron.play_phase, matched_index);
        if numel(pre_phase) < MIN_SPIKES_PER_CONDITION || ...
                numel(play_phase) < MIN_SPIKES_PER_CONDITION
            continue
        end
        pre_ppc = local_ppc(pre_phase);
        play_ppc = local_ppc(play_phase);
        equal_count_difference = mean_equal_count_difference( ...
            pre_phase, play_phase, N_EQUAL_COUNT_DRAWS);

        amplitude_cutoff = channel.amplitude_reference( ...
            max(1, min(101, round(AMPLITUDE_PERCENTILE) + 1)));
        pre_phase_gated = concatenate_gated(neuron.pre_phase, ...
            neuron.pre_amplitude, matched_index, amplitude_cutoff);
        play_phase_gated = concatenate_gated(neuron.play_phase, ...
            neuron.play_amplitude, matched_index, amplitude_cutoff);
        if numel(pre_phase_gated) >= MIN_SPIKES_PER_CONDITION && ...
                numel(play_phase_gated) >= MIN_SPIKES_PER_CONDITION
            amplitude_gated_difference = local_ppc(play_phase_gated) - ...
                local_ppc(pre_phase_gated);
        else
            amplitude_gated_difference = NaN;
        end

        cell_type = preferred_phase_type(neuron.entire_preferred_angle);
        neuron_key = session + "_ch" + string(neuron.channel) + ...
            "_cl" + string(neuron.cluster_id);
        channel_key = session + "_ch" + string(neuron.channel);
        observed_rows(end+1,:) = {animal, session, neuron_key, channel_key, ... %#ok<SAGROW>
            neuron.cluster_id, string(neuron.area), cell_type, neuron.channel, ...
            numel(matched_index), threshold_db, ...
            mean(pre_power_db(matched), 'omitnan'), ...
            mean(play_power_db(matched), 'omitnan'), ...
            mean(power_difference_db(matched), 'omitnan'), ...
            numel(pre_phase), numel(play_phase), ...
            pre_ppc, play_ppc, play_ppc-pre_ppc, ...
            equal_count_difference, amplitude_gated_difference, ...
            NaN, NaN};
        observed_index = size(observed_rows,1);

        valid_index = find(valid_power);
        subset_size = numel(matched_index);
        if subset_size > numel(valid_index), continue; end
        this_monte_carlo = nan(N_MONTE_CARLO,9);
        accepted_samples = 0;
        for rr = 1:N_MONTE_CARLO
            subset = valid_index(randperm(numel(valid_index), subset_size));
            mc_pre_phase = concatenate_cells(neuron.pre_phase, subset);
            mc_play_phase = concatenate_cells(neuron.play_phase, subset);
            if numel(mc_pre_phase) < MIN_SPIKES_PER_CONDITION || ...
                    numel(mc_play_phase) < MIN_SPIKES_PER_CONDITION
                continue
            end
            % Use the same pooled PPC estimator as PPCDifference above.
            % Equal-count PPC remains a separate spike-count sensitivity.
            mc_ppc_difference = local_ppc(mc_play_phase) - ...
                local_ppc(mc_pre_phase);
            accepted_samples = accepted_samples + 1;
            this_monte_carlo(accepted_samples,:) = [ ...
                observed_index, rr, subset_size, ...
                mean(pre_power_db(subset), 'omitnan'), ...
                mean(play_power_db(subset), 'omitnan'), ...
                mean(power_difference_db(subset), 'omitnan'), ...
                numel(mc_pre_phase), numel(mc_play_phase), ...
                mc_ppc_difference];
        end
        monte_carlo_blocks{observed_index} = ...
            this_monte_carlo(1:accepted_samples,:);
    end
end

observed_names = {'Animal','Session','NeuronKey','ChannelKey','ClusterID', ...
    'Area','CellType','Channel','MatchedPairCount','PowerThresholdDB', ...
    'PrePowerDB','PlayPowerDB','PowerDifferenceDB','PreSpikeCount', ...
    'PlaySpikeCount','PrePPC','PlayPPC','PPCDifference', ...
    'EqualCountPPCDifference','AmplitudeGatedPPCDifference', ...
    'ConditionalExpectedPPCDifference','ConditionalPercentile'};
mc_names = {'ObservedRow','Iteration','PairCount','PrePowerDB', ...
    'PlayPowerDB','PowerDifferenceDB','PreSpikeCount','PlaySpikeCount', ...
    'PPCDifference'};
if isempty(observed_rows)
    error(['No neurons passed the power, pair-count, and spike-count criteria. ' ...
        'Inspect the saved sessions or relax the hardcoded sensitivity parameters.']);
end
if ~any(~cellfun(@isempty, monte_carlo_blocks))
    error('No Monte Carlo samples passed the minimum spike-count criterion.');
end
observed_table = cell2table(observed_rows, 'VariableNames', observed_names);
monte_carlo_matrix = vertcat(monte_carlo_blocks{:});
monte_carlo_table = array2table(monte_carlo_matrix, ...
    'VariableNames', mc_names);

%% Conditional Monte Carlo reference for each observed neuron
% Use equal-sized bout combinations with the nearest mean power differences
% from the same neuron. The percentile is a descriptive conditional rank,
% not a permutation p-value: these resamples preserve the observed labels.
for ii = 1:height(observed_table)
    same_neuron = monte_carlo_table.ObservedRow == ii;
    this_mc = monte_carlo_table(same_neuron,:);
    if height(this_mc) < 20, continue; end
    [~, order] = sort(abs(this_mc.PowerDifferenceDB - ...
        observed_table.PowerDifferenceDB(ii)));
    n_nearest = min(height(this_mc), max(50, ceil(0.05*height(this_mc))));
    reference = this_mc.PPCDifference(order(1:n_nearest));
    expected = median(reference, 'omitnan');
    observed_table.ConditionalExpectedPPCDifference(ii) = expected;
    observed_table.ConditionalPercentile(ii) = ...
        (1 + sum(reference <= observed_table.PPCDifference(ii))) / ...
        (numel(reference) + 1);
end

%% Animal-level summaries and paired tests
% Neurons are not treated as independent replicates for the area/type
% summaries. First average within animal, area, and cell type.
[summary_group, summary_animal, summary_area, summary_type] = findgroups( ...
    observed_table.Animal, observed_table.Area, observed_table.CellType);
animal_mean_ppc = splitapply(@(x) mean(x, 'omitnan'), ...
    observed_table.EqualCountPPCDifference, summary_group);
animal_mean_power = splitapply(@(x) mean(x, 'omitnan'), ...
    observed_table.PowerDifferenceDB, summary_group);
animal_summary = table(summary_animal, summary_area, summary_type, ...
    animal_mean_power, animal_mean_ppc, ...
    'VariableNames', {'Animal','Area','CellType', ...
    'PowerDifferenceDB','PPCDifference'});

[test_group, test_area, test_type] = findgroups( ...
    animal_summary.Area, animal_summary.CellType);
n_groups = max(test_group);
group_statistics = table('Size', [n_groups 6], ...
    'VariableTypes', {'string','string','double','double','double','double'}, ...
    'VariableNames', {'Area','CellType','NAnimals','MedianPPCDifference', ...
    'SignrankP','FDRP'});
for gg = 1:n_groups
    values = animal_summary.PPCDifference(test_group == gg);
    group_statistics.Area(gg) = test_area(gg);
    group_statistics.CellType(gg) = test_type(gg);
    group_statistics.NAnimals(gg) = sum(isfinite(values));
    group_statistics.MedianPPCDifference(gg) = median(values, 'omitnan');
    if sum(isfinite(values)) >= 3
        try
            group_statistics.SignrankP(gg) = signrank(values);
        catch
            group_statistics.SignrankP(gg) = NaN;
        end
    else
        group_statistics.SignrankP(gg) = NaN;
    end
end
group_statistics.FDRP = bh_adjust(group_statistics.SignrankP);

%% Mixed-effects sensitivity model
% The coefficient of PowerDifferenceDB directly tests whether PPC change
% covaries with residual power difference. Random channel effects account
% for neurons sharing the same LFP power measurements.
observed_lme = [];
if height(observed_table) >= 20 && numel(unique(observed_table.Animal)) >= 3
    model_table = observed_table;
    model_table.SpikeCountDifference = ...
        model_table.PlaySpikeCount - model_table.PreSpikeCount;
    model_table.Animal = categorical(model_table.Animal);
    model_table.ChannelKey = categorical(model_table.ChannelKey);
    model_table.Area = categorical(model_table.Area);
    model_table.CellType = categorical(model_table.CellType);
    try
        observed_lme = fitlme(model_table, ...
            ['EqualCountPPCDifference ~ PowerDifferenceDB + PrePowerDB + ' ...
             'SpikeCountDifference + Area + CellType + (1|Animal) + ' ...
             '(1|ChannelKey)']);
        disp(observed_lme)
    catch exception
        warning('Mixed-effects model did not converge: %s', exception.message);
    end
end

%% Diagnostic plots
figure('Name', 'PPC delta-power controls');
subplot(1,3,1)
histogram(observed_table.PowerDifferenceDB, 30)
xline(0, 'r-');
xlabel('Matched play - pre-play delta power (dB)')
ylabel('Neurons')

subplot(1,3,2)
gscatter(observed_table.PowerDifferenceDB, ...
    observed_table.EqualCountPPCDifference, observed_table.CellType)
xline(0, 'k:'); yline(0, 'k:');
xlabel('\Delta delta power (dB)')
ylabel('\Delta PPC (equal spike count)')

subplot(1,3,3)
plot(observed_table.EqualCountPPCDifference, ...
    observed_table.AmplitudeGatedPPCDifference, 'k.')
hold on
x_limits = xlim;
y_limits = ylim;
limits = [min(x_limits(1), y_limits(1)), max(x_limits(2), y_limits(2))];
plot(limits, limits, 'r-')
axis([limits limits]); axis square
xlabel('Power-matched \Delta PPC')
ylabel(sprintf('\\Delta PPC above amplitude P%d', AMPLITUDE_PERCENTILE))

analysis_parameters = struct( ...
    'power_match_sd_multiplier', POWER_MATCH_SD_MULTIPLIER, ...
    'min_matched_pairs', MIN_MATCHED_PAIRS, ...
    'min_spikes_per_condition', MIN_SPIKES_PER_CONDITION, ...
    'artifact_robust_z_limit', ARTIFACT_ROBUST_Z_LIMIT, ...
    'amplitude_percentile', AMPLITUDE_PERCENTILE, ...
    'n_equal_count_draws', N_EQUAL_COUNT_DRAWS, ...
    'n_monte_carlo', N_MONTE_CARLO, ...
    'locked_cells_only', LOCKED_CELLS_ONLY, ...
    'random_seed', RANDOM_SEED);
save(fullfile(saving_folder, 'phase_locking_power_control_analysis.mat'), ...
    'observed_table', 'monte_carlo_table', 'animal_summary', ...
    'group_statistics', 'observed_lme', 'analysis_parameters', '-v7.3');

%% Local functions
function z = robust_z(x)
x = x(:);
center = median(x, 'omitnan');
scale = 1.4826 * median(abs(x-center), 'omitnan');
if ~isfinite(scale) || scale == 0
    z = zeros(size(x));
    z(isfinite(x) & x ~= center) = Inf;
else
    z = (x-center) / scale;
end
z(~isfinite(x)) = NaN;
end

function phase = concatenate_cells(phase_cells, index)
phase = [];
for ii = reshape(index,1,[])
    phase = [phase; phase_cells{ii}(:)]; %#ok<AGROW>
end
phase = phase(isfinite(phase));
end

function phase = concatenate_gated(phase_cells, amplitude_cells, index, cutoff)
phase = [];
for ii = reshape(index,1,[])
    this_phase = phase_cells{ii}(:);
    this_amplitude = amplitude_cells{ii}(:);
    keep = isfinite(this_phase) & isfinite(this_amplitude) & ...
        this_amplitude >= cutoff;
    phase = [phase; this_phase(keep)]; %#ok<AGROW>
end
end

function value = local_ppc(phase)
phase = phase(isfinite(phase));
n = numel(phase);
if n < 2
    value = NaN;
    return
end
resultant = sum(exp(1i*phase));
value = (abs(resultant)^2-n) / (n*(n-1));
end

function difference = one_equal_count_difference(pre_phase, play_phase)
n = min(numel(pre_phase), numel(play_phase));
pre_sample = pre_phase(randperm(numel(pre_phase), n));
play_sample = play_phase(randperm(numel(play_phase), n));
difference = local_ppc(play_sample) - local_ppc(pre_sample);
end

function difference = mean_equal_count_difference(pre_phase, play_phase, n_draws)
values = nan(n_draws,1);
for rr = 1:n_draws
    values(rr) = one_equal_count_difference(pre_phase, play_phase);
end
difference = mean(values, 'omitnan');
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

function adjusted = bh_adjust(p)
adjusted = nan(size(p));
valid = find(isfinite(p));
if isempty(valid), return; end
[sorted_p, order] = sort(p(valid));
m = numel(sorted_p);
sorted_adjusted = sorted_p .* m ./ (1:m)';
sorted_adjusted = flipud(cummin(flipud(sorted_adjusted)));
sorted_adjusted = min(sorted_adjusted, 1);
adjusted(valid(order)) = sorted_adjusted;
end
