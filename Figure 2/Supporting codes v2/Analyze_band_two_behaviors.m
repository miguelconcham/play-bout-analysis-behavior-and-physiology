%% Analyze_band_two_behaviors
% Compare delta, theta, or gamma power for two selectable behaviors.
% Defaults: gamma, Pounce/PounceI.
%
% The input is NOT averaged per session. In
% GENERATE_PSTH_ALL_BEHAVIORS_BANDS.m, every row of
% all_behavior_psth_self/other{alignment, behavior} is one annotated
% behavior event. psth_structure and animal_names retain its recording
% session and electrode.
%
% This script therefore fits the event-level model separately at each time:
%
%   BandPower ~ 1 + Behavior + (1|Animal) + (1|Session)
%
% Behavior is the fixed effect. Animal and recording day are random effects,
% with recording days uniquely nested in animals. Only sessions containing
% both selected behaviors are included.
% This multilevel approach retains the event-level information while
% accounting for events nested in recording sessions (Aarts et al., 2014,
% Nature Neuroscience, doi:10.1038/nn.3648; Yu et al., 2022, eNeuro,
% doi:10.1523/ENEURO.0330-21.2021).
%
% A separate LME is needed at each time point. Benjamini-Hochberg FDR is
% applied across time points within each alignment; the gray regions in the
% figures show q < alpha. This gives pointwise error control and avoids
% interpreting a cluster p-value as evidence for every time point (König
% et al., 2024, NeuroImage, doi:10.1016/j.neuroimage.2024.120557).
%
% Two figures are produced: onset-aligned and time-wrapped. Each includes
% every individual event trace, model-estimated means with 95% CIs, every
% paired session-mean difference trace, and the LME behavior contrast.

%% Settings
band_name = 'gamma';                % 'delta', 'theta', or 'gamma'
behavior_names = {'Pounce', 'PounceI'};
behavior_role = 'self';             % 'self' or 'other'
electrode_to_use = 1;
normalization = 'zscore';           % 'zscore' (existing analysis) or 'none'
alpha = 0.05;
session_name_filters = {'B1D1', 'B1S3', 'B2S2', 'B3D2', 'B4S2', 'B4D4'};

data_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';
output_folder = fullfile(data_folder, 'two behavior band comparison');

all_behaviors = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', ...
    'CB', 'CD', 'Grooming', 'Rearing', 'Sniffing', 'Scratching', ...
    'PounceI', 'Bite'};

%% Validate settings and load the event-level PSTHs
valid_bands = {'delta', 'theta', 'gamma'};
if ~ismember(lower(band_name), valid_bands)
    error('Analyze_band_two_behaviors:Band', ...
        'band_name must be ''delta'', ''theta'', or ''gamma''.');
end
band_name = lower(band_name);
structure_file = sprintf('psth_structure_%s_all_behaviors_self_other_partner.mat', band_name);
names_file = sprintf('animal_names_%s_all_behaviors_self_other_partner.mat', band_name);

if numel(behavior_names) ~= 2 || strcmpi(behavior_names{1}, behavior_names{2})
    error('Analyze_band_two_behaviors:Behaviors', ...
        'behavior_names must contain two different behavior names.');
end

behavior_idx = nan(1, 2);
for b = 1:2
    this_idx = find(strcmpi(all_behaviors, behavior_names{b}), 1);
    if isempty(this_idx)
        error('Analyze_band_two_behaviors:BehaviorName', ...
            'Unknown behavior "%s". Choose from: %s.', ...
            behavior_names{b}, strjoin(all_behaviors, ', '));
    end
    behavior_idx(b) = this_idx;
    behavior_names{b} = all_behaviors{behavior_idx(b)};
end

if ~ismember(lower(behavior_role), {'self', 'other'})
    error('Analyze_band_two_behaviors:Role', ...
        'behavior_role must be ''self'' or ''other''.');
end
if ~ismember(lower(normalization), {'zscore', 'none'})
    error('Analyze_band_two_behaviors:Normalization', ...
        'normalization must be ''zscore'' or ''none''.');
end

fprintf('Loading individual-event %s PSTHs\n', band_name)
loaded_psth = load(fullfile(data_folder, structure_file), 'psth_structure');
loaded_names = load(fullfile(data_folder, names_file), 'animal_names');
psth_structure = loaded_psth.psth_structure;
animal_names = loaded_names.animal_names;

if numel(psth_structure) ~= size(animal_names, 1)
    error('Analyze_band_two_behaviors:Names', ...
        'animal_names has %d rows but psth_structure has %d entries.', ...
        size(animal_names, 1), numel(psth_structure));
end
if isempty(psth_structure)
    error('Analyze_band_two_behaviors:Empty', 'The loaded PSTH structure is empty.');
end
if exist(output_folder, 'dir') ~= 7
    mkdir(output_folder)
end

bin_size = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
if ~(isfinite(bin_size) && bin_size > 0)
    error('Analyze_band_two_behaviors:BinSize', 'The saved spectrogram bin size is invalid.');
end

alignments = struct( ...
    'name', {'onset', 'time_warped'}, ...
    'cell_row', {1, 3}, ...
    'baseline_range', {[-2 0], [-5 0]}, ...
    'analysis_range', {[-1 2], [-5 10]});

results = struct;
results.band_name = band_name;
results.behavior_names = behavior_names;
results.behavior_role = behavior_role;
results.electrode = electrode_to_use;
results.normalization = normalization;
results.alpha = alpha;
results.formula = 'BandPower ~ 1 + Behavior + (1|Animal) + (1|Session)';
results.multiple_comparison = 'Benjamini-Hochberg FDR across time within alignment';
results.input_level = 'one row per annotated behavior event; not session-averaged';

%% Analyze onset-aligned and time-wrapped PSTHs
for a = 1:numel(alignments)
    alignment = alignments(a);
    [full_time, expected_n] = alignment_time(psth_structure(1), alignment.name, bin_size);

    traces = cell(1, 2);
    session_ids = cell(1, 2);
    for b = 1:2
        [traces{b}, session_ids{b}] = collect_behavior_events( ...
            psth_structure, animal_names, alignment.cell_row, ...
            behavior_idx(b), behavior_role, electrode_to_use, ...
            session_name_filters, expected_n);
        [traces{b}, session_ids{b}] = normalize_events( ...
            traces{b}, session_ids{b}, full_time, ...
            alignment.baseline_range, normalization);
    end

    % A repeated-measures contrast requires both behaviors in each session.
    paired_sessions = intersect(unique(session_ids{1}), unique(session_ids{2}));
    if numel(paired_sessions) < 2
        error('Analyze_band_two_behaviors:PairedSessions', ...
            ['%s has only %d session(s) containing both %s and %s. ' ...
             'At least two are required for the mixed model.'], ...
            alignment.name, numel(paired_sessions), ...
            behavior_names{1}, behavior_names{2});
    end
    paired_animals = unique(animal_from_sessions(paired_sessions));
    if numel(paired_animals) < 2
        error('Analyze_band_two_behaviors:PairedAnimals', ...
            '%s has fewer than two animals represented by the paired sessions.', ...
            alignment.name);
    end
    for b = 1:2
        keep = ismember(session_ids{b}, paired_sessions);
        traces{b} = traces{b}(keep, :);
        session_ids{b} = session_ids{b}(keep);
    end

    analysis_mask = full_time >= alignment.analysis_range(1) & ...
        full_time <= alignment.analysis_range(2);
    time = full_time(analysis_mask);
    for b = 1:2
        traces{b} = traces{b}(:, analysis_mask);
    end

    fprintf('%s: %d paired sessions; %s n=%d events, %s n=%d events\n', ...
        alignment.name, numel(paired_sessions), ...
        behavior_names{1}, size(traces{1}, 1), ...
        behavior_names{2}, size(traces{2}, 1));

    model = fit_timepoint_lmes(traces, session_ids, behavior_names, time);
    model.q = bh_fdr(model.p);
    [session_difference, difference_sessions] = paired_session_differences( ...
        traces, session_ids, paired_sessions);

    this_result.time = time;
    this_result.individual_traces = traces;
    this_result.event_session = session_ids;
    this_result.paired_sessions = paired_sessions;
    this_result.paired_animals = paired_animals;
    this_result.session_difference = session_difference;
    this_result.difference_sessions = difference_sessions;
    this_result.model = model;
    this_result.baseline_range = alignment.baseline_range;
    this_result.analysis_range = alignment.analysis_range;
    results.(alignment.name) = this_result;

    fig = plot_alignment_result(this_result, alignment.name, band_name, ...
        behavior_names, behavior_role, normalization, alpha);
    file_stem = sprintf('%s_%s_vs_%s_%s_%s', band_name, ...
        safe_name(behavior_names{1}), safe_name(behavior_names{2}), ...
        safe_name(behavior_role), alignment.name);
    print(fig, '-dpng', '-r300', fullfile(output_folder, [file_stem '.png']));
    print(fig, '-vector', '-dsvg', fullfile(output_folder, [file_stem '.svg']));
end

%% Save numerical results
result_stem = sprintf('%s_%s_vs_%s_%s_lme', band_name, ...
    safe_name(behavior_names{1}), safe_name(behavior_names{2}), ...
    safe_name(behavior_role));
save(fullfile(output_folder, [result_stem '.mat']), 'results', '-v7.3');
fprintf('Saved figures and model results in %s\n', output_folder);

%% Local functions
function [time, expected_n] = alignment_time(psth, alignment_name, bin_size)
if strcmp(alignment_name, 'onset')
    expected_n = round(range(psth.hist_range) / bin_size);
    time = psth.hist_range(1) + (0:expected_n - 1) * bin_size;
else
    n_pre = round(abs(psth.range_time_wrap(1)) / bin_size);
    n_post = round(abs(psth.range_time_wrap(2)) / bin_size);
    n_warp = psth.n_bins_time_wrap;
    expected_n = n_pre + n_warp + n_post;
    pre_time = psth.range_time_wrap(1) + (0:n_pre - 1) * bin_size;
    warp_time = linspace(bin_size, 5, n_warp);
    post_time = 5 + (1:n_post) * bin_size;
    time = [pre_time, warp_time, post_time];
end
end

function [traces, session_ids] = collect_behavior_events( ...
    psth_structure, animal_names, alignment_row, behavior_idx, ...
    behavior_role, electrode_to_use, session_filters, expected_n)
traces = [];
session_ids = {};

for j = 1:numel(psth_structure)
    session_name = animal_names{j, 1};
    if ~isempty(session_filters) && ~any(contains(session_name, session_filters))
        continue
    end
    if animal_names{j, 2} ~= electrode_to_use
        continue
    end

    if strcmpi(behavior_role, 'self')
        event_matrix = psth_structure(j).all_behavior_psth_self{alignment_row, behavior_idx};
    else
        event_matrix = psth_structure(j).all_behavior_psth_other{alignment_row, behavior_idx};
    end
    if isempty(event_matrix)
        continue
    end
    if size(event_matrix, 2) < expected_n
        warning('Analyze_band_two_behaviors:ShortPSTH', ...
            '%s has %d rather than %d time bins and was skipped.', ...
            session_name, size(event_matrix, 2), expected_n);
        continue
    end

    event_matrix = event_matrix(:, 1:expected_n);
    traces = [traces; event_matrix]; %#ok<AGROW>
    session_ids = [session_ids; repmat({session_name}, size(event_matrix, 1), 1)]; %#ok<AGROW>
end
end

function [traces, session_ids] = normalize_events( ...
    traces, session_ids, time, baseline_range, normalization)
if isempty(traces)
    return
end

if strcmpi(normalization, 'zscore')
    baseline_mask = time >= baseline_range(1) & time < baseline_range(2);
    if ~any(baseline_mask)
        error('Analyze_band_two_behaviors:Baseline', ...
            'No samples fall in baseline [%.2f, %.2f).', ...
            baseline_range(1), baseline_range(2));
    end
    baseline = traces(:, baseline_mask);
    mu = mean(baseline, 2, 'omitmissing');
    sigma = std(baseline, 0, 2, 'omitmissing');
    valid = isfinite(mu) & isfinite(sigma) & sigma > 0;
    traces = (traces(valid, :) - mu(valid)) ./ sigma(valid);
    session_ids = session_ids(valid);
else
    valid = any(isfinite(traces), 2);
    traces = traces(valid, :);
    session_ids = session_ids(valid);
end
end

function model = fit_timepoint_lmes(traces, session_ids, behavior_names, time)
n_time = numel(time);
model.mean = nan(2, n_time);
model.mean_ci_low = nan(2, n_time);
model.mean_ci_high = nan(2, n_time);
model.difference = nan(1, n_time);
model.difference_se = nan(1, n_time);
model.difference_ci_low = nan(1, n_time);
model.difference_ci_high = nan(1, n_time);
model.p = nan(1, n_time);
model.n_events = nan(2, n_time);
model.n_sessions = nan(1, n_time);
model.n_animals = nan(1, n_time);

session_column = [session_ids{1}; session_ids{2}];
animal_column = animal_from_sessions(session_column);
behavior_column = [repmat(behavior_names(1), size(traces{1}, 1), 1); ...
    repmat(behavior_names(2), size(traces{2}, 1), 1)];

for t = 1:n_time
    y = [traces{1}(:, t); traces{2}(:, t)];
    valid = isfinite(y);
    if ~any(valid)
        continue
    end

    tbl = table;
    tbl.BandPower = y(valid);
    tbl.Animal = categorical(animal_column(valid));
    tbl.Session = categorical(session_column(valid));
    tbl.Behavior = categorical(behavior_column(valid), behavior_names, behavior_names);
    tbl.Session = removecats(tbl.Session);
    tbl.Animal = removecats(tbl.Animal);
    tbl.Behavior = removecats(tbl.Behavior);

    if numel(categories(tbl.Behavior)) < 2
        continue
    end
    sessions_b1 = unique(tbl.Session(tbl.Behavior == behavior_names{1}));
    sessions_b2 = unique(tbl.Session(tbl.Behavior == behavior_names{2}));
    if numel(intersect(sessions_b1, sessions_b2)) < 2
        continue
    end

    try
        lme = fitlme(tbl, ...
            'BandPower ~ 1 + Behavior + (1|Animal) + (1|Session)', ...
            'FitMethod', 'REML');
    catch fit_error
        warning('Analyze_band_two_behaviors:Fit', ...
            'Time %.4f s was not fit: %s', time(t), fit_error.message);
        continue
    end

    coefficient_names = cellstr(string(lme.Coefficients.Name));
    behavior_row = find(contains(coefficient_names, 'Behavior'), 1);
    if isempty(behavior_row)
        continue
    end
    coefficient_ci = coefCI(lme);
    model.difference(t) = lme.Coefficients.Estimate(behavior_row);
    model.difference_se(t) = lme.Coefficients.SE(behavior_row);
    model.difference_ci_low(t) = coefficient_ci(behavior_row, 1);
    model.difference_ci_high(t) = coefficient_ci(behavior_row, 2);
    model.p(t) = lme.Coefficients.pValue(behavior_row);

    prediction_table = table;
    prediction_table.Animal = repmat(tbl.Animal(1), 2, 1);
    prediction_table.Session = repmat(tbl.Session(1), 2, 1);
    prediction_table.Behavior = categorical( ...
        behavior_names(:), behavior_names, behavior_names);
    [predicted_mean, predicted_ci] = predict(lme, prediction_table, ...
        'Conditional', false);
    model.mean(:, t) = predicted_mean;
    model.mean_ci_low(:, t) = predicted_ci(:, 1);
    model.mean_ci_high(:, t) = predicted_ci(:, 2);
    model.n_events(1, t) = sum(tbl.Behavior == behavior_names{1});
    model.n_events(2, t) = sum(tbl.Behavior == behavior_names{2});
    model.n_sessions(t) = numel(unique(tbl.Session));
    model.n_animals(t) = numel(unique(tbl.Animal));
end
end

function animal_ids = animal_from_sessions(session_ids)
animal_ids = cell(size(session_ids));
for i = 1:numel(session_ids)
    parts = strsplit(strtrim(session_ids{i}));
    animal_ids{i} = parts{1};
end
end

function q = bh_fdr(p)
% Benjamini-Hochberg adjusted p-values without requiring mafdr.
q = nan(size(p));
valid_idx = find(isfinite(p));
if isempty(valid_idx)
    return
end
[sorted_p, order] = sort(p(valid_idx));
m = numel(sorted_p);
adjusted = sorted_p .* m ./ (1:m);
for k = m - 1:-1:1
    adjusted(k) = min(adjusted(k), adjusted(k + 1));
end
adjusted = min(adjusted, 1);
restored = nan(1, m);
restored(order) = adjusted;
q(valid_idx) = restored;
end

function [difference, sessions] = paired_session_differences( ...
    traces, session_ids, paired_sessions)
sessions = paired_sessions(:);
difference = nan(numel(sessions), size(traces{1}, 2));
for s = 1:numel(sessions)
    mean_1 = mean(traces{1}(strcmp(session_ids{1}, sessions{s}), :), ...
        1, 'omitmissing');
    mean_2 = mean(traces{2}(strcmp(session_ids{2}, sessions{s}), :), ...
        1, 'omitmissing');
    difference(s, :) = mean_2 - mean_1;
end
end

function fig = plot_alignment_result(result, alignment_name, band_name, ...
    behavior_names, behavior_role, normalization, alpha)
colors = [0.10 0.35 0.80; 0.85 0.20 0.15];
light_colors = 0.78 + 0.22 * colors;
time = result.time;
model = result.model;

fig = figure('Color', 'w', 'Units', 'normalized', ...
    'OuterPosition', [0 0 1 1], 'NumberTitle', 'off', ...
    'Name', sprintf('%s %s vs %s %s', ...
    band_name, behavior_names{1}, behavior_names{2}, alignment_name));
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(layout);
hold(ax1, 'on')
for b = 1:2
    if ~isempty(result.individual_traces{b})
        individual_lines = plot(ax1, time, result.individual_traces{b}', ...
            'Color', light_colors(b, :), 'LineWidth', 0.25);
        set(individual_lines, 'HandleVisibility', 'off')
    end
    fill_ci(ax1, time, model.mean_ci_low(b, :), ...
        model.mean_ci_high(b, :), colors(b, :));
end
mean_lines = gobjects(1, 2);
for b = 1:2
    mean_lines(b) = plot(ax1, time, model.mean(b, :), ...
        'Color', colors(b, :), 'LineWidth', 2);
end
yline(ax1, 0, 'k:', 'HandleVisibility', 'off')
xline(ax1, 0, 'k-', 'HandleVisibility', 'off')
shade_significant(ax1, time, model.q < alpha, [0.35 0.35 0.35]);
ylabel(ax1, power_label(band_name, normalization))
title(ax1, sprintf(['All events and LME means (95%% CI): %s n=%d, ' ...
    '%s n=%d; %d sessions from %d animals'], ...
    behavior_names{1}, size(result.individual_traces{1}, 1), ...
    behavior_names{2}, size(result.individual_traces{2}, 1), ...
    numel(result.paired_sessions), numel(result.paired_animals)))
legend(ax1, mean_lines, behavior_names, 'Location', 'best', 'Box', 'off')
set(ax1, 'TickDir', 'out', 'Layer', 'top')

ax2 = nexttile(layout);
hold(ax2, 'on')
if ~isempty(result.session_difference)
    session_lines = plot(ax2, time, result.session_difference', ...
        'Color', [0.72 0.72 0.72], 'LineWidth', 0.75);
    set(session_lines, 'HandleVisibility', 'off')
end
fill_ci(ax2, time, model.difference_ci_low, ...
    model.difference_ci_high, [0.15 0.15 0.15]);
effect_line = plot(ax2, time, model.difference, 'k', 'LineWidth', 2);
yline(ax2, 0, 'k:', 'HandleVisibility', 'off')
xline(ax2, 0, 'k-', 'HandleVisibility', 'off')
shade_significant(ax2, time, model.q < alpha, [0.35 0.35 0.35]);
xlabel(ax2, alignment_xlabel(alignment_name))
ylabel(ax2, sprintf('%s - %s', behavior_names{2}, behavior_names{1}))
title(ax2, sprintf(['Gray: every paired session mean difference; black: LME ' ...
    'contrast (95%% CI); shading: BH-FDR q < %.3g'], alpha))
legend(ax2, effect_line, 'Behavior contrast', 'Location', 'best', 'Box', 'off')
set(ax2, 'TickDir', 'out', 'Layer', 'top')
linkaxes([ax1 ax2], 'x')
xlim(ax1, [time(1), time(end)])
title(layout, sprintf('%s, %s behaviors (%s), %s', ...
    upper_first(band_name), behavior_role, normalization, ...
    strrep(alignment_name, '_', ' ')))
end

function fill_ci(ax, time, low, high, color)
valid = isfinite(time) & isfinite(low) & isfinite(high);
segments = contiguous_segments(valid);
for s = 1:size(segments, 1)
    idx = segments(s, 1):segments(s, 2);
    patch(ax, [time(idx), fliplr(time(idx))], ...
        [low(idx), fliplr(high(idx))], color, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end
end

function shade_significant(ax, time, mask, color)
if numel(time) < 2 || ~any(mask)
    return
end
yl = ylim(ax);
segments = contiguous_segments(mask);
dt = median(diff(time));
for s = 1:size(segments, 1)
    x0 = time(segments(s, 1)) - dt / 2;
    x1 = time(segments(s, 2)) + dt / 2;
    h = patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], color, ...
        'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    uistack(h, 'bottom')
end
ylim(ax, yl)
end

function segments = contiguous_segments(mask)
edges = diff([false, logical(mask(:)'), false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
segments = [starts(:), stops(:)];
end

function label = power_label(band_name, normalization)
if strcmpi(normalization, 'zscore')
    label = sprintf('%s power (within-event baseline z score)', ...
        upper_first(band_name));
else
    label = sprintf('%s power (saved log_{10} spectral amplitude)', ...
        upper_first(band_name));
end
end

function value = upper_first(value)
value = char(value);
value(1) = upper(value(1));
end

function label = alignment_xlabel(alignment_name)
if strcmp(alignment_name, 'onset')
    label = 'Time from behavior onset (s)';
else
    label = 'Normalized time (onset = 0, offset = 5 s)';
end
end

function name = safe_name(name)
name = regexprep(char(name), '[^A-Za-z0-9_-]', '_');
end
