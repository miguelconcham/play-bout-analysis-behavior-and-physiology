%% Analyze_acg_delta_component_phase
% Are the spikes that make the delta side peak of the autocorrelogram also
% locked to the delta trough, in trough neurons, or to the delta peak, in
% peak neurons?
%
% Estimate_acg_delta_component_phase writes one file per animal. Each pair
% inside a firing run whose lag is in delta_component_lag is counted once.
% The symmetric autocorrelogram lists that pair twice, and an interior
% spike of a rhythmic train is listed twice, so the saved phase sample is
% the unique spikes. Complement spikes are the other spikes in the same
% runs. Spikes in the gaps between runs are in neither sample.
%
% Lock labels come from EntireSession in delta_all_neurons_v2.mat, with the
% same split as the locking-period autocorrelograms. The target phase is
% the trough (+/- pi) for trough neurons and the peak (0) for peak neurons.
% Each neuron has one vote. The signed-rank tests ask whether the
% component spikes project onto that target, and whether they do so more
% than the complement spikes of the same neuron.
%
% delta_component_lag must match the range used in Estimate. figure_folder
% empty keeps figures on screen. A folder also receives the printed tests
% and SVG/JPG copies of the figures.

%% Paths and parameters
data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
run('\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

result_folder = fullfile(data_root, 'Analysis results', 'acg delta component phase');
phase_folder  = fullfile(data_root, 'Analysis results', 'phase locking data');

delta_component_lag = [0.2 0.5]; % s, must match Estimate
freq_range          = [1 5];
alpha               = 0.01;
non_entrained_p     = 0.1;
min_spikes          = 8;   % component spikes required for a neuron to vote
areas_to_keep       = {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'};
figure_folder       = '';

source_order = {'trough', 'peak', 'nonentrained'};
source_title = {'Trough-locked', 'Peak-locked', 'Non-entrained'};
target_angle = [pi, 0, NaN];

log_fid = open_statistics_log(figure_folder);

say(log_fid, 'Delta-component spikes and delta phase\n');
say(log_fid, 'Lag window %.0f-%.0f ms. Phase band %g-%g Hz.\n', ...
    1000 * delta_component_lag(1), 1000 * delta_component_lag(2), freq_range(1), freq_range(2));
say(log_fid, ['Each autocorrelogram pair is stored once. ', ...
    'The phase sample is the unique spikes in those pairs.\n']);

%% Labels and saved sessions
load(fullfile(phase_folder, 'delta_all_neurons_v2.mat'), 'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) = {'isRt'};
id_all = double(cluster_id_column(all_neurons));

files = dir(fullfile(result_folder, '*acg_delta_component_phase.mat'));
files = files(~strcmp({files.name}, 'acg_delta_component_phase_animal_names.mat'));
if isempty(files)
    error('Analyze_acg_delta_component_phase:empty', ...
        ['No component-phase files in %s. ', ...
        'Run Estimate_acg_delta_component_phase first.'], result_folder);
end

session = {};
cluster_id = [];
area = {};
group_name = {};
n_component = [];
n_complement = [];
n_pairs = [];
n_acg_entries = [];
comp_proj = [];
rest_proj = [];
comp_frac = [];
comp_angle = [];
comp_mvl = [];
comp_hist = [];
phase_edges = [];

for f = 1:numel(files)
    loaded = load(fullfile(files(f).folder, files(f).name), 'session');
    s = loaded.session;
    if max(abs(s.lag_range(:)' - delta_component_lag)) > 1e-6
        error('Analyze_acg_delta_component_phase:lag', ...
            '%s used lag [%.3f %.3f]. Set the same delta_component_lag as Estimate.', ...
            files(f).name, s.lag_range(1), s.lag_range(2));
    end
    if max(abs(s.osc_range(:)' - freq_range)) > 1e-6
        error('Analyze_acg_delta_component_phase:band', ...
            '%s used %g-%g Hz. Set the same freq_range as Estimate.', ...
            files(f).name, s.osc_range(1), s.osc_range(2));
    end
    if isempty(phase_edges)
        phase_edges = s.phase_edges(:)';
    elseif numel(s.phase_edges) ~= numel(phase_edges) || max(abs(s.phase_edges(:)' - phase_edges)) > 1e-9
        error('Analyze_acg_delta_component_phase:edges', ...
            '%s has a different phase histogram.', files(f).name);
    end
    neurons = s.neurons;
    n = height(neurons);
    if size(s.component_hist, 1) ~= n
        error('Analyze_acg_delta_component_phase:hist', ...
            '%s histogram rows do not match the neuron table.', files(f).name);
    end
    animal_code = char(string(s.animal_code));
    this_group = repmat({''}, n, 1);
    this_area = repmat({''}, n, 1);
    id_col = double(neurons.cluster_id);
    for i = 1:n
        this_area{i} = canonical_kept_area(neurons.area{i}, areas_to_keep);
        if isempty(this_area{i})
            continue
        end
        hit = find(session_matches(all_neurons.session, animal_code) & id_all == id_col(i), 1);
        if isempty(hit)
            this_group{i} = 'unlabeled';
            continue
        end
        this_group{i} = delta_lock_group( ...
            all_neurons.EntireSession.PPCPval(hit), ...
            all_neurons.EntireSession.PreferedAngle(hit), ...
            alpha, non_entrained_p);
        if isempty(this_group{i})
            this_group{i} = 'intermediate';
        end
    end
    keep = ~cellfun(@isempty, this_area) & ~cellfun(@isempty, this_group);
    if ~any(keep)
        continue
    end
    session = [session; repmat({animal_code}, sum(keep), 1)]; %#ok<AGROW>
    cluster_id = [cluster_id; id_col(keep)]; %#ok<AGROW>
    area = [area; this_area(keep)]; %#ok<AGROW>
    group_name = [group_name; this_group(keep)]; %#ok<AGROW>
    n_component = [n_component; neurons.n_component(keep)]; %#ok<AGROW>
    n_complement = [n_complement; neurons.n_complement(keep)]; %#ok<AGROW>
    n_pairs = [n_pairs; neurons.n_pairs(keep)]; %#ok<AGROW>
    n_acg_entries = [n_acg_entries; neurons.n_acg_entries(keep)]; %#ok<AGROW>
    comp_angle = [comp_angle; neurons.component_angle(keep)]; %#ok<AGROW>
    comp_mvl = [comp_mvl; neurons.component_mvl(keep)]; %#ok<AGROW>
    comp_hist = [comp_hist; s.component_hist(keep, :)]; %#ok<AGROW>
    proj_peak = neurons.component_proj_peak(keep);
    proj_trough = neurons.component_proj_trough(keep);
    rest_peak = neurons.complement_proj_peak(keep);
    rest_trough = neurons.complement_proj_trough(keep);
    frac_peak = neurons.component_frac_peak(keep);
    % Stored below as the target-specific columns. Peak and trough are
    % filled after the loop from these raw columns, so keep both.
    comp_proj = [comp_proj; [proj_trough, proj_peak]]; %#ok<AGROW>
    rest_proj = [rest_proj; [rest_trough, rest_peak]]; %#ok<AGROW>
    comp_frac = [comp_frac; [1 - frac_peak, frac_peak]]; %#ok<AGROW>
end

if isempty(group_name)
    error('Analyze_acg_delta_component_phase:none', ...
        'No neurons in %s were matched to a lock label.', strjoin(areas_to_keep, ', '));
end

say(log_fid, 'Areas: %s\n', strjoin(areas_to_keep, ', '));
for g = 1:numel(source_order)
    say(log_fid, '%s: %d\n', source_title{g}, sum(strcmp(group_name, source_order{g})));
end
say(log_fid, 'Between the lock thresholds: %d. Unlabeled: %d.\n', ...
    sum(strcmp(group_name, 'intermediate')), sum(strcmp(group_name, 'unlabeled')));

has_component = n_component > 0;
ratio = n_acg_entries(has_component) ./ n_component(has_component);
if isempty(ratio)
    say(log_fid, 'No neuron had a pair inside the lag window.\n\n');
else
    say(log_fid, ['Symmetric autocorrelogram entries per component spike: ', ...
        'median %.2f (n = %d neurons with at least one pair).\n'], ...
        median(ratio), numel(ratio));
    say(log_fid, ['A value near 2 means interior spikes were listed twice ', ...
        'and were kept once in the phase sample.\n\n']);
end

centers = (phase_edges(1:end - 1) + phase_edges(2:end)) / 2;
curve_colors = [0.20 0.40 0.75; 0.85 0.25 0.20; 0.45 0.45 0.45];

for g = 1:numel(source_order)
    enough = strcmp(group_name, source_order{g}) & n_component >= min_spikes & isfinite(comp_angle);
    if g < 3
        paired = enough & n_complement >= min_spikes;
        report_locked_group(log_fid, source_title{g}, target_angle(g), ...
            comp_proj(enough, g), comp_frac(enough, g), comp_angle(enough), comp_mvl(enough), ...
            comp_proj(paired, g), rest_proj(paired, g), min_spikes);
    else
        report_nonentrained(log_fid, comp_proj(enough, 2), comp_frac(enough, 2), min_spikes);
    end
end

plot_component_histograms(comp_hist, centers, group_name, source_order, source_title, ...
    curve_colors, target_angle, min_spikes);
plot_component_against_complement(comp_proj, rest_proj, n_component, n_complement, ...
    group_name, source_order, source_title, curve_colors, min_spikes);
plot_preferred_angles(comp_angle, n_component, group_name, source_order, source_title, ...
    curve_colors, target_angle, min_spikes);

if ~isempty(figure_folder)
    save_open_figures(figure_folder);
end
if log_fid > 0
    fclose(log_fid);
end


function report_locked_group(fid, title_text, target, projection, fraction, angles, spike_mvl, paired_component, paired_complement, min_spikes)
say(fid, '\n%s, target %.0f deg\n', title_text, target * 180 / pi);
say(fid, '  Neurons with at least %d component spikes: %d\n', min_spikes, numel(projection));
[mu, resultant] = circular_mean_angle(angles);
if isfinite(mu)
    say(fid, '  Circular mean of the component preferred phase: %.0f deg, resultant %.2f\n', ...
        mu * 180 / pi, resultant);
end
finite_mvl = spike_mvl(isfinite(spike_mvl));
if ~isempty(finite_mvl)
    say(fid, '  Median within-neuron MVL of component spikes: %.3f\n', median(finite_mvl));
end
report_signrank(fid, 'Component projection on the target', projection, 0, 'right');
report_signrank(fid, 'Fraction of component spikes on the target half', fraction, 0.5, 'right');
say(fid, '  Paired neurons (complement also has >= %d spikes): %d\n', ...
    min_spikes, numel(paired_component));
report_signrank(fid, 'Component minus complement, projection on the target', ...
    paired_component - paired_complement, 0, 'right');
end


function report_nonentrained(fid, proj_peak, frac_peak, min_spikes)
say(fid, '\nNon-entrained, no target phase\n');
say(fid, '  Neurons with at least %d component spikes: %d\n', min_spikes, numel(proj_peak));
report_signrank(fid, 'Component projection on the peak (two sided)', proj_peak, 0, 'both');
report_signrank(fid, 'Fraction of component spikes on the peak half (two sided)', frac_peak, 0.5, 'both');
end


function report_signrank(fid, label, values, null_value, tail_name)
values = values(isfinite(values));
if numel(values) < 5
    say(fid, '  %s: n = %d, too few neurons for a signed-rank test\n', label, numel(values));
    return
end
p = signrank(values, null_value, 'tail', tail_name);
say(fid, '  %s: n = %d, median = %.3f, signed-rank p = %.4g (%s)\n', ...
    label, numel(values), median(values), p, tail_name);
end


function plot_component_histograms(hist_rows, centers, group_name, source_order, source_title, curve_colors, target_angle, min_spikes)
figure('Name', 'Component spike phase', 'Color', 'w');
for g = 1:numel(source_order)
    subplot(1, numel(source_order), g)
    in_group = strcmp(group_name, source_order{g});
    [mean_hist, n] = mean_unit_histogram(hist_rows(in_group, :), min_spikes);
    hold on
    if g < 3
        on_target = target_half(centers, target_angle(g));
        yl = [0, max(0.05, max(mean_hist) * 1.15)];
        bar(centers(on_target), mean_hist(on_target), 1, ...
            'FaceColor', curve_colors(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.95);
        bar(centers(~on_target), mean_hist(~on_target), 1, ...
            'FaceColor', curve_colors(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.35);
        ylim(yl)
    else
        bar(centers, mean_hist, 1, 'FaceColor', curve_colors(g, :), 'EdgeColor', 'none');
    end
    xlim([-pi, pi])
    xticks([-pi, -pi / 2, 0, pi / 2, pi])
    xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
    xlabel('Delta phase')
    ylabel('Mean fraction of spikes')
    title(sprintf('%s, n = %d', source_title{g}, n))
end
end


function plot_component_against_complement(comp_proj, rest_proj, n_component, n_complement, ...
    group_name, source_order, source_title, curve_colors, min_spikes)
figure('Name', 'Component versus complement', 'Color', 'w');
for g = 1:2
    subplot(1, 2, g)
    paired = strcmp(group_name, source_order{g}) ...
        & n_component >= min_spikes & n_complement >= min_spikes;
    x = rest_proj(paired, g);
    y = comp_proj(paired, g);
    hold on
    scatter(x, y, 22, curve_colors(g, :), 'filled', 'MarkerFaceAlpha', 0.65)
    finite = isfinite(x) & isfinite(y);
    if any(finite)
        lim = [min([x(finite); y(finite); -1]), max([x(finite); y(finite); 1])];
        plot(lim, lim, 'k--')
        xlim(lim)
        ylim(lim)
    end
    axis square
    xlabel('Complement, projection on target')
    ylabel('Component, projection on target')
    title(sprintf('%s, n = %d', source_title{g}, sum(finite)))
end
end


function plot_preferred_angles(angles, n_component, group_name, source_order, source_title, ...
    curve_colors, target_angle, min_spikes)
figure('Name', 'Component preferred phase', 'Color', 'w');
edges = -pi:pi / 8:pi;
for g = 1:numel(source_order)
    subplot(1, numel(source_order), g)
    keep = strcmp(group_name, source_order{g}) & n_component >= min_spikes & isfinite(angles);
    sample = angles(keep);
    if isempty(sample)
        title(sprintf('%s, n = 0', source_title{g}))
        continue
    end
    polarhistogram(sample, edges, 'FaceColor', curve_colors(g, :), 'EdgeColor', 'none');
    hold on
    if g < 3
        rmax = max(1, max(histcounts(sample, edges)));
        polarplot([target_angle(g), target_angle(g)], [0, rmax], 'k-', 'LineWidth', 1.5);
    end
    title(sprintf('%s, n = %d', source_title{g}, numel(sample)))
end
end


function [mean_hist, n] = mean_unit_histogram(hist_rows, min_spikes)
n_bins = size(hist_rows, 2);
mean_hist = zeros(1, n_bins);
n = 0;
if isempty(hist_rows)
    return
end
for row = 1:size(hist_rows, 1)
    total = sum(hist_rows(row, :));
    if total < min_spikes
        continue
    end
    mean_hist = mean_hist + hist_rows(row, :) / total;
    n = n + 1;
end
if n > 0
    mean_hist = mean_hist / n;
end
end


function on_target = target_half(centers, target)
if abs(target) < pi / 2
    on_target = centers >= -pi / 2 & centers <= pi / 2;
else
    on_target = centers < -pi / 2 | centers > pi / 2;
end
end


function [mu, mvl] = circular_mean_angle(angles)
angles = angles(isfinite(angles));
if isempty(angles)
    mu = NaN;
    mvl = NaN;
    return
end
z = mean(exp(1i * angles));
mu = atan2(imag(z), real(z));
mvl = abs(z);
end


function say(fid, varargin)
fprintf(varargin{:});
if fid > 0
    fprintf(fid, varargin{:});
end
end


function fid = open_statistics_log(folder)
fid = -1;
if isempty(folder)
    return
end
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
fid = fopen(fullfile(folder, 'statistics.txt'), 'w');
end


function save_open_figures(folder)
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
figs = findall(0, 'Type', 'figure');
for i = 1:numel(figs)
    name = get(figs(i), 'Name');
    if isempty(name)
        name = sprintf('figure_%d', figs(i).Number);
    end
    name = regexprep(name, '[^\w\-]+', '_');
    saveas(figs(i), fullfile(folder, [name, '.svg']));
    saveas(figs(i), fullfile(folder, [name, '.jpg']));
end
end


function id_col = cluster_id_column(ci)
if ismember('cluster_id', ci.Properties.VariableNames)
    id_col = ci.cluster_id;
elseif ismember('id', ci.Properties.VariableNames)
    id_col = ci.id;
else
    error('Analyze_acg_delta_component_phase:id', ...
        'The delta table has neither cluster_id nor id.');
end
end


function tf = session_matches(session_col, animal_code)
if iscell(session_col)
    tf = strcmp(session_col, animal_code);
    empty = cellfun(@isempty, session_col);
    tf(empty) = false;
    nested = cellfun(@iscell, session_col);
    if any(nested)
        tf(nested) = cellfun(@(x) strcmp(x{1}, animal_code), session_col(nested));
    end
else
    tf = strcmp(string(session_col), animal_code);
end
end


function name = canonical_kept_area(label, areas_to_keep)
name = '';
text = strtrim(char(string(label)));
if isempty(text)
    return
end
tokens = regexp(text, '[,;/]+', 'split');
for k = 1:numel(tokens)
    token = upper(regexprep(strtrim(tokens{k}), '\s+', ''));
    if any(strcmp(token, {'SC', 'SUPCOL', 'SUPERIORCOLLICULUS'}))
        token = 'SC';
    end
    hit = strcmpi(areas_to_keep, token);
    if any(hit)
        name = areas_to_keep{find(hit, 1)};
        return
    end
end
end
