%% Estimate_psth_play_bout_lme
% Mean play-bout power for each recording day and partner, then the mixed model
% loaded by Figure 2.m (section 11) and Play partner.m (section 6).
%
% A session is one recording day (one NPX folder), electrode 1. Each day
% contributes one mean PSTH per partner. At every time bin the model is
%
%   Power ~ 1 + Partner + (1|Session)
%
% Partner is the fixed effect. The recording day is the repeated measure.
% results.est is the Partner1 level. results.partner_est is Partner2 minus
% Partner1. Both series, their SEs, CIs, Cohen d, and p values are saved in
% the same fields Figure 2 already reads.
%
% Samples at or after the bout end are left out. |z| above zscore_limit is
% left out before the day mean. z uses the baseline (-2, 0) and a 20-bin
% moving mean. Partner 3 is left out.
%
% Partner number comes from Behavior already stored in the PSTH, with the
% same rule as the original delta script. Session_structure rows are kept,
% Tickling rows are removed, and a bout takes the first block that contains
% it completely (bout start >= block start and bout end < block end).
%
% Input mats come from Estimate_psth_band_power.m. The names below are the
% delta play-bout pair. Swap them for theta or gamma (see the combinations file).

%% Paths
this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

saving_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';

psth_file = 'psth_structure_delta_updated.mat';
names_file = 'animal_names_delta_updated.mat';
results_file = 'results_play_bout_PBonly_zscore4_updated.mat';

%% Model settings
smooth_wind = 20;
baseline_range = [-2 0];
zscore_limit = 4;
time_range = [-5 5];
length_limit = 0;
electrode_to_use = 1;
animal_label = {'B1D1', 'B1S3', 'B2S2', 'B3D2', 'B4S2', 'B4D4'};

%% Load the play-bout PSTH
disp('loading')
load(fullfile(saving_folder, psth_file), 'psth_structure');
load(fullfile(saving_folder, names_file), 'animal_names');
if numel(psth_structure) ~= size(animal_names, 1)
    error('Estimate_psth_play_bout_lme:Names', ...
        'animal_names has %d rows and psth_structure has %d entries.', ...
        size(animal_names, 1), numel(psth_structure));
end

%% Z-score each bout, then average within recording day and partner
bin_size = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
hist_range = psth_structure(1).hist_range;
n_bins = size(psth_structure(1).play_bout_onset, 2);
time = hist_range(1) + (0:n_bins - 1) * bin_size;
baseline_index = time > baseline_range(1) & time < baseline_range(2);
if ~any(baseline_index)
    error('Estimate_psth_play_bout_lme:Baseline', ...
        'No time bins fall inside (%.2f, %.2f) s.', baseline_range(1), baseline_range(2));
end

power = [];
partner = [];
bout_length = [];
session_id = {};

for j = 1:numel(psth_structure)
    session_name = animal_names{j, 1};
    if ~any(contains(session_name, animal_label))
        continue
    end
    if animal_names{j, 2} ~= electrode_to_use
        continue
    end

    onset = psth_structure(j).play_bout_onset;
    play_bouts = bout_times(psth_structure(j).play_bouts_table);
    if isempty(onset) || isempty(play_bouts)
        continue
    end
    if size(onset, 2) ~= numel(time)
        error('Estimate_psth_play_bout_lme:Time', ...
            '%s has %d bins; the time axis has %d.', session_name, size(onset, 2), numel(time));
    end
    if size(play_bouts, 1) ~= size(onset, 1)
        error('Estimate_psth_play_bout_lme:Bouts', ...
            '%s has %d bouts and %d PSTH rows.', session_name, size(play_bouts, 1), size(onset, 1));
    end

    this_partner = partner_index_from_behavior(play_bouts, psth_structure(j).Behavior);
    fprintf('%s: %d bouts, partner 1 = %d, partner 2 = %d\n', ...
        session_name, numel(this_partner), sum(this_partner == 1), sum(this_partner == 2));

    z = baseline_zscore(onset, baseline_index, smooth_wind);
    power = [power; z]; %#ok<AGROW>
    partner = [partner; this_partner]; %#ok<AGROW>
    bout_length = [bout_length; play_bouts(:, 2) - play_bouts(:, 1)]; %#ok<AGROW>
    session_id = [session_id; repmat({session_name}, size(z, 1), 1)]; %#ok<AGROW>
end

if isempty(power)
    error('Estimate_psth_play_bout_lme:Empty', ...
        'No electrode-%d bouts matched the animal list.', electrode_to_use);
end

obs = session_partner_mean_psth(power, partner, session_id, bout_length, time, ...
    zscore_limit, length_limit);

%% Fit Power ~ 1 + Partner + (1|Session) at each time bin
valid_bins = find(time >= time_range(1) & time <= time_range(2));
if isempty(valid_bins)
    error('Estimate_psth_play_bout_lme:Time', ...
        'No time bins fall inside [%.2f, %.2f] s.', time_range(1), time_range(2));
end
limited_time = time(valid_bins);
power_window = obs.power(:, valid_bins);
n_obs = size(power_window, 1);
n_time = numel(limited_time);

[obs_idx, time_idx] = ndgrid(1:n_obs, 1:n_time);
tbl = table;
tbl.Power = power_window(:);
tbl.Session = categorical(obs.session_id(obs_idx(:)));
tbl.Partner = categorical(obs.partner(obs_idx(:)), [1 2], {'Partner1', 'Partner2'});
tbl.TimeIndex = time_idx(:);
tbl = tbl(isfinite(tbl.Power), :);

[est, se, pvals, d, partner_est, partner_se, partner_pvals, partner_d] = deal(nan(n_time, 1));
ci = nan(n_time, 2);
partner_ci = nan(n_time, 2);

for i = 1:n_time
    tbl_t = tbl(tbl.TimeIndex == i, :);
    if numel(unique(tbl_t.Session)) < 2
        continue
    end
    tbl_t.Partner = removecats(tbl_t.Partner);
    tbl_t.Session = removecats(tbl_t.Session);
    if numel(categories(tbl_t.Partner)) < 2
        continue
    end

    try
        lme = fitlme(tbl_t, 'Power ~ 1 + Partner + (1|Session)');
    catch fit_error
        warning('Estimate_psth_play_bout_lme:Fit', ...
            'Time %.3f s was not fit: %s', limited_time(i), fit_error.message);
        continue
    end

    coef = lme.Coefficients;
    intercept_row = coefficient_row(coef, '(Intercept)');
    partner_row = coefficient_row(coef, 'Partner');
    if isempty(intercept_row) || isempty(partner_row)
        warning('Estimate_psth_play_bout_lme:Coefficients', ...
            'Time %.3f s did not return an intercept and a partner term.', limited_time(i));
        continue
    end

    ci_t = coefCI(lme);
    est(i) = coef.Estimate(intercept_row);
    se(i) = coef.SE(intercept_row);
    pvals(i) = coef.pValue(intercept_row);
    ci(i, :) = ci_t(intercept_row, :);
    partner_est(i) = coef.Estimate(partner_row);
    partner_se(i) = coef.SE(partner_row);
    partner_pvals(i) = coef.pValue(partner_row);
    partner_ci(i, :) = ci_t(partner_row, :);

    sd_within = std(tbl_t.Power, 'omitmissing');
    d(i) = est(i) / sd_within;
    partner_d(i) = partner_est(i) / sd_within;
end

results.time = limited_time(:);
results.est = est;
results.se = se;
results.ci = ci;
results.d = d;
results.pvals = pvals;
results.pvals_fdr = fdr_bh(pvals);
results.partner_est = partner_est;
results.partner_se = partner_se;
results.partner_ci = partner_ci;
results.partner_d = partner_d;
results.partner_pvals = partner_pvals;
results.partner_pvals_fdr = fdr_bh(partner_pvals);
results.zscore_limit = zscore_limit;
results.time_range = time_range;
results.formula = 'Power ~ 1 + Partner + (1|Session)';
results.session = 'recording day';
results.electrode = electrode_to_use;

disp('saving')
save(fullfile(saving_folder, results_file), 'results');
fprintf('Recording days: %d. Day-partner means: %d. Time bins: %d.\n', ...
    numel(unique(obs.session_id)), n_obs, n_time);
fprintf('Saved %s\n', fullfile(saving_folder, results_file));

%% Local functions
function z = baseline_zscore(onset, baseline_index, smooth_wind)
baseline = onset(:, baseline_index);
mu = mean(baseline, 2, 'omitmissing');
sd = std(baseline, 0, 2, 'omitmissing');
sd(~isfinite(sd) | sd == 0) = NaN;
z = (onset - mu) ./ sd;
z = movmean(z, smooth_wind, 2);
end

function play_bouts = bout_times(play_bouts)
if is_table(play_bouts)
    play_bouts = play_bouts{:, 1:min(2, size(play_bouts, 2))};
end
play_bouts = double(play_bouts);
end

function partner_index = partner_index_from_behavior(play_bouts, behavior)
partner_sessions = behavior(ismember(behavior.Animal, 'Session_structure'), :);
partner_sessions(ismember(partner_sessions.Type, 'Tickling'), :) = [];
partner_index = partner_index_in_windows( ...
    play_bouts(:, 1), play_bouts(:, 2), partner_sessions.Start, partner_sessions.End);
end

function tf = is_table(value)
tf = false;
if exist('istable', 'builtin') == 5 || exist('istable', 'file') == 2
    tf = istable(value);
end
end

function row = coefficient_row(coef, term)
names = coef.Name;
if isstring(names)
    names = cellstr(names);
elseif ~iscell(names)
    names = cellstr(string(names));
end
if strcmp(term, '(Intercept)')
    row = find(strcmp(names, '(Intercept)'), 1);
else
    row = find(contains(names, term) & ~strcmp(names, '(Intercept)'), 1);
end
end

function p_fdr = fdr_bh(pvals)
p_fdr = nan(size(pvals));
ok = isfinite(pvals);
if ~any(ok)
    return
end
corrected = mafdr(pvals(ok), 'BHFDR', true);
p_fdr(ok) = corrected(:);
end
