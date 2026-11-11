%% Analyze_locking_acg_slope
% Within-neuron slopes of each autocorrelogram measure on LFP delta power.
%
% This is the simplified model. It does not estimate autocorrelograms and
% it does not fit the session, sequential, PPC, or multi-predictor strata
% models. Those stay in Analyze_locking_period_acg. Run that script once
% so the strata file exists, then run this script.
%
% One regression uses every stratum from every neuron that has enough
% strata. Each neuron has its own intercept. Trough-locked, peak-locked,
% and non-entrained neurons each have their own slope:
%
%   y_ij = alpha_i + beta_g * P_ij + error_ij
%
% y is one autocorrelogram measure, P is LFP delta power, i is the neuron,
% and j is the stratum. The three outcomes are fit separately: folded ACG
% delta power, first-peak prominence, and second-cycle prominence. The
% same fit is then repeated with the stratum firing rate beside LFP power,
% with its own slope in each lock group. The rate-adjusted slope is the
% test of whether the positive trough association remains after firing
% rate. The ratio of the two slopes is the fraction of the slope that
% firing rate removes. A second number, the unique R^2, is the variation
% LFP power still explains after firing rate. Those two numbers differ
% when firing rate and LFP power rise together.
%
% Beside the common slope, each qualifying neuron contributes one slope
% and one Spearman correlation, with and without firing rate. The paired
% change is that neuron's value without firing rate minus its value with
% firing rate. Violin plots keep the per-neuron distributions. One bar
% plot gives the median correlation of each lock group: the bar is filled
% when the mean correlation differs from zero and empty otherwise.
% The comparison with firing rate is a square scatter: the original slope
% on x, the rate-adjusted slope on y, the identity line, and the mean as
% a cross. The title is the paired test of that difference.
%
% For trough-locked neurons, one figure shows five neurons near the median
% correlation of each autocorrelogram measure with LFP delta power, and
% another shows up to five neurons whose correlation is significant. The
% axis of each example follows that neuron's values. A normalized
% autocorrelogram has a small delta-band power, so those values can sit
% near 1e-13; that scale is reported and is not stretched to -1 and 1.
% A matrix for each lock group shows whether the same neurons have high
% correlations on more than one measure. For trough-locked neurons, the
% scatter of those correlation values sits beside the scatter of the
% measures themselves. The first principal component of the three
% standardized measures is correlated with LFP delta power in one more figure.

data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';

parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
acg_method = 'eligible';             % 'eligible' or 'overlap'. Loads that summary.
maximize_strata_number = false;      % true: locking_period_acg_summary[_eligible]_max_strata.mat
min_strata = 10;                     % finite strata required before a neuron enters
min_clean = 6;                       % strata required after within-neuron outlier removal
slope_alpha = 0.05;                  % per-neuron threshold for a significant slope or correlation
show_figures = true;
figure_folder = '';                  % empty: figures stay on screen and nothing is written

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
addpath(script_dir);

summary_path = fullfile(parsing_folder, ...
    locking_summary_filename(maximize_strata_number, acg_method));
if exist(summary_path, 'file') ~= 2
    error('Analyze_locking_acg_slope:summary', ...
        ['No strata file at %s. Run Analyze_locking_period_acg once ' ...
        'with acg_method ''%s''.'], summary_path, acg_method);
end
loaded = load(summary_path, 'locking_acg');
if ~isfield(loaded, 'locking_acg') || ~isfield(loaded.locking_acg, 'power_strata')
    error('Analyze_locking_acg_slope:summary', ...
        'The file %s has no power strata.', summary_path);
end
options = struct( ...
    'min_strata', min_strata, ...
    'min_clean', min_clean, ...
    'alpha', slope_alpha, ...
    'do_plot', show_figures, ...
    'figure_folder', figure_folder, ...
    'acg_method', acg_method);
locking_acg_within_neuron(loaded.locking_acg, options);
