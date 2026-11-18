%% Analyze_distance_freq_relation
% Is delta power higher during play than at the same inter-animal distance
% outside play?
%
% Loads the delta-band speed-power mats and keeps electrode 1. Distance is
% the saved within-session z-score of the separation between the implanted
% animal and its partner. Positive z is farther than that session's mean.
% Power is the saved log10 spectrogram amplitude in the delta file.
%
% Each distance bin keeps sessions that have both play and non-play there.
% The first table and the first figure compare those two powers. A positive
% difference means higher delta power during play. The later output counts
% how many sessions supply each distance, including the far tail (z >= 1)
% and the close tail (z <= -1).
%
% test_distance_session_coverage checks the comparison without these mats.

script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end

this_file = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes\Speed and peak frequency\Analyze_distance_freq_relation.m';
repo_root = fileparts(fileparts(fileparts(fileparts(this_file))));
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');

%% Saving folder for figures and the stats text file
figure_folder = fullfile(data_root, 'Analysis results', 'distance freq relation');
[figure_folder, stats_file, figures_before] = start_analysis_saving( ...
    figure_folder, 'distance_freq_relation_stats.txt');

load([saving_folder, '\psth_structure_speed_delta.mat'], 'psth_structure')
load([saving_folder, '\animal_names_speed_delta.mat'], 'animal_names')

electrode_id = cellfun(@double, animal_names(:, 3));
session_rows = find(electrode_id == 1);
if ~isempty(session_rows) && isfield(psth_structure(session_rows(1)), 'freq_range') ...
        && ~isempty(psth_structure(session_rows(1)).freq_range)
    band = psth_structure(session_rows(1)).freq_range;
    fprintf('Band in this file: %.1f-%.1f Hz\n', band(1), band(2));
end
fprintf(['Distance is the saved within-session z-score of inter-animal separation.\n' ...
    'Positive z is farther than that session mean. Power is log10 spectrogram amplitude.\n']);
print_distance_column_check(psth_structure, session_rows);

options = distance_coverage_options();
sessions = distance_coverage_sessions(psth_structure, animal_names);
report = distance_coverage_report(sessions, options);
print_distance_coverage(report);
plot_distance_coverage(report);

%% Save figures and close the stats text file
finish_analysis_saving(figure_folder, stats_file, figures_before);

function print_distance_column_check(psth_structure, session_rows)
r_relative = nan(numel(session_rows), 1);
r_speed = nan(numel(session_rows), 1);
for s = 1:numel(session_rows)
    mdl = psth_structure(session_rows(s)).model_data;
    if ~has_column(mdl, 'Distance') || ~has_column(mdl, 'RelativeSpeed') ...
            || ~has_column(mdl, 'Speed')
        continue
    end
    distance = column_values(mdl, 'Distance');
    relative = column_values(mdl, 'RelativeSpeed');
    speed = column_values(mdl, 'Speed');
    ok = isfinite(distance(:)) & isfinite(relative(:)) & isfinite(speed(:));
    if nnz(ok) < 10
        continue
    end
    r_relative(s) = corr(distance(ok), relative(ok));
    r_speed(s) = corr(distance(ok), speed(ok));
end

if ~any(isfinite(r_relative)) && ~any(isfinite(r_speed))
    fprintf('Speed or RelativeSpeed is absent, so the column check was skipped.\n');
    return
end
fprintf('Median |corr(Distance, Speed)| = %.3f\n', median(abs(r_speed), 'omitnan'));
fprintf('Median |corr(Distance, RelativeSpeed)| = %.3f\n', median(abs(r_relative), 'omitnan'));
if median(abs(r_relative), 'omitnan') > 0.95
    fprintf(['Distance and RelativeSpeed agree almost perfectly in this file.\n' ...
        'Distance should be the separation and RelativeSpeed its change.\n' ...
        'The coverage below still reads Distance.\n']);
end
end

function present = has_column(mdl, name)
if isstruct(mdl)
    present = isfield(mdl, name);
else
    present = any(strcmp(name, mdl.Properties.VariableNames));
end
end

function values = column_values(mdl, name)
values = double(mdl.(name)(:));
end
