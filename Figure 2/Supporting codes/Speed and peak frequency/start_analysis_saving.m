function [figure_folder, stats_file, figures_before] = start_analysis_saving(figure_folder, stats_name)
%START_ANALYSIS_SAVING  Create the output folder and start the stats diary.
%
% Command-window text written after this call is copied into the text file.
% figures_before is passed to save_analysis_figures at the end so only
% figures from this run are written.

if exist(figure_folder, 'dir') ~= 7
    mkdir(figure_folder);
end
stats_file = fullfile(figure_folder, stats_name);
if exist(stats_file, 'file') == 2
    delete(stats_file);
end
if strcmp(get(0, 'Diary'), 'on')
    diary off
end
figures_before = findall(0, 'Type', 'figure');
diary(stats_file);
diary on
fprintf('Saving figures and stats to %s\n', figure_folder);
end
