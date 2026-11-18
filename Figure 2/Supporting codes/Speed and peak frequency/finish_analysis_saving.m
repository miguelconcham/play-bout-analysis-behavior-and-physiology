function finish_analysis_saving(figure_folder, stats_file, figures_before)
%FINISH_ANALYSIS_SAVING  Save this run's figures and close the stats diary.

try
    saved = save_analysis_figures(figure_folder, figures_before);
    fprintf('Saved %d figures in %s\n', numel(saved), figure_folder);
    for k = 1:numel(saved)
        fprintf('  %s\n', saved{k});
    end
    fprintf('Stats file: %s\n', stats_file);
catch ME
    diary off
    rethrow(ME)
end
diary off
end
