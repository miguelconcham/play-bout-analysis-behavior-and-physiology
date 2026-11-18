function animal_counts = add_population_logit_counts(animal_counts, plot_counts)
%ADD_POPULATION_LOGIT_COUNTS  Attach NPopulationLogit to the good-unit table.
%
% plot_counts empty, and not a table, leaves NPopulationLogit missing
% because the shuffle checkpoint was not read. A table, including one
% with no rows, fills the column with the population-logit count. An
% animal that appears only in Plot_acg_population_logit is added with
% NNeurons 0.

if ~istable(animal_counts)
    animal_counts = table(cell(0, 1), cell(0, 1), zeros(0, 1), zeros(0, 1), ...
        cell(0, 1), 'VariableNames', ...
        {'Animal', 'Recording', 'NNeurons', 'NSessions', 'Sessions'});
end
if ~istable(plot_counts)
    animal_counts.NPopulationLogit = nan(height(animal_counts), 1);
    return
end
n_plot = zeros(height(animal_counts), 1);
used = false(height(plot_counts), 1);
for row = 1:height(animal_counts)
    hit = plot_row(plot_counts, animal_counts.Animal{row}, animal_counts.Recording{row});
    if ~isempty(hit)
        n_plot(row) = plot_counts.NPopulationLogit(hit);
        used(hit) = true;
    end
end
animal_counts.NPopulationLogit = n_plot;
for row = 1:height(plot_counts)
    if used(row)
        continue
    end
    extra = table({plot_counts.Animal{row}}, {plot_counts.Recording{row}}, ...
        0, plot_counts.NSessions(row), {plot_counts.Sessions{row}}, ...
        plot_counts.NPopulationLogit(row), 'VariableNames', ...
        {'Animal', 'Recording', 'NNeurons', 'NSessions', 'Sessions', ...
        'NPopulationLogit'});
    animal_counts = [animal_counts; extra]; %#ok<AGROW>
end
if height(animal_counts) > 1
    [~, order] = sort(strcat(animal_counts.Animal, ' ', animal_counts.Recording));
    animal_counts = animal_counts(order, :);
end
end


function hit = plot_row(plot_counts, animal, recording)
hit = [];
for row = 1:height(plot_counts)
    if strcmp(plot_counts.Animal{row}, animal) ...
            && strcmp(plot_counts.Recording{row}, recording)
        hit = row;
        return
    end
end
end
