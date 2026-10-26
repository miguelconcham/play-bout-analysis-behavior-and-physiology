function figure_handle = plot_region_gaussian_centers( ...
        region_results, frequency_limits, selection_criterion, ...
        threshold_hz, histogram_bin_hz, spectral_method)
%PLOT_REGION_GAUSSIAN_CENTERS  Centers, BIC/AIC, and the selected mixture.
%
% One row per region. The left panel is the information criterion for
% every fitted component count. The right panel is the distribution of
% Gaussian centers and, when a model was eligible, the selected mixture.

if nargin < 4
    threshold_hz = [];
end
if nargin < 5 || isempty(histogram_bin_hz)
    histogram_bin_hz = 1;
end
if nargin < 6 || isempty(spectral_method)
    spectral_method = '';
end

n_regions = numel(region_results);
criterion = upper(strtrim(char(string(selection_criterion))));
figure_name = 'Regional Gaussian-center mixtures';
if ~isempty(spectral_method)
    figure_name = sprintf('%s | %s', figure_name, spectral_method);
end
figure_handle = figure('Name', figure_name, 'Color', 'w', ...
    'Position', [40 40 1100 max(280, 230 * n_regions)]);

for region_n = 1:n_regions
    this_region = region_results(region_n);
    subplot(n_regions, 2, (region_n - 1) * 2 + 1)
    plot_criterion(this_region, criterion, region_n == n_regions)

    subplot(n_regions, 2, (region_n - 1) * 2 + 2)
    plot_distribution(this_region, frequency_limits, threshold_hz, ...
        histogram_bin_hz, region_n == n_regions)
end

if isempty(spectral_method)
    title_text = sprintf( ...
        'Gaussian bump centers by region | %s selects the mixture', ...
        criterion);
else
    title_text = sprintf( ...
        '%s Gaussian bump centers by region | %s selects the mixture', ...
        spectral_method, criterion);
end
sgtitle(title_text)
end

function plot_criterion(this_region, criterion, label_x)
model_table = this_region.fit.model_selection;
if isempty(model_table) || height(model_table) == 0
    axis off
    text(0.5, 0.5, this_region.fit.status, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Interpreter', 'none');
    title(this_region.name)
    return
end

hold on
plot(model_table.NComponents, model_table.BIC, '-', ...
    'Color', [0.15 0.15 0.15], 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
plot(model_table.NComponents, model_table.AIC, '-', ...
    'Color', [0.55 0.55 0.55], 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
eligible = model_table.Eligible;
legend_handles = gobjects(0, 1);
legend_labels = {};
if any(eligible)
    legend_handles(end + 1, 1) = plot(model_table.NComponents(eligible), ...
        model_table.BIC(eligible), 'o', 'Color', [0.15 0.15 0.15], ...
        'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerSize', 6);
    legend_labels{end + 1} = 'BIC'; %#ok<AGROW>
    legend_handles(end + 1, 1) = plot(model_table.NComponents(eligible), ...
        model_table.AIC(eligible), 'o', 'Color', [0.55 0.55 0.55], ...
        'MarkerFaceColor', [0.55 0.55 0.55], 'MarkerSize', 6);
    legend_labels{end + 1} = 'AIC'; %#ok<AGROW>
end
ineligible = ~eligible & (isfinite(model_table.BIC) | ...
    isfinite(model_table.AIC));
if any(ineligible)
    plot(model_table.NComponents(ineligible), ...
        model_table.BIC(ineligible), 'o', 'Color', [0.15 0.15 0.15], ...
        'MarkerSize', 6, 'HandleVisibility', 'off');
    plot(model_table.NComponents(ineligible), ...
        model_table.AIC(ineligible), 'o', 'Color', [0.55 0.55 0.55], ...
        'MarkerSize', 6, 'HandleVisibility', 'off');
end

best_n = this_region.fit.best_n_components;
if isfinite(best_n)
    best_row = model_table.NComponents == best_n;
    if strcmp(criterion, 'AIC')
        selected_value = model_table.AIC(best_row);
    else
        selected_value = model_table.BIC(best_row);
    end
    legend_handles(end + 1, 1) = plot(best_n, selected_value, 'o', ...
        'Color', [0.75 0.15 0.55], 'MarkerSize', 11, 'LineWidth', 1.8);
    legend_labels{end + 1} = sprintf('Selected %s', criterion); %#ok<AGROW>
end

xlabel_if(label_x, 'Number of Gaussian components')
ylabel('Information criterion')
if isfinite(best_n)
    title(sprintf('%s | %s selects %d', this_region.name, ...
        criterion, best_n))
else
    title(sprintf('%s | %s', this_region.name, this_region.fit.status), ...
        'Interpreter', 'none')
end
xticks(model_table.NComponents)
if ~isempty(legend_handles)
    legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 7)
end
grid on
box off
end

function plot_distribution(this_region, frequency_limits, threshold_hz, ...
        histogram_bin_hz, label_x)
centers = this_region.fit.center_frequencies(:);
centers = centers(isfinite(centers));
hold on
if ~isempty(centers)
    histogram(centers, 'Normalization', 'pdf', ...
        'BinWidth', histogram_bin_hz, ...
        'FaceColor', region_mixture_color(this_region.name), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.45, ...
        'HandleVisibility', 'off');
end

fitted = strcmp(this_region.fit.status, 'fitted');
component_handles = gobjects(0);
mixture_handle = gobjects(0);
if fitted
    component_table = this_region.fit.best_component_table;
    component_colors = lines(height(component_table));
    component_handles = gobjects(height(component_table), 1);
    for component_n = 1:height(component_table)
        component_handles(component_n) = plot( ...
            this_region.fit.density_frequency, ...
            this_region.fit.component_density(:, component_n), ':', ...
            'Color', component_colors(component_n, :), 'LineWidth', 1.5);
    end
    mixture_handle = plot(this_region.fit.density_frequency, ...
        this_region.fit.mixture_density, 'k', 'LineWidth', 2.2);
end

threshold_handle = gobjects(0);
if ~isempty(threshold_hz) && isfinite(threshold_hz)
    xline(threshold_hz, '--', 'Color', [0.75 0.2 0.2], ...
        'LineWidth', 1.1, 'HandleVisibility', 'off');
    threshold_handle = plot(nan, nan, '--', ...
        'Color', [0.75 0.2 0.2], 'LineWidth', 1.1);
end

if ~fitted
    message = this_region.fit.status;
    if isempty(centers)
        message = sprintf('No bumps\n%s', message);
    end
    text(mean(frequency_limits), 0, message, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Interpreter', 'none', 'FontSize', 8);
end

xlabel_if(label_x, 'Gaussian center frequency (Hz)')
ylabel('Probability density')
title(sprintf('%s | %d centers, %d channels', this_region.name, ...
    this_region.n_peaks, this_region.n_channels))
xlim(frequency_limits(:).')

if fitted
    component_labels = arrayfun(@(component_n) sprintf( ...
        'C%d: %.1f Hz', component_n, ...
        component_table.MeanHz(component_n)), ...
        component_table.Component, 'UniformOutput', false);
    if numel(component_labels) > 4
        legend_handles = mixture_handle;
        legend_labels = {'Mixture'};
    else
        legend_handles = [component_handles; mixture_handle];
        legend_labels = [component_labels; {'Mixture'}];
    end
    if ~isempty(threshold_handle) && all(isgraphics(threshold_handle))
        legend_handles = [legend_handles; threshold_handle];
        legend_labels = [legend_labels; ...
            {sprintf('%g Hz reference', threshold_hz)}];
    end
    legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 7)
end
box off
end

function xlabel_if(label_x, label)
if label_x
    xlabel(label)
end
end

function color = region_mixture_color(name)
switch name
    case 'SupCol'
        color = [0.90 0.45 0.05];
    case 'DLPAG'
        color = [0.16 0.42 0.78];
    case 'LPAG'
        color = [0.13 0.55 0.36];
    case 'VLPAG'
        color = [0.55 0.72 0.16];
    case 'DR'
        color = [0.45 0.28 0.72];
    case 'DMPAG'
        color = [0.05 0.55 0.58];
    case 'PAG'
        color = [0.20 0.38 0.30];
    otherwise
        color = [0.45 0.45 0.45];
end
end
