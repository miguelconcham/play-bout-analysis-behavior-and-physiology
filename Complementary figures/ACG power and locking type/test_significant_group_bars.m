function test_significant_group_bars()
%TEST_SIGNIFICANT_GROUP_BARS  Unpaired summary keeps only significant features.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

y_osc = [1; 1.2; 0.8; 1.1];
y_non = [4; 4.2; 3.9; 4.1];
values = [y_osc, y_osc + 10, y_osc, y_non; ...
    y_non, y_non + 10, y_non, y_osc];
group_index = [ones(4, 1); 2 * ones(4, 1)];
titles = {'Modulation1', 'PeakFrequency', 'Gamma R2', 'AUC'};
families = {'oscillatory', 'oscillatory', 'non-oscillatory', 'non-oscillatory'};
holm_p = [0.0004; 0.20; 0.004; 0.03];
panels = significant_group_bars(values, group_index, titles, families, holm_p);

assert(numel(panels) == 2, 'The summary did not keep two feature families.');
assert(isequal(panels(1).titles, {'Modulation1'}), ...
    'The oscillatory panel did not keep only its significant feature.');
assert(isequal(panels(2).titles, {'Gamma R2', 'AUC'}), ...
    'The non-oscillatory panel did not keep both significant features.');
assert(strcmp(panels(1).stars{1}, '***'), ...
    'A Holm-Sidak p below 0.001 was not drawn as three stars.');
assert(strcmp(panels(2).stars{1}, '**'), ...
    'A Holm-Sidak p below 0.01 was not drawn as two stars.');
assert(strcmp(panels(2).stars{2}, '*'), ...
    'A Holm-Sidak p below 0.05 was not drawn as one star.');

pooled = [y_osc; y_non];
expected_osc = mean((y_osc - mean(pooled)) / std(pooled, 0));
expected_non = mean((y_non - mean(pooled)) / std(pooled, 0));
assert(abs(panels(1).mean1(1) - expected_osc) < 1e-12, ...
    'The oscillatory bar is not the mean z-score of that group.');
assert(abs(panels(1).mean2(1) - expected_non) < 1e-12, ...
    'The non-oscillatory bar is not the mean z-score of that group.');
weight = numel(panels(1).z1{1}) * panels(1).mean1(1) ...
    + numel(panels(1).z2{1}) * panels(1).mean2(1);
assert(abs(weight) < 1e-10, ...
    'The two group means are not centered on the pooled mean.');

same = significant_group_bars([1; 2; 1; 2], [1; 1; 2; 2], ...
    {'Baseline'}, {'oscillatory'}, 0.04);
assert(abs(same(1).mean1) < 1e-12 && abs(same(1).mean2) < 1e-12, ...
    'Equal groups did not land on a zero z-score.');
assert(isempty(same(2).titles), ...
    'An empty non-oscillatory family still kept a feature.');

boundary = significant_group_bars([1 1 1; 2 2 2; 8 8 8; 9 9 9], [1; 1; 2; 2], ...
    {'Peak1', 'Peak2', 'Trough1'}, {'oscillatory', 'oscillatory', 'oscillatory'}, ...
    [0.05; 0.01; 0.001]);
assert(isequal(boundary(1).titles, {'Peak2', 'Trough1'}), ...
    'A Holm-Sidak p of 0.05 was kept.');
assert(strcmp(boundary(1).stars{1}, '*'), ...
    'A Holm-Sidak p of 0.01 was not one star.');
assert(strcmp(boundary(1).stars{2}, '**'), ...
    'A Holm-Sidak p of 0.001 was not two stars.');

script_text = fileread(fullfile(this_dir, ...
    'Plot_acg_population_nonoscillatory_description.m'));
assert(~isempty(strfind(script_text, 'significant_group_bars(')), ...
    'Section 2 does not build the unpaired summary.');
assert(~isempty(strfind(script_text, 'wide-bars')), ...
    'The summary panels would be forced into a square.');
draw_text = fileread(fullfile(this_dir, 'draw_significant_group_bars.m'));
assert(~isempty(strfind(draw_text, 'subplot(2, 1, g)')), ...
    'The summary does not use one subplot per feature family.');
assert(isempty(strfind(draw_text, 'ksdensity')), ...
    'The summary still draws a violin.');
assert(isempty(strfind(draw_text, 'draw_one_violin')), ...
    'The summary still draws a violin.');

fig = draw_significant_group_bars(panels, 'Ward');
tagged = findall(fig, 'Tag', 'wide-bars');
assert(numel(tagged) == 2, 'The summary did not draw two wide panels.');
star_text = findall(fig, 'Type', 'text', 'String', '***');
assert(~isempty(star_text), 'The summary figure did not place the stars.');
print(fig, '/tmp/group_t_summary_test.jpg', '-djpeg', '-r120')
close(fig)

disp('test_significant_group_bars passed')
end
