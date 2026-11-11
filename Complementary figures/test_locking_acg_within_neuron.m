function test_locking_acg_within_neuron()
%TEST_LOCKING_ACG_WITHIN_NEURON  Synthetic strata for the simplified slope.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

options = struct('min_strata', 10, 'min_clean', 6, 'alpha', 0.05, ...
    'do_plot', false, 'figure_folder', '', 'acg_method', '');
rng(1);
strata = planted_strata();
result = locking_acg_within_neuron(strata, options);

trough = group_index(result, 'trough');
peak = group_index(result, 'peak');
direct = outcome_index(result, 'acg_delta');
through_rate = outcome_index(result, 'first_prominence');
noise = outcome_index(result, 'second_prominence');

assert(result.model(direct).simple.n_neurons == 60, ...
    'The common slope did not keep the 60 neurons with enough strata.');
assert(result.model(direct).simple.n_by_group(trough) == 20, ...
    'A trough neuron with too few strata entered the common slope.');
assert(result.per_neuron(direct).n_other_group == 1, ...
    'A neuron outside the three lock groups was not left out.');
assert(numel(result.per_neuron(direct).rho{panel_index(result, 'delta')}) == 40, ...
    'Delta-locked correlations were not trough and peak together.');

direct_simple = result.model(direct).simple.beta(trough);
direct_rate = result.model(direct).with_rate.beta(trough);
assert(direct_simple > 0.95 && direct_simple < 1.25, ...
    'The trough slope of the direct autocorrelogram effect was not about 1.1.');
assert(direct_rate > 0.85 && direct_rate < 1.15, ...
    'The trough slope did not remain near 1 after firing rate.');
assert(result.model(direct).share(trough).fraction_of_slope_remaining > 0.8, ...
    'Most of the direct trough slope did not remain after firing rate.');
assert(result.model(direct).simple.p(trough) < 0.01, ...
    'The direct trough slope was not above zero.');
assert(result.model(direct).with_rate.p(trough) < 0.01, ...
    'The direct trough slope was not above zero after firing rate.');
assert(abs(result.model(direct).with_rate.beta(n_type_rate(result, trough)) - 0.2) < 0.08, ...
    'The firing-rate slope in the direct trough model was not about 0.2.');

via_simple = result.model(through_rate).simple.beta(trough);
via_rate = result.model(through_rate).with_rate.beta(trough);
assert(via_simple > 0.35 && via_simple < 0.7, ...
    'An autocorrelogram measure carried only by firing rate had no simple trough slope.');
assert(abs(via_rate) < 0.1, ...
    'The trough slope carried by firing rate remained after firing rate was included.');
assert(abs(result.model(through_rate).share(trough).fraction_of_slope_remaining) < 0.2, ...
    'The fraction of the firing-rate slope that remained was not near zero.');
assert(result.model(through_rate).share(trough).fraction_due_to_rate > 0.9, ...
    'Firing rate did not account for the within-neuron R^2 it generated.');
assert(result.model(direct).share(trough).fraction_remaining_r2 < 0.15, ...
    'The unique R^2 was treated as if it were the remaining slope.');

assert(abs(result.model(direct).simple.beta(peak)) < 0.08, ...
    'Peak-locked neurons showed an LFP slope that was not planted.');
assert(abs(result.model(noise).simple.beta(trough)) < 0.08, ...
    'A noise outcome had a trough slope.');
assert(result.model(direct).simple.contrast_p(1) < 0.01, ...
    'The trough slope was not larger than the peak slope.');
assert(finite_values(result.per_neuron(direct).slope{panel_index(result, 'trough')}) > 0.9, ...
    'The median trough neuron slope was not positive.');
assert(finite_values(result.per_neuron(direct).partial_slope{panel_index(result, 'trough')}) > 0.85, ...
    'The median rate-adjusted trough slope did not remain positive.');
assert(abs(finite_values(result.per_neuron(through_rate).partial_slope{panel_index(result, 'trough')})) < 0.1, ...
    'The median neuron slope carried by firing rate remained after the adjustment.');
assert(result.model(direct).simple.n_neurons == result.model(direct).with_rate.n_neurons, ...
    'The rate-adjusted fit did not use the same neurons as the simple fit.');
trough_cell = find(strcmp(result.cell_order, 'trough'), 1);
peak_cell = find(strcmp(result.cell_order, 'peak'), 1);
direct_bar = result.correlation_bars(direct);
assert(direct_bar.filled(trough_cell) && direct_bar.median_r(trough_cell) > 0.8, ...
    'The trough median correlation was not a filled bar above zero.');
assert(~direct_bar.filled(peak_cell), ...
    'The peak correlation bar was filled when its mean did not differ from zero.');
rate_change = result.slope_change(through_rate);
assert(rate_change.filled(trough_cell) && rate_change.median(trough_cell) > 0.3, ...
    'Adding firing rate did not change the trough slope that was carried by rate.');
assert(~result.slope_change(noise).filled(trough_cell), ...
    'A noise slope was marked as changed by firing rate.');
assert(result.joint_slope_change(trough_cell).changes, ...
    'Trough-locked neurons were not marked when a slope changed after firing rate.');

constant = constant_rate_strata();
constant_result = locking_acg_within_neuron(constant, options);
assert(abs(constant_result.model(direct).simple.beta(trough) - 2) < 1e-8, ...
    'A constant firing rate changed the simple slope.');
assert(abs(constant_result.model(direct).with_rate.beta(trough) - 2) < 1e-8, ...
    'A constant firing rate changed the rate-adjusted LFP slope.');

summary = struct('parameters', struct('acg_method', 'overlap'), 'power_strata', strata);
bad_options = options;
bad_options.acg_method = 'eligible';
failed = false;
try
    locking_acg_within_neuron(summary, bad_options);
catch err
    failed = ~isempty(strfind(err.message, 'acg_method')); %#ok<STREMP>
end
assert(failed, 'A summary from the other autocorrelogram method was accepted.');

if exist('graphics_toolkit', 'file') == 2 || exist('graphics_toolkit', 'builtin') == 5
    try
        graphics_toolkit('gnuplot');
    catch
    end
end
folder = fullfile(tempdir, 'acg_slope_test');
if exist(folder, 'dir') == 7
    rmdir(folder, 's');
end
plot_options = options;
plot_options.do_plot = true;
plot_options.figure_visible = 'off';
plot_options.figure_folder = folder;
locking_acg_within_neuron(strata, plot_options);
assert(exist(fullfile(folder, 'within_neuron_acg_slopes.txt'), 'file') == 2, ...
    'The slope report was not written.');
assert(exist(fullfile(folder, 'Neuron slopes, violins.png'), 'file') == 2, ...
    'The slope violin figure was not written.');
assert(exist(fullfile(folder, 'Neuron slopes, rank violins.png'), 'file') == 2, ...
    'The rank violin figure was not written.');
assert(exist(fullfile(folder, 'Median Spearman by lock group.png'), 'file') == 2, ...
    'The median correlation bars were not written.');
assert(exist(fullfile(folder, 'Slope change after firing rate.png'), 'file') == 2, ...
    'The slope-change figure was not written.');
assert(exist(fullfile(folder, 'Common within-neuron slopes.png'), 'file') == 2, ...
    'The common-slope figure was not written.');
close all
disp('locking_acg_within_neuron tests passed.')
end


function value = finite_values(values)
values = values(isfinite(values));
value = median(values);
end


function k = group_index(result, name)
k = find(strcmp(result.type_order, name), 1);
end


function k = outcome_index(result, name)
k = find(strcmp(result.outcome_fields, name), 1);
end


function k = panel_index(result, name)
k = find(strcmp(result.panel_names, name), 1);
end


function k = n_type_rate(result, group)
k = numel(result.type_order) + group;
end


function strata = planted_strata()
P = (1:12)' - 6.5;
groups = {'trough', 'peak', 'nonentrained'};
n_per = 20;
n = n_per * numel(groups) * 12 + 3 + 12;
strata = blank_strata(n);
row = 0;
for g = 1:numel(groups)
    for i = 1:n_per
        rate = 0.5 * P + 0.35 * randn(12, 1);
        direct = 1.0 * P + 0.2 * rate + 0.05 * randn(12, 1);
        via_rate = 1.0 * rate + 0.05 * randn(12, 1);
        if strcmp(groups{g}, 'trough')
            acg = direct;
            first = via_rate;
        else
            acg = 0.05 * randn(12, 1);
            first = 0.05 * randn(12, 1);
        end
        second = 0.05 * randn(12, 1);
        [strata, row] = write_neuron(strata, row, sprintf('%s-%02d', groups{g}, i), ...
            groups{g}, P, rate, acg, first, second);
    end
end
[strata, row] = write_neuron(strata, row, 'short', 'trough', P(1:3), ...
    0.5 * P(1:3), P(1:3), P(1:3), P(1:3));
strata = write_neuron(strata, row, 'other', 'delta', P, 0.5 * P, P, P, P);
end


function strata = constant_rate_strata()
P = (1:12)' - 6.5;
strata = blank_strata(12);
strata = write_neuron(strata, 0, 'constant', 'trough', P, 3 * ones(12, 1), ...
    2 * P, zeros(12, 1), zeros(12, 1));
end


function strata = blank_strata(n)
strata.n = n;
strata.neuron = cell(n, 1);
strata.lock_group = cell(n, 1);
strata.power = nan(n, 1);
strata.rate = nan(n, 1);
strata.acg_delta = nan(n, 1);
strata.first_prominence = nan(n, 1);
strata.second_prominence = nan(n, 1);
end


function [strata, row] = write_neuron(strata, row, id, group, power, rate, acg, first, second)
n = numel(power);
for k = 1:n
    strata.neuron{row + k} = id;
    strata.lock_group{row + k} = group;
    strata.power(row + k) = power(k);
    strata.rate(row + k) = rate(k);
    strata.acg_delta(row + k) = acg(k);
    strata.first_prominence(row + k) = first(k);
    strata.second_prominence(row + k) = second(k);
end
row = row + n;
end
