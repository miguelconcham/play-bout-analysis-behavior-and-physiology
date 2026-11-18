function test_acg_spectrum_stack()
%TEST_ACG_SPECTRUM_STACK  Folded autocorrelogram spectrum on one frequency grid.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
    addpath(fileparts(this_dir));
end

dt = 0.01;
centers = -1:dt:1;
f0 = 4;
wave = cos(2 * pi * f0 * centers);
f_grid = 0.25:0.25:12;
stack = acg_spectrum_stack(wave, centers, [0.2 0.5], f_grid);
assert(isequal(size(stack), [1, numel(f_grid)]), ...
    'One autocorrelogram did not return one spectrum row.');
finite = isfinite(stack);
assert(any(finite), 'A 4 Hz autocorrelogram produced no spectrum.');
finite_grid = f_grid(finite);
finite_spec = stack(finite);
[~, peak_at] = max(finite_spec);
assert(abs(finite_grid(peak_at) - f0) <= 0.75, ...
    'A 4 Hz autocorrelogram did not peak near 4 Hz.');

silent = zeros(1, numel(centers));
both = acg_spectrum_stack([wave; silent], centers, [0.2 0.5], f_grid);
assert(isequal(size(both), [2, numel(f_grid)]), ...
    'Two autocorrelograms did not return two spectrum rows.');
wave_spec = both(1, isfinite(both(1, :)));
flat_spec = both(2, isfinite(both(2, :)));
assert(~isempty(wave_spec) && ~isempty(flat_spec), ...
    'One of the two spectra was empty.');
assert(max(flat_spec) < max(wave_spec), ...
    'The flat autocorrelogram was not weaker than the 4 Hz rhythm.');

blank = acg_spectrum_stack(nan(1, numel(centers)), centers, [0.2 0.5], f_grid);
assert(all(~isfinite(blank)), 'A missing autocorrelogram still received a spectrum.');

empty_stack = acg_spectrum_stack([], centers, [0.2 0.5], f_grid);
assert(isequal(size(empty_stack), [0, numel(f_grid)]), ...
    'An empty autocorrelogram did not keep the frequency grid.');

positive = centers > 0;
[reflected, reflected_lag] = reflect_positive_acg(wave(positive), centers(positive));
reflected_spec = acg_spectrum_stack(reflected, reflected_lag, [0.2 0.5], f_grid);
assert(max(abs(stack - reflected_spec)) < 1e-9, ...
    'Reflecting the positive side changed the spectrum of an even autocorrelogram.');

script_text = fileread(fullfile(this_dir, ...
    'Plot_acg_population_nonoscillatory_description.m'));
assert(~isempty(strfind(script_text, ...
    'draw_oscillatory_acg_spectra(osc_rows, spectrum_curves,')), ...
    'Section 4 does not draw the oscillatory autocorrelogram spectra.');
assert(~isempty(strfind(script_text, 'spectrum_curves.all = imagesc_logit.sum.all')), ...
    'The percentile-logit spectrum does not use the displayed all-runs curve.');
assert(~isempty(strfind(script_text, 'reflect_positive_acg(')), ...
    'The observed spectrum does not use the displayed positive-lag curve.');
assert(~isempty(strfind(script_text, 'ACG spectrum, oscillatory, high and low power')), ...
    'Section 4 does not draw the high and low power spectra.');
assert(isempty(strfind(script_text, 'aligned_two_sided_acgs')), ...
    'The spectrum still loads a separate autocorrelogram.');

disp('test_acg_spectrum_stack passed')
end
