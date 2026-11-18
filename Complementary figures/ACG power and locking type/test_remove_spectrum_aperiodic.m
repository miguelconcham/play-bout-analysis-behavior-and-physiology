function test_remove_spectrum_aperiodic()
%TEST_REMOVE_SPECTRUM_APERIODIC  A power law flattens and a peak remains.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

f = 0.25:0.25:12;
offset = 0.4;
exponent = 1.5;
law = 10 .^ (offset - exponent .* log10(f));
flat_law = remove_spectrum_aperiodic(law, f);
assert(max(abs(flat_law)) < 1e-8, ...
    'A pure power-law spectrum was not removed.');

bump = 0.8 * exp(-((f - 4) .^ 2) ./ (2 * 0.35 ^ 2));
peaked = 10 .^ (offset - exponent .* log10(f) + bump);
flat_peak = remove_spectrum_aperiodic(peaked, f);
[~, peak_at] = max(bump);
outside = bump < 0.05;
assert(flat_peak(peak_at) > 0.4, ...
    'The oscillatory bump was removed with the aperiodic fit.');
assert(max(abs(flat_peak(outside))) < 0.08, ...
    'Frequencies away from the bump still carry the power law.');

missing = remove_spectrum_aperiodic([nan(1, numel(f)); law], f);
assert(all(~isfinite(missing(1, :))), 'A missing spectrum received an aperiodic residual.');
assert(max(abs(missing(2, :))) < 1e-8, ...
    'A valid row was changed by a missing row above it.');

empty_flat = remove_spectrum_aperiodic([], f);
assert(isempty(empty_flat), 'An empty spectrum stack was filled in.');

script_text = fileread(fullfile(this_dir, ...
    'Plot_acg_population_nonoscillatory_description.m'));
assert(~isempty(strfind(script_text, 'remove_spectrum_aperiodic(all_power, f_grid)')), ...
    'The all-runs aperiodic column does not use the unnormalized power.');
assert(~isempty(strfind(script_text, 'remove_spectrum_aperiodic(high_power, f_grid)')), ...
    'The high-power aperiodic column does not use the unnormalized power.');
assert(~isempty(strfind(script_text, 'remove_spectrum_aperiodic(low_power, f_grid)')), ...
    'The low-power aperiodic column does not use the unnormalized power.');
assert(~isempty(strfind(script_text, 'subplot(2, 2,')), ...
    'The all-runs spectrum is not drawn in two columns.');
assert(~isempty(strfind(script_text, 'subplot(3, 2,')), ...
    'The high and low spectra are not drawn in two columns.');

disp('test_remove_spectrum_aperiodic passed')
end
