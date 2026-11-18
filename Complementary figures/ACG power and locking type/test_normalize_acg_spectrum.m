function test_normalize_acg_spectrum()
%TEST_NORMALIZE_ACG_SPECTRUM  Row totals and the shared high-low reference.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
    addpath(fileparts(this_dir));
end

stack = [1 3 nan; 0 0 0; nan nan nan];
normalized = normalize_spectrum_rows(stack);
assert(abs(sum(normalized(1, isfinite(normalized(1, :)))) - 1) < 1e-12, ...
    'A spectrum was not divided by its own total power.');
assert(abs(normalized(1, 1) - 0.25) < 1e-12, ...
    'The first bin is not its share of the row total.');
assert(all(~isfinite(normalized(2, :))), 'A zero spectrum was still scaled.');
assert(all(~isfinite(normalized(3, :))), 'A missing spectrum was still scaled.');

reference = [2 2 0];
high_spec = [1 1 0];
low_spec = [0 2 0];
[high_out, low_out] = normalize_spectrum_pair(high_spec, low_spec, reference);
assert(abs(high_out(1) - 0.25) < 1e-12 && abs(high_out(2) - 0.25) < 1e-12, ...
    'High power was not divided by the all-runs total.');
assert(abs(low_out(2) - 0.5) < 1e-12, ...
    'Low power was not divided by the same all-runs total.');
assert(abs(high_out(1) / high_spec(1) - low_out(2) / low_spec(2)) < 1e-12, ...
    'High and low did not share one reference.');

[high_fallback, low_fallback] = normalize_spectrum_pair([2 0], [2 0], [nan nan]);
assert(abs(high_fallback(1) - 0.5) < 1e-12 && abs(low_fallback(1) - 0.5) < 1e-12, ...
    'A missing all-runs spectrum did not use the combined high and low total.');

dt = 0.01;
centers = -1:dt:1;
f_grid = 0.25:0.25:12;
wave = cos(2 * pi * 4 * centers);
broad = wave + cos(2 * pi * 20 * centers);
[wave_norm, wave_range] = normalize_spectrum_like_score(wave, centers, [0.2 0.5], f_grid);
[broad_norm, broad_range] = normalize_spectrum_like_score(broad, centers, [0.2 0.5], f_grid);
[scaled_norm, ~] = normalize_spectrum_like_score(3 * wave, centers, [0.2 0.5], f_grid);
assert(wave_range(2) > 12 && broad_range(2) > 12, ...
    'The OscillationScore mean did not include frequencies above the plotted axis.');
[~, at_four] = min(abs(f_grid - 4));
assert(wave_norm(at_four) > broad_norm(at_four), ...
    'Power above 12 Hz did not reduce the relative 4 Hz peak.');
assert(max(abs(scaled_norm - wave_norm)) < 1e-8, ...
    'Scaling the autocorrelogram changed the normalized spectrum.');
missing = normalize_spectrum_like_score(nan(1, numel(centers)), centers, [0.2 0.5], f_grid);
assert(all(~isfinite(missing)), 'A missing autocorrelogram received a normalized spectrum.');

script_text = fileread(fullfile(this_dir, ...
    'Plot_acg_population_nonoscillatory_description.m'));
assert(~isempty(strfind(script_text, 'normalize_spectrum_like_score(')), ...
    'The left column does not divide by the full-spectrum mean.');
assert(isempty(strfind(script_text, 'normalize_spectrum_rows(all_power)')), ...
    'The left column still divides by the plotted-bin total.');
assert(isempty(strfind(script_text, 'normalize_spectrum_pair(high_power, low_power, all_power)')), ...
    'High and low still share the all-runs total.');
assert(~isempty(strfind(script_text, 'mark_normalization_range')), ...
    'The mean panel does not mark the frequencies used for normalization.');

disp('test_normalize_acg_spectrum passed')
end
