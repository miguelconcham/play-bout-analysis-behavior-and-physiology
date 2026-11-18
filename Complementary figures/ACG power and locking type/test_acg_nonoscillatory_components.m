function test_acg_nonoscillatory_components()
%TEST_ACG_NONOSCILLATORY_COMPONENTS  Bump height stays finite.
%
% DiffExpAmplitude is the height of the positive-lag bump. A decay that
% sits close to the rise used to make that coefficient run away, and the
% power-gap regression then returned zeros and NaN standard errors.
% Run this file by name.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

t = 0.005:0.005:0.6;
separated = difference_curve(t, 0.02, 0.12, 0.8, 0.15);
close_times = difference_curve(t, 0.03, 0.0303, 0.5, 0.1);
alpha_curve = 0.1 + 0.6 * alpha_kernel(t, 0.04);
gamma_curve = gamma_curve_of(t, 0.7, 2, 0.05, 0.1);
curves = [separated; close_times; alpha_curve; gamma_curve];
expected_height = [0.8; 0.5; 0.6; 0.7];

features = acg_nonoscillatory_components(curves, t);
heights = features.DiffExpAmplitude;
assert(all(isfinite(heights(1:3))), ...
    'A difference-of-exponentials curve did not return a finite amplitude.');
assert(all(heights(1:3) < 20), ...
    'DiffExpAmplitude left the scale of the bump.');
assert(all(abs(heights(1:3) - expected_height(1:3)) < 0.25), ...
    'DiffExpAmplitude was not the bump height.');
assert(all(features.DiffExpDecay(1:3) > features.DiffExpRise(1:3)), ...
    'Decay was not longer than rise.');
assert(all(features.DiffExpR2(1:3) > 0.9), ...
    'The difference-of-exponentials fit missed a noiseless curve.');
assert(isfinite(features.GammaAmplitude(4)) && features.GammaAmplitude(4) < 20, ...
    'GammaAmplitude left the scale of the bump.');
assert(abs(features.GammaAmplitude(4) - expected_height(4)) < 0.25, ...
    'GammaAmplitude was not the bump height.');
assert(all(isfinite(features.PeakAmplitude)) && all(isfinite(features.AUC)), ...
    'Peak amplitude or area was missing.');

short_lags = [0.01 0.02];
short_curve = [0.2 0.4];
short_features = acg_nonoscillatory_components(short_curve, short_lags);
assert(isnan(short_features.DiffExpAmplitude) && isfinite(short_features.PeakAmplitude), ...
    'A curve with too few lags still received a difference-of-exponentials fit.');
end


function y = difference_curve(t, rise, decay, amplitude, baseline)
rel = (decay - rise) / rise;
raw = exp(-t ./ decay) - exp(-t ./ rise);
t_peak = rise * decay * log(1 + rel) / (decay - rise);
height = exp(-t_peak / rise) * rel;
kernel = raw ./ height;
kernel = kernel ./ max(kernel);
y = baseline + amplitude * kernel;
end


function y = gamma_curve_of(t, amplitude, shape, scale, baseline)
log_g = (shape - 1) * log(t) - t ./ scale;
g = exp(log_g);
g = g ./ max(g);
y = baseline + amplitude * g;
end


function k = alpha_kernel(t, rise)
u = t ./ rise;
k = u .* exp(1 - u);
k = k ./ max(k);
end
