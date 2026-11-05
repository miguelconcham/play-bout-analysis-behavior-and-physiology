function signal = irasa_simulate_fractal(n_samples, sampling_rate, exponent)
%IRASA_SIMULATE_FRACTAL  Power-law series used to test IRASA.
%
% Section 2.5 of Wen & Liu (2016). The amplitude spectrum follows
% f^(-exponent/2), so the power follows f^(-exponent). Phases are drawn
% independently from a uniform distribution on [0, 2*pi], except DC and
% Nyquist, which stay real. The inverse FFT is then scaled to unit
% standard deviation. Set the random-number generator before calling
% this if a repeatable series is needed.

if nargin < 3
    error('irasa_simulate_fractal:input', ...
        'n_samples, sampling_rate, and exponent are required.');
end
if ~isscalar(n_samples) || n_samples < 32 || n_samples ~= round(n_samples)
    error('irasa_simulate_fractal:length', ...
        'n_samples must be an integer of at least 32.');
end
if ~isscalar(sampling_rate) || ~isfinite(sampling_rate) || sampling_rate <= 0
    error('irasa_simulate_fractal:rate', ...
        'sampling_rate must be a positive finite scalar.');
end
if ~isscalar(exponent) || ~isfinite(exponent)
    error('irasa_simulate_fractal:exponent', ...
        'exponent must be a finite scalar.');
end

nfft = 2 ^ nextpow2(n_samples);
if nfft < n_samples
    nfft = nfft * 2;
end
n_positive = nfft / 2 + 1;
frequency = (0:n_positive - 1).' * (sampling_rate / nfft);
amplitude = zeros(n_positive, 1);
positive = frequency > 0;
amplitude(positive) = frequency(positive) .^ (-exponent / 2);
phase = 2 * pi * rand(n_positive, 1);
phase(1) = 0;
if mod(nfft, 2) == 0
    phase(end) = 0;
end
positive_spectrum = amplitude .* exp(1i * phase);
spectrum = zeros(nfft, 1);
spectrum(1:n_positive) = positive_spectrum;
spectrum(n_positive + 1:end) = conj(positive_spectrum(end - 1:-1:2));
signal = real(ifft(spectrum));
signal = signal(1:n_samples);
signal = signal / std(signal);
end
