function y = theta_envelope_phase_randomize(x)
%THETA_ENVELOPE_PHASE_RANDOMIZE  Real surrogate with the same Fourier magnitudes.
%
% DC and Nyquist stay real. Positive-frequency phases are uniform and the
% negative frequencies are their conjugates, so the power spectrum is kept
% and the amplitude modulation is not.

x = double(x(:));
n = numel(x);
spectrum = fft(x);
if mod(n, 2) == 0
    positive = 2:(n / 2);
    nyquist = n / 2 + 1;
else
    positive = 2:((n + 1) / 2);
    nyquist = [];
end
if ~isempty(positive)
    phi = 2 * pi * rand(numel(positive), 1);
    spectrum(positive) = abs(spectrum(positive)) .* exp(1i * phi);
    spectrum(n - positive + 2) = conj(spectrum(positive));
end
spectrum(1) = real(spectrum(1));
if ~isempty(nyquist)
    spectrum(nyquist) = real(spectrum(nyquist));
end
y = real(ifft(spectrum));
end
