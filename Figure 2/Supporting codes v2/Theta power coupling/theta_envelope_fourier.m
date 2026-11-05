function coefficient = theta_envelope_fourier(epochs, time, frequencies)
%THETA_ENVELOPE_FOURIER  Least-squares sinusoid at each center frequency.
%
% The epoch and the sine and cosine are demeaned, so a constant level is
% not scored as a slow oscillation. The coefficient is a+ib for
% a*cos(2*pi*f*t) + b*sin(2*pi*f*t). Frequencies that are ill-conditioned
% on a short window, including 0.1 Hz inside 1-2 s, come back NaN.
% Neighboring 0.2 Hz fits on a 2 s window are not independent estimates.

frequencies = frequencies(:).';
n_trials = size(epochs, 1);
n_freq = numel(frequencies);
coefficient = nan(n_trials, n_freq);
valid = all(isfinite(epochs), 2);
if ~any(valid) || isempty(time)
    return
end
y = epochs(valid, :).';
y = y - mean(y, 1);
t = time(:);
for k = 1:n_freq
    cosine = cos(2 * pi * frequencies(k) * t);
    sine = sin(2 * pi * frequencies(k) * t);
    cosine = cosine - mean(cosine);
    sine = sine - mean(sine);
    basis = [cosine, sine];
    gram = basis' * basis;
    if rcond(gram) < 1e-8
        continue
    end
    beta = gram \ (basis' * y);
    coefficient(valid, k) = beta(1, :).' + 1i * beta(2, :).';
end
end
