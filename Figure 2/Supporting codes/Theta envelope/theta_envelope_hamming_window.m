function window = theta_envelope_hamming_window(n_samples)
%THETA_ENVELOPE_HAMMING_WINDOW  Symmetric Hamming window, as a column.

n_samples = round(n_samples);
if n_samples < 2
    window = 1;
    return
end
k = (0:n_samples - 1).';
window = 0.54 - 0.46 * cos(2 * pi * k / (n_samples - 1));
end
