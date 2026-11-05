function ppc = theta_envelope_phase_consistency(coefficients)
%THETA_ENVELOPE_PHASE_CONSISTENCY  Pairwise phase consistency across behaviors.
%
% PPC = (n*R^2 - 1)/(n - 1), with R the mean resultant length of unit
% phasors. It does not grow just because a recording has more onsets.
% Fewer than two finite trials leaves NaN.

n_freq = size(coefficients, 2);
ppc = nan(1, n_freq);
for f = 1:n_freq
    trial = coefficients(:, f);
    trial = trial(isfinite(trial) & abs(trial) > 0);
    n = numel(trial);
    if n < 2
        continue
    end
    unit = trial ./ abs(trial);
    resultant = abs(mean(unit)) ^ 2;
    ppc(f) = (n * resultant - 1) / (n - 1);
end
end
