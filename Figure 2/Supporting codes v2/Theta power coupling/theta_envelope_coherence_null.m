function result = theta_envelope_coherence_null(envelopes, lfps, observed, rate, grids, welch, n_perm)
%THETA_ENVELOPE_COHERENCE_NULL  Segment-shift null for envelope-LFP coherence.
%
% The theta envelope stays fixed. Each Welch segment of the LFP is
% circularly shifted by a new lag, then coherence is recomputed. Shifting
% the whole recording by one lag changes the cross-spectrum phase and leaves
% magnitude-squared coherence in place, so that is not a null for this measure.

n_test = numel(grids.test.centers);
n_fine = numel(grids.fine.centers);
null_test = nan(n_perm, n_test);
null_fine = nan(n_perm, n_fine);
for p = 1:n_perm
    coherence = theta_envelope_session_coherence( ...
        envelopes, lfps, rate, grids, welch, true);
    null_test(p, :) = coherence.band_coherence_test;
    null_fine(p, :) = coherence.band_coherence_fine;
end
result.test = theta_envelope_surrogate_stats(observed.band_coherence_test, null_test);
result.fine = theta_envelope_surrogate_stats(observed.band_coherence_fine, null_fine);
result.n_perm = n_perm;
end
