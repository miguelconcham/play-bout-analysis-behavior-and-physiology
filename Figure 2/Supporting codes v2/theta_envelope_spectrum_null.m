function result = theta_envelope_spectrum_null(lfp_segments, coefficients, observed, rate, grids, welch, n_perm)
%THETA_ENVELOPE_SPECTRUM_NULL  Phase-randomize the LFP and rebuild the envelope.
%
% Each draw keeps the LFP spectrum, applies the same theta filter, and
% recomputes the envelope spectrum. Excess power is modulation beyond the
% spectrum and the filter. A time shift is not used: it does not change a
% spectrum.

if ~iscell(lfp_segments)
    lfp_segments = {lfp_segments};
end
n_test = numel(grids.test.centers);
n_fine = numel(grids.fine.centers);
null_test = nan(n_perm, n_test);
null_fine = nan(n_perm, n_fine);
for p = 1:n_perm
    surrogates = cell(size(lfp_segments));
    for s = 1:numel(lfp_segments)
        segment = lfp_segments{s};
        segment = segment(isfinite(segment));
        if numel(segment) <= 3 * numel(coefficients)
            surrogates{s} = nan(size(lfp_segments{s}));
            continue
        end
        surrogates{s} = theta_envelope_signal( ...
            theta_envelope_phase_randomize(segment), coefficients);
    end
    spectrum = theta_envelope_session_spectrum(surrogates, rate, grids, welch);
    null_test(p, :) = spectrum.band_power_test;
    null_fine(p, :) = spectrum.band_power_fine;
end
result.test = theta_envelope_surrogate_stats(observed.band_power_test, null_test);
result.fine = theta_envelope_surrogate_stats(observed.band_power_fine, null_fine);
result.n_perm = n_perm;
end
