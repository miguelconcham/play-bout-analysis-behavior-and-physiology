function result = theta_envelope_onset_null(envelope, rate_out, onsets, window, frequencies, confirmatory, intervals, rate_in, null_options)
%THETA_ENVELOPE_ONSET_NULL  Pre-onset Fourier coefficients and an onset-shift null.
%
% Onsets that share an interval are shifted by one common lag, so the gaps
% between behaviors stay as they were. The lag is at least min_shift_seconds
% when the interval is long enough. Windows that land on non-finite samples
% are dropped for that draw.
%
% A common lag rotates a stationary sinusoid without changing the spread of
% phases across behaviors, so this null is for a pattern that is local to
% onset. A rhythm that runs through the whole session is the session test.
%
% When compute_trial_summaries is true, the result also holds the
% leave-one-out correlation of the pre-onset waveform and the fraction of
% behaviors whose confirmatory-band phase lies within a quarter cycle of the
% other behaviors in the recording.

if nargin < 8 || isempty(rate_in)
    rate_in = rate_out;
end
onsets = double(onsets(:));
frequencies = frequencies(:).';
confirmatory = logical(confirmatory(:).');
[epochs, time] = theta_envelope_epochs(envelope, rate_out, onsets, window, rate_in);
coefficient = theta_envelope_fourier(epochs, time, frequencies);

result.window = window(:).';
result.time = time(:).';
result.coefficient = coefficient;
result.n_valid = sum(all(isfinite(epochs), 2));
result.mean_coefficient = mean_coefficient(coefficient);
result.ppc = theta_envelope_phase_consistency(coefficient);
result.loo_r = [];
result.loo_mean_r = NaN;
result.phase_fraction = NaN;
result.null_mean_coefficient = zeros(0, numel(frequencies));
result.null_ppc = zeros(0, numel(frequencies));
result.null_loo_mean_r = zeros(0, 1);
result.null_phase_fraction = zeros(0, 1);

do_trials = isfield(null_options, 'compute_trial_summaries') && ...
    null_options.compute_trial_summaries;
if do_trials
    result.loo_r = theta_envelope_leave_one_out(epochs);
    finite_loo = result.loo_r(isfinite(result.loo_r));
    if ~isempty(finite_loo)
        result.loo_mean_r = mean(finite_loo);
    end
    result.phase_fraction = locked_phase_fraction(coefficient, confirmatory);
end

n_perm = null_options.n_perm;
if n_perm < 1
    return
end
result.null_mean_coefficient = nan(n_perm, numel(frequencies));
result.null_ppc = nan(n_perm, numel(frequencies));
if do_trials
    result.null_loo_mean_r = nan(n_perm, 1);
    result.null_phase_fraction = nan(n_perm, 1);
end
if isempty(intervals)
    intervals = default_interval(envelope, rate_out, rate_in);
end
for p = 1:n_perm
    shifted = shift_onsets(onsets, intervals, null_options.min_shift_seconds);
    shifted_epochs = theta_envelope_epochs( ...
        envelope, rate_out, shifted, window, rate_in);
    shifted_coefficient = theta_envelope_fourier(shifted_epochs, time, frequencies);
    result.null_mean_coefficient(p, :) = mean_coefficient(shifted_coefficient);
    result.null_ppc(p, :) = theta_envelope_phase_consistency(shifted_coefficient);
    if do_trials
        loo = theta_envelope_leave_one_out(shifted_epochs);
        finite_loo = loo(isfinite(loo));
        if ~isempty(finite_loo)
            result.null_loo_mean_r(p) = mean(finite_loo);
        end
        result.null_phase_fraction(p) = locked_phase_fraction( ...
            shifted_coefficient, confirmatory);
    end
end
end

function shifted = shift_onsets(onsets, intervals, min_shift)
shifted = onsets;
for interval_n = 1:size(intervals, 1)
    t0 = intervals(interval_n, 1);
    t1 = intervals(interval_n, 2);
    span = t1 - t0;
    inside = onsets >= t0 & onsets < t1;
    if ~any(inside) || span <= 0
        continue
    end
    usable = min_shift;
    if span <= 2 * usable
        usable = span / 4;
    end
    if span <= 2 * usable
        continue
    end
    lag = usable + rand * (span - 2 * usable);
    shifted(inside) = t0 + mod((onsets(inside) - t0) + lag, span);
end
end

function intervals = default_interval(envelope, rate_out, rate_in)
if rate_in == rate_out
    t1 = numel(envelope) / rate_out;
else
    factor = rate_in / rate_out;
    t1 = ((numel(envelope) - 1) * factor + 1) / rate_in;
end
intervals = [1 / rate_in, t1];
end

function value = mean_coefficient(coefficient)
value = nan(1, size(coefficient, 2));
for f = 1:size(coefficient, 2)
    trial = coefficient(:, f);
    trial = trial(isfinite(trial));
    if ~isempty(trial)
        value(f) = mean(trial);
    end
end
end

function fraction = locked_phase_fraction(coefficient, mask)
fraction = NaN;
if isempty(mask) || ~any(mask) || isempty(coefficient)
    return
end
selected = coefficient(:, mask);
trial = nan(size(selected, 1), 1);
for row = 1:size(selected, 1)
    values = selected(row, :);
    values = values(isfinite(values));
    if ~isempty(values)
        trial(row) = mean(values);
    end
end
trial = trial(isfinite(trial) & abs(trial) > 0);
n = numel(trial);
if n < 2
    return
end
total = sum(trial);
count = 0;
for row = 1:n
    others = total - trial(row);
    difference = angle(trial(row) * conj(others));
    if abs(difference) <= pi / 2
        count = count + 1;
    end
end
fraction = count / n;
end
