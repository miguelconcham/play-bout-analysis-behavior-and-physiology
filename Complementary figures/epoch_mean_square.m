function power = epoch_mean_square(values, sample_rate, epochs)
%EPOCH_MEAN_SQUARE  Mean of values squared inside each epoch.
%
% Sample k is at time k/sample_rate, the same clock as the LFP time vector
% in the phase-locking power control. An epoch is [start, end). Pass the
% delta-band envelope and the result is mean(envelope.^2). An epoch with
% no sample is NaN.

power = nan(size(epochs, 1), 1);
if isempty(epochs) || isempty(values) || ~isscalar(sample_rate) || sample_rate <= 0
    return
end
values = double(values(:)) .^ 2;
n = numel(values);
epochs = double(epochs(:, 1:2));
for i = 1:size(epochs, 1)
    if ~all(isfinite(epochs(i, :))) || epochs(i, 2) <= epochs(i, 1)
        continue
    end
    i0 = max(1, ceil(epochs(i, 1) * sample_rate - 1e-9));
    i1 = min(n, floor(epochs(i, 2) * sample_rate - 1e-12));
    if i1 >= i0
        power(i) = mean(values(i0:i1));
    end
end
end
