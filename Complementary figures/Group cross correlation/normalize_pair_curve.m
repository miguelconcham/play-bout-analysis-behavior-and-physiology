function [ratio, normalized] = normalize_pair_curve(conditional_rate, target_rate)
%NORMALIZE_PAIR_CURVE  Put one session correlogram on a rate-free scale.
%
% conditional_rate is target spikes per second around a reference spike.
% target_rate is that same target group's own rate in the session.
% ratio is conditional_rate divided by target_rate. One is independence,
% whatever the firing rate and the session length were. normalized is
% ratio minus 1, so zero is independence and a negative value is the
% fraction by which the target fell below its own rate. Session curves
% are normalized before they are averaged, and the confidence interval
% is the spread of these normalized values.

conditional_rate = double(conditional_rate(:))';
ratio = nan(size(conditional_rate));
normalized = ratio;
if ~(isscalar(target_rate) && isfinite(double(target_rate)) && double(target_rate) > 0)
    return
end
ratio = conditional_rate / double(target_rate);
ratio(~isfinite(conditional_rate)) = NaN;
normalized = ratio - 1;
end
