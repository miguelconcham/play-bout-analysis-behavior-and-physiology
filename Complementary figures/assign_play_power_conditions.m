function conditions = assign_play_power_conditions(play_power, outside_power, options)
%ASSIGN_PLAY_POWER_CONDITIONS  Play tails and power-matched outside periods.
%
% play_power and outside_power are linear delta power, one row per epoch.
% Low and high play are the outer percentiles of play power, the same cut
% as the locking-period analysis. Matching is done on 10*log10 power.
% The caliper is options.caliper_sd times the standard deviation of the
% finite play-epoch decibel values on this channel, and that same caliper
% is used for the whole play set, the low tail, and the high tail.
%
% The low-versus-high comparison uses every play epoch in those tails.
% Each play-versus-no-play comparison uses only the epochs that received
% a partner, so the two sides stay inside the caliper.

if nargin < 3 || isempty(options)
    options = struct();
end
if ~isfield(options, 'low_percent') || isempty(options.low_percent)
    options.low_percent = 25;
end
if ~isfield(options, 'high_percent') || isempty(options.high_percent)
    options.high_percent = 75;
end
if ~isfield(options, 'caliper_sd') || isempty(options.caliper_sd)
    options.caliper_sd = 0.25;
end

play_power = double(play_power(:));
outside_power = double(outside_power(:));
play_db = linear_power_to_db(play_power);
outside_db = linear_power_to_db(outside_power);
tails = play_power_tails(play_power, options.low_percent, options.high_percent);
caliper = power_caliper_db(play_db, options.caliper_sd);

conditions.play_low = find(tails.low);
conditions.play_high = find(tails.high);
conditions.threshold_low = tails.threshold_low;
conditions.threshold_high = tails.threshold_high;
conditions.caliper_db = caliper;
[conditions.play_matched, conditions.outside_matched, gap_all] = ...
    match_db(play_db, outside_db, find(isfinite(play_db)), caliper);
[conditions.play_low_matched, conditions.outside_low_matched, gap_low] = ...
    match_db(play_db, outside_db, find(tails.low & isfinite(play_db)), caliper);
[conditions.play_high_matched, conditions.outside_high_matched, gap_high] = ...
    match_db(play_db, outside_db, find(tails.high & isfinite(play_db)), caliper);
conditions.gap_db = struct( ...
    'play_matched', gap_all, ...
    'play_low_matched', gap_low, ...
    'play_high_matched', gap_high);
end


function db = linear_power_to_db(power)
db = nan(size(power));
ok = isfinite(power) & power > 0;
db(ok) = 10 * log10(power(ok));
end


function caliper = power_caliper_db(play_db, caliper_sd)
values = play_db(isfinite(play_db));
if numel(values) < 2
    scale = 0;
else
    scale = sqrt(sum((values - mean(values)) .^ 2) / (numel(values) - 1));
end
if ~isfinite(scale) || scale == 0
    if isempty(values)
        scale = eps;
    else
        scale = eps(max(1, max(abs(values))));
    end
end
caliper = caliper_sd * scale;
end


function [reference_index, candidate_index, gap] = match_db(play_db, outside_db, reference_rows, caliper)
reference_index = zeros(0, 1);
candidate_index = zeros(0, 1);
gap = NaN;
if isempty(reference_rows) || ~any(isfinite(outside_db))
    return
end
pairs = match_power_periods(play_db(reference_rows), outside_db, caliper);
if isempty(pairs)
    return
end
reference_index = reference_rows(pairs(:, 1));
candidate_index = pairs(:, 2);
gap = median(abs(play_db(reference_index) - outside_db(candidate_index)));
end
