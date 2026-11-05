function detected = detect_broadband_cycles(lfp, sampling_rate, options)
%DETECT_BROADBAND_CYCLES  Trough cycles on a broad low-pass.
%
% Local peaks and troughs are found on the low-pass itself. Consecutive
% extrema of the same type are collapsed to the most extreme sample, so
% the retained sequence alternates. A cycle runs from one trough to the
% next trough when one peak lies between them. Its frequency is
% 1 / period. A cycle is accepted when its frequency is inside
% frequency_limits_hz and it meets the amplitude, consistency, and
% monotonicity gates.

if nargin < 3 || isempty(options)
    options = struct();
end
options = broadband_cycle_options(options);
if isempty(options.lowpass_coefficients)
    options.lowpass_coefficients = broadband_lowpass_coefficients(options);
end

lfp = double(lfp(:).');
lfp = lfp - mean(lfp, 'omitnan');
if any(~isfinite(lfp))
    error('detect_broadband_cycles:lfp', 'LFP contains non-finite samples.');
end
order = numel(options.lowpass_coefficients) - 1;
if numel(lfp) <= 3 * order
    error('detect_broadband_cycles:length', ...
        'The interval is shorter than three times the low-pass order.');
end

broad = filtfilt(options.lowpass_coefficients, 1, lfp);
edge = order;
usable = false(size(broad));
usable(edge + 1:numel(broad) - edge) = true;
extrema = alternating_extrema(broad);
extrema = extrema(usable(extrema));

frequency = zeros(0, 1);
amplitude = zeros(0, 1);
rdsym = zeros(0, 1);
monotonicity = zeros(0, 1);
period_consistency = zeros(0, 1);
amplitude_consistency = zeros(0, 1);
trough_start = zeros(0, 1);
peak_sample = zeros(0, 1);
trough_end = zeros(0, 1);

is_trough = false(size(extrema));
for extremum_n = 1:numel(extrema)
    sample = extrema(extremum_n);
    is_trough(extremum_n) = broad(sample) <= broad(max(1, sample - 1)) && ...
        broad(sample) <= broad(min(numel(broad), sample + 1));
end

for trough_n = 1:numel(extrema) - 2
    if ~is_trough(trough_n) || ~is_trough(trough_n + 2) || is_trough(trough_n + 1)
        continue
    end
    start_sample = extrema(trough_n);
    peak = extrema(trough_n + 1);
    stop_sample = extrema(trough_n + 2);
    period_s = (stop_sample - start_sample) / sampling_rate;
    if period_s <= 0
        continue
    end
    frequency(end + 1, 1) = 1 / period_s; %#ok<AGROW>
    rise = broad(peak) - broad(start_sample);
    decay = broad(peak) - broad(stop_sample);
    amplitude(end + 1, 1) = (rise + decay) / 2; %#ok<AGROW>
    rdsym(end + 1, 1) = ((peak - start_sample) / sampling_rate) / period_s; %#ok<AGROW>
    monotonicity(end + 1, 1) = flank_monotonicity( ...
        broad, start_sample, peak, stop_sample); %#ok<AGROW>
    amplitude_consistency(end + 1, 1) = min(rise, decay) / max(rise, decay); %#ok<AGROW>
    period_consistency(end + 1, 1) = NaN; %#ok<AGROW>
    trough_start(end + 1, 1) = start_sample; %#ok<AGROW>
    peak_sample(end + 1, 1) = peak; %#ok<AGROW>
    trough_end(end + 1, 1) = stop_sample; %#ok<AGROW>
end

for cycle_n = 1:numel(frequency)
    ratios = [];
    if cycle_n > 1
        ratios(end + 1) = adjacent_ratio(frequency(cycle_n), frequency(cycle_n - 1)); %#ok<AGROW>
    end
    if cycle_n < numel(frequency)
        ratios(end + 1) = adjacent_ratio(frequency(cycle_n), frequency(cycle_n + 1)); %#ok<AGROW>
    end
    if ~isempty(ratios)
        period_consistency(cycle_n) = min(ratios);
    end
end

if isempty(amplitude)
    cutoff = NaN;
else
    cutoff = linear_percentile(amplitude, options.amplitude_percentile);
end
in_band = frequency >= options.frequency_limits_hz(1) & ...
    frequency <= options.frequency_limits_hz(2);
accepted = in_band & amplitude >= cutoff & ...
    amplitude_consistency >= options.amplitude_consistency_threshold & ...
    period_consistency >= options.period_consistency_threshold & ...
    monotonicity >= options.monotonicity_threshold;

detected = struct();
detected.broad = broad;
detected.frequency_hz = frequency;
detected.amplitude = amplitude;
detected.rdsym = rdsym;
detected.monotonicity = monotonicity;
detected.period_consistency = period_consistency;
detected.amplitude_consistency = amplitude_consistency;
detected.trough_start = trough_start;
detected.peak_sample = peak_sample;
detected.trough_end = trough_end;
detected.accepted = accepted;
detected.amplitude_cutoff = cutoff;
end

function extrema = alternating_extrema(signal)
n = numel(signal);
is_peak = false(1, n);
is_trough = false(1, n);
is_peak(2:n - 1) = signal(2:n - 1) > signal(1:n - 2) & ...
    signal(2:n - 1) >= signal(3:n);
is_trough(2:n - 1) = signal(2:n - 1) < signal(1:n - 2) & ...
    signal(2:n - 1) <= signal(3:n);
indexes = find(is_peak | is_trough);
extrema = zeros(0, 1);
cursor = 1;
while cursor <= numel(indexes)
    stop = cursor;
    this_is_peak = is_peak(indexes(cursor));
    while stop < numel(indexes) && is_peak(indexes(stop + 1)) == this_is_peak
        stop = stop + 1;
    end
    group = indexes(cursor:stop);
    if this_is_peak
        [~, relative] = max(signal(group));
    else
        [~, relative] = min(signal(group));
    end
    extrema(end + 1, 1) = group(relative); %#ok<AGROW>
    cursor = stop + 1;
end
end

function fraction = flank_monotonicity(signal, trough_start, peak, trough_end)
rise = diff(signal(trough_start:peak));
decay = diff(signal(peak:trough_end));
steps = [rise(:); decay(:)];
if isempty(steps)
    fraction = NaN;
    return
end
correct = [rise(:) > 0; decay(:) < 0];
fraction = mean(correct);
end

function ratio = adjacent_ratio(frequency_a, frequency_b)
period_a = 1 / frequency_a;
period_b = 1 / frequency_b;
ratio = min(period_a, period_b) / max(period_a, period_b);
end
