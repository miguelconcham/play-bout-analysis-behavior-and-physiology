function index = theta_envelope_time_index(time_s, rate_out, rate_in)
%THETA_ENVELOPE_TIME_INDEX  Decimated sample nearest a time in seconds.
%
% Sample k of the LFP is at time k / sampling_rate. After integer
% decimation that keeps sample 1, (k-1)*factor+1 of the original trace,
% the decimated index follows that same clock.

if nargin < 3 || isempty(rate_in)
    rate_in = rate_out;
end
time_s = double(time_s);
if rate_in == rate_out
    index = round(time_s .* rate_out);
    return
end
factor = rate_in / rate_out;
if factor ~= round(factor)
    error('theta_envelope_time_index:rate', ...
        'rate_in must be an integer multiple of rate_out.');
end
factor = round(factor);
original_index = round(time_s .* rate_in);
index = round((original_index - 1) ./ factor) + 1;
end
