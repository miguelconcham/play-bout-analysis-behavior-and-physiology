function method = canonical_count_method(method)
%CANONICAL_COUNT_METHOD  'pairs', 'fft', or 'crosscorr'.
%
% 'pairs' counts autocorrelogram pairs and divides by T-|lag|. 'fft'
% uses that same correction on counts from the FFT of the binned window.
% 'crosscorr' calls crosscorrelogram. Reference spikes are inset from
% each edge by the lag limit, and the count is divided by that inner
% duration. An empty value is 'pairs'.

if nargin < 1 || isempty(method)
    method = 'pairs';
end
method = lower(strtrim(char(method)));
switch method
    case {'pairs', 'pair', 'autocorrelogram'}
        method = 'pairs';
    case {'fft', 'binned'}
        method = 'fft';
    case {'crosscorr', 'crosscorrelogram', 'cross'}
        method = 'crosscorr';
    otherwise
        error('canonical_count_method:name', ...
            'ACG count must be ''pairs'', ''fft'', or ''crosscorr''.');
end
end
