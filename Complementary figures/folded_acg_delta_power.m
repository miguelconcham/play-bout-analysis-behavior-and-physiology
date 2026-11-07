function power = folded_acg_delta_power(acg, centers, delta_lag)
%FOLDED_ACG_DELTA_POWER  Delta-band power of the half-mirrored autocorrelogram.
%
% Positive and negative lags are averaged. Only lags after 100 ms are
% kept, so the burst does not enter the spectrum. delta_lag is the delta
% period window in seconds, [short long]. The band is 1/long to 1/short
% Hz, which is 2-5 Hz when delta_lag is [0.2 0.5]. Power is the variance
% carried by the one-sided Fourier components inside that band.

if nargin < 3
    delta_lag = [];
end
[~, ~, power] = folded_acg_spectrum(acg, centers, delta_lag);
end
