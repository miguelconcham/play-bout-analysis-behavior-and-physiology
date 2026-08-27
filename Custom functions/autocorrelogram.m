function all_lags = autocorrelogram(x, range_limits)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if size(x,2)>1
    x = x';
end
all_lags = [];

for j=1:numel(x)
    current_spike = x(j);
    lags_in_range = x(x>=current_spike+range_limits(1) & x<=current_spike+range_limits(2))-current_spike;
    lags_in_range(lags_in_range==0) = [];
    all_lags = [all_lags;lags_in_range];
end

end