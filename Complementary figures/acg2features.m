function feature_table = acg2features(logit, lag)
%ACG2FEATURES  Table of autocorrelogram features, one neuron per row.
%
% Intervals use the lag values in seconds. A lag axis that starts at
% 0.05 s does not change PeakInterval, TroughInterval, or OscillationScore.

[values, names] = acg_feature_values(logit, lag);
feature_table = array2table(values, 'VariableNames', names);
end
