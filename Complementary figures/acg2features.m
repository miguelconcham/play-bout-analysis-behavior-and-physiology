function ACG_features = acg2features(imagesc_sum_power_high, lags)
%% Input
ACG_matrix = imagesc_sum_power_high;

% Lags in seconds (one value per ACG column)
lags_s = lags(:);

nNeurons = size(ACG_matrix, 1);

assert(size(ACG_matrix, 2) == numel(lags_s), ...
    'Number of ACG columns must match number of lags.');

%% Preallocate
featureNames = {
    'Baseline'
    'Peak1'
    'Trough1'
    'Peak2'
    'Trough2'
    'PeakLag1'
    'TroughLag1'
    'PeakLag2'
    'TroughLag2'
    'Modulation1'
    'Modulation2'
    'Persistence'
    'NormModulation1'
    'NormModulation2'
    'PeakInterval'
    'TroughInterval'
    'PeakFrequency'
    'TroughFrequency'
    'OscillationScore'
    'DominantFrequency'
    'FrequencyResolution'
};

features = nan(nNeurons, numel(featureNames));

%% Analyze each neuron

for i = 1:nNeurons

    acg = ACG_matrix(i, :)';

    R = quantify_acg_rhythmicity(acg, lags_s, ...
        'TroughWindows', [0.05 0.15; 0.35 0.45], ...
        'PeakWindows', [0.25 0.35; 0.55 0.65], ...
        'BaselineWindow', [0.8 1.0], ...
        'FrequencyBand', [2 5]);

    features(i,:) = [
        R.Baseline
        R.PeakValues(1)
        R.TroughValues(1)
        R.PeakValues(2)
        R.TroughValues(2)
        R.PeakLags(1)
        R.TroughLags(1)
        R.PeakLags(2)
        R.TroughLags(2)
        R.FirstModulation
        R.SecondModulation
        R.Persistence
        R.BaselineNormalizedFirstModulation
        R.BaselineNormalizedSecondModulation
        R.PeakInterval
        R.TroughInterval
        R.PeakBasedFrequencyHz
        R.TroughBasedFrequencyHz
        R.OscillationScore
        R.DominantFrequencyHz
        R.FrequencyResolutionHz
    ];

end

%% Create feature table

ACG_features = array2table(features, ...
    'VariableNames', featureNames);

ACG_features.NeuronIndex = (1:nNeurons)';

% Put neuron index first
ACG_features = movevars(ACG_features, ...
    'NeuronIndex', 'Before', 1);
end
