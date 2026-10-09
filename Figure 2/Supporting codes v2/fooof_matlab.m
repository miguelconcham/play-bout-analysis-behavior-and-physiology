function out = fooof_matlab(f, P, varargin)
% FOOOF_MATLAB
%
% MATLAB implementation closely following FOOOF 1.1.1.
%
% Reference:
% Donoghue et al. (2020). Parameterizing neural power spectra into
% periodic and aperiodic components. Nature Neuroscience.
% DOI: 10.1038/s41593-020-00744-x
%
% This implementation follows the algorithm and parameter conventions
% of FOOOF 1.1.1, while using MATLAB's LSQCURVEFIT as the nonlinear
% least-squares optimizer.
%
% INPUT
%   f : frequency vector in Hz, LINEAR scale
%   P : power vector, LINEAR scale
%
% IMPORTANT:
%   FOOOF expects linear power as input. Internally the spectrum is
%   converted to log10(P).
%
% NAME-VALUE OPTIONS
%
%   'PeakWidthLimits' : [lower upper] Hz
%                       Default: [0.5 12]
%
%   'MaxNPeaks'       : maximum number of peaks
%                       Default: Inf
%
%   'MinPeakHeight'   : minimum peak height above aperiodic component
%                       in log10 power
%                       Default: 0
%
%   'PeakThreshold'   : relative peak threshold in SD
%                       Default: 2
%
%   'AperiodicMode'   : 'fixed' or 'knee'
%                       Default: 'fixed'
%
% OUTPUT
%
%   out.aperiodic_params
%       fixed: [offset exponent]
%       knee : [offset knee exponent]
%
%   out.gaussian_params
%       [center frequency, height, standard deviation]
%
%   out.peak_params
%       [CF, PW, BW]
%
%   out.aperiodic_fit
%   out.peak_fit
%   out.full_fit
%   out.flat_spectrum
%   out.peak_removed_spectrum
%
%   out.r_squared
%   out.error
%
%   out.settings
%
% -------------------------------------------------------------------------

%% Parse options
p = inputParser;
addParameter(p, 'PeakWidthLimits', [0.5 12]);
addParameter(p, 'MaxNPeaks', Inf);
addParameter(p, 'MinPeakHeight', 0);
addParameter(p, 'PeakThreshold', 2);
addParameter(p, 'AperiodicMode', 'fixed');
parse(p, varargin{:});

peakWidthLimits = p.Results.PeakWidthLimits;
maxNPeaks       = p.Results.MaxNPeaks;
minPeakHeight   = p.Results.MinPeakHeight;
peakThreshold   = p.Results.PeakThreshold;
aperiodicMode   = lower(p.Results.AperiodicMode);

%% Validate input
f = f(:);
P = P(:);

if numel(f) ~= numel(P)
    error('f and P must have the same number of elements.');
end

if any(~isfinite(f)) || any(~isfinite(P))
    error('f and P must not contain NaN or Inf.');
end

if any(f <= 0)
    error('Frequencies must be > 0.');
end

if any(P <= 0)
    error('Power must be > 0.');
end

% Sort frequency vector
[f, order] = sort(f);
P = P(order);

%% FOOOF requires evenly spaced frequencies
freqDiff = diff(f);
freqRes  = freqDiff(1);

if ~all(abs(freqDiff - freqRes) < 1e-10 * max(1, abs(freqRes)))
    error(['FOOOF requires evenly spaced frequencies. ' ...
           'Interpolate/resample the spectrum before fitting.']);
end

freqRange = [f(1), f(end)];

%% Convert power to log10 scale
logP = log10(P);

%% Settings corresponding to FOOOF internals
% FOOOF uses 0.025 as percentile value.
% MATLAB prctile uses 0-100.
apPercentileThresh = 0.025;

% Peak width limits are specified as TWO-SIDED bandwidth.
% Internally FOOOF uses Gaussian SD = BW / 2.
gaussStdLimits = peakWidthLimits / 2;

% Gaussian candidate edge threshold in units of SD
bwStdEdge = 1.0;

% Gaussian overlap threshold in units of SD
gaussOverlapThresh = 0.75;

% Center-frequency fitting bound
cfBound = 1.5;

% Maximum optimizer evaluations
maxFEv = 5000;

%% ================================================================
% STEP 1
% Initial robust aperiodic fit
% ================================================================
initialAperiodicParams = robust_aperiodic_fit( ...
    f, logP, aperiodicMode, ...
    apPercentileThresh, maxFEv);

initialAperiodicFit = aperiodic_function( ...
    f, initialAperiodicParams, aperiodicMode);

%% ================================================================
% STEP 2
% Flatten spectrum
% ================================================================
flatSpectrum = logP - initialAperiodicFit;

%% ================================================================
% STEP 3
% Iterative peak detection
% ================================================================
[gaussianGuesses, flatIter] = iterative_peak_search( ...
    f, ...
    flatSpectrum, ...
    freqRes, ...
    freqRange, ...
    peakWidthLimits, ...
    gaussStdLimits, ...
    maxNPeaks, ...
    minPeakHeight, ...
    peakThreshold); %#ok<ASGLU>

%% ================================================================
% Remove peaks too close to spectral edges
% and peaks that overlap too strongly
% ================================================================
gaussianGuesses = drop_edge_peaks( ...
    gaussianGuesses, ...
    freqRange, ...
    bwStdEdge);

gaussianGuesses = drop_overlapping_peaks( ...
    gaussianGuesses, ...
    gaussOverlapThresh);

%% ================================================================
% STEP 4
% Jointly fit all Gaussian peaks
% ================================================================
if isempty(gaussianGuesses)
    gaussianParams = zeros(0, 3);
else
    gaussianParams = fit_gaussian_peaks( ...
        f, ...
        flatSpectrum, ...
        gaussianGuesses, ...
        gaussStdLimits, ...
        cfBound, ...
        freqRange, ...
        maxFEv);

    % Sort by center frequency
    gaussianParams = sortrows(gaussianParams, 1);
end

%% ================================================================
% Construct periodic fit
% ================================================================
peakFit = zeros(size(f));
for ii = 1:size(gaussianParams, 1)
    peakFit = peakFit + gaussian_function( ...
        f, ...
        gaussianParams(ii, 1), ...
        gaussianParams(ii, 2), ...
        gaussianParams(ii, 3));
end

%% ================================================================
% STEP 5
% Remove periodic component
% ================================================================
peakRemovedSpectrum = logP - peakFit;

%% ================================================================
% STEP 6
% Final aperiodic fit
%
% IMPORTANT:
% This is the SIMPLE aperiodic fit in FOOOF, not the robust fit.
% ================================================================
aperiodicParams = simple_aperiodic_fit( ...
    f, ...
    peakRemovedSpectrum, ...
    aperiodicMode, ...
    maxFEv);

aperiodicFit = aperiodic_function( ...
    f, ...
    aperiodicParams, ...
    aperiodicMode);

%% ================================================================
% STEP 7
% Full model
% ================================================================
fullFit = aperiodicFit + peakFit;

%% ================================================================
% Convert Gaussian parameters to FOOOF peak parameters
%
% FOOOF:
%
%   CF = Gaussian center
%   PW = height of full model above aperiodic at nearest CF
%   BW = 2 * Gaussian SD
% ================================================================
peakParams = zeros(size(gaussianParams));
for ii = 1:size(gaussianParams, 1)
    cf = gaussianParams(ii, 1);
    sd = gaussianParams(ii, 3);

    [~, ind] = min(abs(f - cf));
    pw = fullFit(ind) - aperiodicFit(ind);
    bw = 2 * sd;

    peakParams(ii, :) = [cf pw bw];
end

%% ================================================================
% Goodness of fit
% ================================================================
R = corrcoef(logP, fullFit);
if numel(R) >= 4
    R2 = R(1, 2)^2;
else
    R2 = NaN;
end

MAE = mean(abs(logP - fullFit));

%% ================================================================
% Output
% ================================================================
out.f                        = f;
out.power                    = P;
out.log_power                = logP;
out.freq_range               = freqRange;
out.freq_res                 = freqRes;

out.initial_aperiodic_params = initialAperiodicParams;
out.initial_aperiodic_fit    = initialAperiodicFit;

out.flat_spectrum            = flatSpectrum;
out.gaussian_guesses         = gaussianGuesses;
out.gaussian_params          = gaussianParams;

out.peak_fit                 = peakFit;
out.peak_removed_spectrum    = peakRemovedSpectrum;

out.aperiodic_params         = aperiodicParams;
out.aperiodic_fit            = aperiodicFit;

out.full_fit                 = fullFit;
out.peak_params              = peakParams;

out.r_squared                = R2;
out.error                    = MAE;

out.settings.PeakWidthLimits = peakWidthLimits;
out.settings.MaxNPeaks       = maxNPeaks;
out.settings.MinPeakHeight   = minPeakHeight;
out.settings.PeakThreshold   = peakThreshold;
out.settings.AperiodicMode   = aperiodicMode;
out.settings.FOOOFVersion    = '1.1.1-compatible';

end

%% ========================================================================
% APERIODIC FUNCTION
% ========================================================================
function y = aperiodic_function(f, params, mode)

switch lower(mode)
    case 'fixed'
        % FOOOF fixed aperiodic:
        %
        %   L(F) = b - chi * log10(F)
        offset   = params(1);
        exponent = params(2);
        y = offset - exponent .* log10(f);

    case 'knee'
        % FOOOF knee aperiodic:
        %
        %   L(F) = b - log10(k + F^chi)
        offset   = params(1);
        knee     = params(2);
        exponent = params(3);
        y = offset - log10(knee + f.^exponent);

    otherwise
        error('AperiodicMode must be ''fixed'' or ''knee''.');
end

end

%% ========================================================================
% GAUSSIAN FUNCTION
% ========================================================================
function y = gaussian_function(f, center, height, sigma)
% FOOOF Gaussian:
%
%   G(F) = a * exp(-(F-c)^2 / (2*w^2))
%
%   center = c
%   height = a
%   sigma  = w

y = height .* exp( ...
    -(f - center).^2 ./ (2 * sigma.^2));

end

%% ========================================================================
% ROBUST APERIODIC FIT
% ========================================================================
function params = robust_aperiodic_fit( ...
    f, logP, mode, percentileThreshold, maxFEv)

% ------------------------------------------------------------
% First quick aperiodic fit
% ------------------------------------------------------------
p0 = simple_aperiodic_fit( ...
    f, logP, mode, maxFEv);

initialFit = aperiodic_function( ...
    f, p0, mode);

% ------------------------------------------------------------
% Flatten
% ------------------------------------------------------------
flatSpec = logP - initialFit;

% FOOOF sets values below zero to zero before percentile
% selection.
flatSpec(flatSpec < 0) = 0;

% ------------------------------------------------------------
% Select low percentile points
%
% IMPORTANT:
% FOOOF internal setting = 0.025
% This means the 0.025th percentile, not 2.5th percentile.
% ------------------------------------------------------------
threshold = prctile( ...
    flatSpec, ...
    percentileThreshold);

mask = flatSpec <= threshold;

fIgnore = f(mask);
pIgnore = logP(mask);

% ------------------------------------------------------------
% Refit
% ------------------------------------------------------------
params = fit_aperiodic( ...
    fIgnore, ...
    pIgnore, ...
    p0, ...
    mode, ...
    maxFEv);

end

%% ========================================================================
% SIMPLE APERIODIC FIT
% ========================================================================
function params = simple_aperiodic_fit( ...
    f, logP, mode, maxFEv)

% FOOOF initial guesses:
%
%   offset = first power value
%
%   exponent =
%       abs((last power - first power) /
%           (log10(last frequency) - log10(first frequency)))
%
%   knee = 0

offsetGuess = logP(1);

exponentGuess = abs( ...
    (logP(end) - logP(1)) ./ ...
    (log10(f(end)) - log10(f(1))));

switch lower(mode)
    case 'fixed'
        p0 = [offsetGuess exponentGuess];

    case 'knee'
        kneeGuess = 0;
        p0 = [ ...
            offsetGuess ...
            kneeGuess ...
            exponentGuess];

    otherwise
        error('Unknown aperiodic mode.');
end

params = fit_aperiodic( ...
    f, logP, p0, mode, maxFEv);

end

%% ========================================================================
% FIT APERIODIC
% ========================================================================
function params = fit_aperiodic( ...
    f, logP, p0, mode, maxFEv)

switch lower(mode)
    case 'fixed'
        model = @(p, x) ...
            p(1) - p(2) .* log10(x);

        % FOOOF itself uses unbounded parameters.
        lb = [-Inf -Inf];
        ub = [ Inf  Inf];

    case 'knee'
        model = @(p, x) ...
            p(1) - log10(p(2) + x.^p(3));

        % FOOOF source uses unbounded parameters by default.
        %
        % This is important: do NOT impose arbitrary positive
        % constraints here if you want FOOOF-like behavior.
        lb = [-Inf -Inf -Inf];
        ub = [ Inf  Inf  Inf];

    otherwise
        error('Unknown aperiodic mode.');
end

options = optimoptions( ...
    'lsqcurvefit', ...
    'Display', 'off', ...
    'MaxFunctionEvaluations', maxFEv);

% Suppress warnings caused by temporary invalid logarithms during
% optimization, analogous to FOOOF suppressing RuntimeWarnings.
warningState = warning;
warning('off', 'all');

try
    params = lsqcurvefit( ...
        model, ...
        p0, ...
        f, ...
        logP, ...
        lb, ...
        ub, ...
        options);
catch ME
    warning(warningState);
    rethrow(ME);
end

warning(warningState);

end

%% ========================================================================
% ITERATIVE PEAK SEARCH
% ========================================================================
function [guess, flatIter] = iterative_peak_search( ...
    f, ...
    flatSpectrum, ...
    freqRes, ...
    ~, ...
    peakWidthLimits, ...
    gaussStdLimits, ...
    maxNPeaks, ...
    minPeakHeight, ...
    peakThreshold)

guess = zeros(0, 3);
flatIter = flatSpectrum;

while size(guess, 1) < maxNPeaks

    %% Find maximum of current flattened spectrum
    [maxHeight, maxInd] = max(flatIter);

    %% Relative threshold
    if maxHeight <= peakThreshold * std(flatIter)
        break;
    end

    %% Candidate peak
    guessFreq   = f(maxInd);
    guessHeight = maxHeight;

    %% Absolute threshold
    if ~(guessHeight > minPeakHeight)
        break;
    end

    %% Estimate Gaussian SD from half height
    halfHeight = 0.5 * maxHeight;

    % FOOOF searches LEFT from maxInd-1 down to index 2
    leInd = [];
    for ind = maxInd-1:-1:2
        if flatIter(ind) <= halfHeight
            leInd = ind;
            break;
        end
    end

    % FOOOF searches RIGHT from maxInd+1
    riInd = [];
    for ind = maxInd+1:length(flatIter)
        if flatIter(ind) <= halfHeight
            riInd = ind;
            break;
        end
    end

    %% Determine shortest side
    sides = [];
    if ~isempty(leInd)
        sides(end+1) = abs(leInd - maxInd); %#ok<AGROW>
    end
    if ~isempty(riInd)
        sides(end+1) = abs(riInd - maxInd); %#ok<AGROW>
    end

    if isempty(sides)
        % FOOOF fallback
        guessStd = mean(peakWidthLimits);
    else
        shortSide = min(sides);

        % FOOOF:
        %
        %   FWHM = shortest_side * 2 * freq_res
        fwhm = shortSide * 2 * freqRes;

        % Gaussian:
        %
        %   SD = FWHM / (2*sqrt(2*ln(2)))
        guessStd = fwhm / ...
            (2 * sqrt(2 * log(2)));
    end

    %% Restrict guess to allowed SD limits
    guessStd = max( ...
        guessStd, ...
        gaussStdLimits(1));

    guessStd = min( ...
        guessStd, ...
        gaussStdLimits(2));

    %% Save candidate
    guess(end+1, :) = [ ...
        guessFreq ...
        guessHeight ...
        guessStd]; %#ok<AGROW>

    %% Subtract candidate Gaussian
    peakGaussian = gaussian_function( ...
        f, ...
        guessFreq, ...
        guessHeight, ...
        guessStd);

    flatIter = flatIter - peakGaussian;
end

end

%% ========================================================================
% DROP EDGE PEAKS
% ========================================================================
function guess = drop_edge_peaks( ...
    guess, freqRange, bwStdEdge)

if isempty(guess)
    return;
end

cf      = guess(:, 1);
stdVals = guess(:, 3);

keep = ...
    abs(cf - freqRange(1)) > stdVals * bwStdEdge & ...
    abs(cf - freqRange(2)) > stdVals * bwStdEdge;

guess = guess(keep, :);

end

%% ========================================================================
% DROP OVERLAPPING PEAKS
% ========================================================================
function guess = drop_overlapping_peaks( ...
    guess, overlapThreshold)

if size(guess, 1) < 2
    return;
end

% Sort by center frequency
guess = sortrows(guess, 1);

bounds = zeros(size(guess, 1), 2);
for ii = 1:size(guess, 1)
    bounds(ii, :) = [ ...
        guess(ii, 1) - guess(ii, 3) * overlapThreshold ...
        guess(ii, 1) + guess(ii, 3) * overlapThreshold];
end

drop = false(size(guess, 1), 1);

for ii = 1:size(bounds, 1)-1
    leftUpper  = bounds(ii, 2);
    rightLower = bounds(ii+1, 1);

    if leftUpper > rightLower
        % Drop lower-height Gaussian
        if guess(ii, 2) < guess(ii+1, 2)
            drop(ii) = true;
        else
            drop(ii+1) = true;
        end
    end
end

guess = guess(~drop, :);

end

%% ========================================================================
% JOINT GAUSSIAN FIT
% ========================================================================
function params = fit_gaussian_peaks( ...
    f, ...
    flatSpectrum, ...
    guess, ...
    gaussStdLimits, ...
    cfBound, ...
    freqRange, ...
    maxFEv)

nPeaks = size(guess, 1);

%% Initial parameter vector
p0 = reshape(guess.', [], 1);

%% Construct parameter bounds
lb = zeros(size(p0));
ub = zeros(size(p0));

for ii = 1:nPeaks
    ind = (ii-1)*3 + (1:3);

    center = guess(ii, 1);
    sigma  = guess(ii, 3);

    % FOOOF:
    %
    %   center bounds:
    %
    %       center +/- 2 * cf_bound * sigma
    %
    %   cf_bound = 1.5
    %
    %   => center +/- 3 sigma
    cfLow  = center - ...
        2 * cfBound * sigma;
    cfHigh = center + ...
        2 * cfBound * sigma;

    % Restrict to spectrum frequency range
    cfLow  = max(cfLow,  freqRange(1));
    cfHigh = min(cfHigh, freqRange(2));

    lb(ind) = [ ...
        cfLow ...
        0 ...
        gaussStdLimits(1)];

    ub(ind) = [ ...
        cfHigh ...
        Inf ...
        gaussStdLimits(2)];
end

%% Joint Gaussian model
model = @(p, x) ...
    sum_gaussians(p, x, nPeaks);

%% Fit
options = optimoptions( ...
    'lsqcurvefit', ...
    'Display', 'off', ...
    'MaxFunctionEvaluations', maxFEv);

warningState = warning;
warning('off', 'all');

try
    pFit = lsqcurvefit( ...
        model, ...
        p0, ...
        f, ...
        flatSpectrum, ...
        lb, ...
        ub, ...
        options);
catch ME
    warning(warningState);
    rethrow(ME);
end

warning(warningState);

%% Reorganize
params = reshape(pFit, 3, []).';

end

%% ========================================================================
% SUM OF GAUSSIANS
% ========================================================================
function y = sum_gaussians(p, f, nPeaks)

y = zeros(size(f));

for ii = 1:nPeaks
    ind = (ii-1)*3 + (1:3);

    center = p(ind(1));
    height = p(ind(2));
    sigma  = p(ind(3));

    y = y + gaussian_function( ...
        f, center, height, sigma);
end

end
