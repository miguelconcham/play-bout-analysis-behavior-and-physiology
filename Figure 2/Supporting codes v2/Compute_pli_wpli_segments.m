function pli_struct = Compute_pli_wpli_segments(...
    PAG_LFP, Hd, n_samples, sample_length, pctls, varargin)
% Compute PLI and wPLI across channel pairs using multiple segments.
%
% OPTIONAL NAME-VALUE:
%   'Overlap'     : fraction overlap to force [0..0.99]. Default auto-overlap.
%   'Selection'   : 'random' (default) or 'sequential'
%   'PowerThresh' : k for artifact rejection (median + k*MAD). [] disables. Default = 5.
%   'Verbose'     : true/false (default true)

%% Parse inputs
p = inputParser;
addParameter(p, 'Overlap', [], @(x) isempty(x) || (isnumeric(x) && x >= 0 && x < 1));
addParameter(p, 'Selection', 'random', @(s) any(strcmpi(s, {'random', 'sequential'})));
addParameter(p, 'PowerThresh', 5, @(x) isempty(x) || (isnumeric(x) && x >= 0));
addParameter(p, 'Verbose', true, @islogical);
parse(p, varargin{:});
userOverlap = p.Results.Overlap;
selection = lower(p.Results.Selection);
power_k = p.Results.PowerThresh;
verbose = p.Results.Verbose;

[nCh, nTime] = size(PAG_LFP);
if sample_length > nTime
    error('sample_length (%d) > total recording length (%d).', sample_length, nTime);
end

%% Filter and compute analytic signals
if verbose
    fprintf('Filtering entire recording (zero-phase FIR)...\n');
end
filtered = filtfilt(Hd.Coefficients, 1, double(PAG_LFP') )';
analytic = hilbert(filtered')';

%% Compute segment starts
n_max_nonoverlap = floor(nTime / sample_length);

if ~isempty(userOverlap)
    ov = userOverlap;
    step = max(1, round(sample_length * (1 - ov)));
    starts_possible = 1:step:(nTime - sample_length + 1);
    if isempty(starts_possible)
        error('No possible windows with the specified Overlap and sample_length.');
    end
    n_possible = numel(starts_possible);
    if strcmp(selection, 'random')
        if n_samples <= n_possible
            rand_idx = randperm(n_possible, n_samples);
            starts = starts_possible(rand_idx);
        else
            perm = randperm(n_possible);
            starts = starts_possible(perm);
            extra = n_samples - n_possible;
            starts = [starts, starts_possible(randi(n_possible, 1, extra))];
            starts = starts(randperm(numel(starts)));
        end
    else
        if n_samples <= n_possible
            starts = starts_possible(1:n_samples);
        else
            reps = ceil(n_samples / n_possible);
            starts = repmat(starts_possible, 1, reps);
            starts = starts(1:n_samples);
        end
    end
    mean_overlap = 1 - (step / sample_length);
else
    if n_samples <= n_max_nonoverlap
        step = sample_length;
        starts = 1:step:(nTime - sample_length + 1);
        if numel(starts) > n_samples
            starts = starts(1:n_samples);
        end
        mean_overlap = 0;
        if numel(starts) > n_samples
            if strcmp(selection, 'random')
                starts = starts(randperm(numel(starts), n_samples));
            else
                starts = starts(1:n_samples);
            end
        end
    else
        step = floor((nTime - sample_length) / (n_samples - 1));
        if step < 1
            step = 1;
        end
        starts = 1:step:(1 + step * (n_samples - 1));
        if numel(starts) >= n_samples
            starts = starts(1:n_samples);
        else
            while numel(starts) < n_samples
                starts(end + 1) = starts(end);
            end
        end
        mean_overlap = 1 - (step / sample_length);
    end
end

%% Optional artifact rejection
if ~isempty(power_k) && power_k > 0
    if verbose, fprintf('Performing simple power-based rejection (k=%.2g)...\n', power_k); end
    nStarts = numel(starts);
    rms_vals = zeros(1, nStarts);
    for s = 1:nStarts
        idx = starts(s):(starts(s) + sample_length - 1);
        seg_filtered = filtered(:, idx);
        rms_vals(s) = median(sqrt(mean(seg_filtered .^ 2, 2)));
    end
    med_rms = median(rms_vals);
    mad_rms = mad(rms_vals, 1);
    thr = med_rms + power_k * mad_rms;
    keep_mask = rms_vals <= thr;
    if sum(keep_mask) < 1
        warning('All segments rejected by power threshold. Keeping all segments instead.');
        keep_mask = true(size(keep_mask));
    end
    starts = starts(keep_mask);
    if numel(starts) > n_samples
        if strcmp(selection, 'random')
            starts = starts(randperm(numel(starts), n_samples));
        else
            starts = starts(1:n_samples);
        end
    end
end

n_used = numel(starts);
if verbose
    fprintf('Using %d segments (requested %d). Each segment length = %d samples. Mean overlap = %.3f\n', ...
        n_used, n_samples, sample_length, mean_overlap);
end

%% Compute PLI/wPLI per segment
PLI_vals = nan(n_used, nCh, nCh);
wPLI_vals = nan(n_used, nCh, nCh);
phaseLag_vals = nan(n_used, nCh, nCh);

if verbose, fprintf('Computing PLI/wPLI over segments...\n'); end
for s = 1:n_used
    if verbose && mod(s, 50) == 1
        fprintf('  segment %d / %d\n', s, n_used);
    end
    idx = starts(s):(starts(s) + sample_length - 1);
    seg_analytic = analytic(:, idx);
    for i = 1:nCh
        xi = seg_analytic(i, :);
        for j = i:nCh
            yj = seg_analytic(j, :);
            cs = xi .* conj(yj);
            phlag = angle(mean(cs));
            im_part = imag(cs);
            pli = abs(mean(sign(im_part)));
            numer = abs(mean(im_part));
            denom = mean(abs(im_part));
            wpli = numer / (denom + eps);
            PLI_vals(s, i, j)      = pli;
            PLI_vals(s, j, i)      = pli;
            wPLI_vals(s, i, j)     = wpli;
            wPLI_vals(s, j, i)     = wpli;
            phaseLag_vals(s, i, j) = phlag;
            phaseLag_vals(s, j, i) = -phlag;
        end
    end
end

%% Aggregate results
PLI_matrix  = squeeze(mean(PLI_vals, 1));
wPLI_matrix = squeeze(mean(wPLI_vals, 1));
mean_phaseLag_matrix = squeeze(angle(mean(exp(1i * phaseLag_vals), 1)));

phase_deg = mod(rad2deg(phaseLag_vals) + 360, 360);
phaseLag_matrix_pctls = nan(numel(pctls), nCh, nCh);
for i = 1:nCh
    for j = 1:nCh
        phaseLag_matrix_pctls(:, i, j) = prctile(phase_deg(:, i, j), pctls);
    end
end
phaseLag_matrix_pctls = deg2rad(phaseLag_matrix_pctls);

PLI_matrix_pctls  = squeeze(prctile(PLI_vals, pctls, 1));
wPLI_matrix_pctls = squeeze(prctile(wPLI_vals, pctls, 1));

pli_struct.PLI_matrix            = PLI_matrix;
pli_struct.wPLI_matrix           = wPLI_matrix;
pli_struct.PLI_matrix_pctls      = PLI_matrix_pctls;
pli_struct.wPLI_matrix_pctls     = wPLI_matrix_pctls;
pli_struct.pctls                 = pctls;
pli_struct.sample_length         = sample_length;
pli_struct.n_samples_requested   = n_samples;
pli_struct.n_used                = n_used;
pli_struct.used_starts           = starts;
pli_struct.Hd                    = Hd;
pli_struct.mean_overlap          = mean_overlap;
pli_struct.mean_phaseLag_matrix  = mean_phaseLag_matrix;
pli_struct.phaseLag_matrix_pctls = phaseLag_matrix_pctls;
try
    pli_struct.freq_range = [Hd.CutoffFrequency1 Hd.CutoffFrequency2];
catch
    pli_struct.freq_range = [];
end

if verbose, fprintf('Done.\n'); end
end
