function pli_struct = COMPUTE_PLI_WPLI_SEGMENTS(...
    PAG_LFP, Hd, n_samples, sample_length, pctls, varargin)
% Compute PLI and wPLI across channel pairs using multiple segments.
%
% INPUTS:
%   PAG_LFP       : (nCh × nTime) LFP data (double recommended)
%   Hd            : FIR filter object (designed outside)
%   n_samples     : integer number of segments to extract
%   sample_length : integer length of each segment in samples
%   pctls         : array of percentiles (values between 0 and 100)
%
% OPTIONAL NAME-VALUE:
%   'Overlap'     : fraction overlap to force [0..0.99]. If empty (default), use auto-overlap when needed.
%   'Selection'   : 'random' (default) or 'sequential'
%   'PowerThresh' : k for artifact rejection (median + k*MAD). [] disables. Default = 5.
%   'Verbose'     : true/false (default true)
%
% OUTPUT (fields in pli_struct):
%   PLI_matrix, wPLI_matrix                : (nCh × nCh) mean across segments
%   PLI_matrix_pctls, wPLI_matrix_pctls    : (numel(pctls) × nCh × nCh)
%   used_starts, n_used, n_samples_requested
%   mean_overlap                            : mean overlap ratio actually used (0..1)
%   pctls, sample_length, Hd, freq_range
%
% Note: memory ~ n_used * nCh^2 doubles; consider block-processing for very large n_samples & nCh.

% ---------- parse inputs ----------
p = inputParser;
addParameter(p,'Overlap',[], @(x) isempty(x) || (isnumeric(x) && x>=0 && x<1));
addParameter(p,'Selection','random',@(s) any(strcmpi(s,{'random','sequential'})));
addParameter(p,'PowerThresh',5, @(x) isempty(x) || (isnumeric(x) && x>=0));
addParameter(p,'Verbose',true,@islogical);
parse(p,varargin{:});
userOverlap = p.Results.Overlap;
selection = lower(p.Results.Selection);
power_k = p.Results.PowerThresh;
verbose = p.Results.Verbose;

[nCh, nTime] = size(PAG_LFP);
if sample_length > nTime
    error('sample_length (%d) > total recording length (%d).', sample_length, nTime);
end

if verbose
    fprintf('Filtering entire recording (zero-phase FIR)...\n');
end

% ---------- filter & analytic ----------
filtered = filtfilt(Hd.Coefficients, 1, double(PAG_LFP') )'; % nCh x nTime
analytic = hilbert(filtered')';   % nCh x nTime (complex analytic signal)

% ---------- compute starts (auto-overlap if needed, unless userOverlap provided) ----------
n_max_nonoverlap = floor(nTime / sample_length);

if ~isempty(userOverlap)
    % user-specified overlap -> compute step and candidate starts, then select
    ov = userOverlap;
    step = max(1, round(sample_length * (1 - ov)));
    starts_possible = 1:step:(nTime - sample_length + 1);
    if isempty(starts_possible)
        error('No possible windows with the specified Overlap and sample_length.');
    end
    n_possible = numel(starts_possible);
    if strcmp(selection,'random')
        if n_samples <= n_possible
            rand_idx = randperm(n_possible, n_samples);
            starts = starts_possible(rand_idx);
        else
            % take all then sample extra with replacement
            perm = randperm(n_possible);
            starts = starts_possible(perm);
            extra = n_samples - n_possible;
            starts = [starts, starts_possible(randi(n_possible,1,extra))];
            starts = starts(randperm(numel(starts)));
        end
    else % sequential
        if n_samples <= n_possible
            starts = starts_possible(1:n_samples);
        else
            % tile
            reps = ceil(n_samples / n_possible);
            starts = repmat(starts_possible, 1, reps);
            starts = starts(1:n_samples);
        end
    end
    mean_overlap = 1 - (step / sample_length);
else
    % automatic: if requested fits non-overlapping -> no overlap; else compute step to fit
    if n_samples <= n_max_nonoverlap
        step = sample_length;
        starts = 1:step:(nTime - sample_length + 1);
        if numel(starts) > n_samples
            starts = starts(1:n_samples);
        end
        mean_overlap = 0;
        % apply selection policy if we generated more than needed
        if numel(starts) > n_samples
            if strcmp(selection,'random')
                starts = starts(randperm(numel(starts), n_samples));
            else
                starts = starts(1:n_samples);
            end
        end
    else
        % need overlap to fit n_samples: compute step so that first..last segment span recording
        step = floor((nTime - sample_length) / (n_samples - 1));
        if step < 1
            step = 1;
        end
        starts = 1 : step : (1 + step*(n_samples-1));
        % ensure exactly n_samples (due to floor issues)
        if numel(starts) >= n_samples
            starts = starts(1:n_samples);
        else
            % pad by repeating last valid start (shouldn't happen often)
            while numel(starts) < n_samples
                starts(end+1) = starts(end); %#ok<AGROW>
            end
        end
        mean_overlap = 1 - (step / sample_length);
    end
end

% ---------- optional artifact rejection (power-based) ----------
if ~isempty(power_k) && power_k > 0
    if verbose, fprintf('Performing simple power-based rejection (k=%.2g)...\n', power_k); end
    nStarts = numel(starts);
    rms_vals = zeros(1,nStarts);
    for s = 1:nStarts
        idx = starts(s) : (starts(s) + sample_length - 1);
        seg_filtered = filtered(:, idx);
        % robust summary across channels
        rms_vals(s) = median(sqrt(mean(seg_filtered.^2,2)));
    end
    med_rms = median(rms_vals);
    mad_rms = mad(rms_vals,1);
    thr = med_rms + power_k * mad_rms;
    keep_mask = rms_vals <= thr;
    if sum(keep_mask) < 1
        warning('All segments rejected by power threshold. Keeping all segments instead.');
        keep_mask = true(size(keep_mask));
    end
    starts = starts(keep_mask);
    % If after rejection we have more starts than requested, trim:
    if numel(starts) > n_samples
        if strcmp(selection,'random')
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

% ---------- allocate arrays ----------
PLI_vals = nan(n_used, nCh, nCh);
wPLI_vals = nan(n_used, nCh, nCh);
phaseLag_vals = nan(n_used, nCh, nCh);

% ---------- compute PLI/wPLI per segment ----------
if verbose, fprintf('Computing PLI/wPLI over segments...\n'); end
for s = 1:n_used
    if verbose && mod(s,50)==1
        fprintf('  segment %d / %d\n', s, n_used);
    end
    idx = starts(s) : (starts(s) + sample_length - 1);
    seg_analytic = analytic(:, idx); % nCh x sample_length
    for i = 1:nCh
        xi = seg_analytic(i,:);
        for j = i:nCh
            yj = seg_analytic(j,:);
            cs = xi .* conj(yj);
            phlag = angle(mean(cs));
            im_part = imag(cs);
            % PLI
            pli = abs(mean(sign(im_part)));
            % wPLI
            numer = abs(mean(im_part));
            denom = mean(abs(im_part));
            wpli = numer / (denom + eps);
            PLI_vals(s,i,j)         = pli;
            PLI_vals(s,j,i)         = pli;
            wPLI_vals(s,i,j)        = wpli;
            wPLI_vals(s,j,i)        = wpli;
            phaseLag_vals(s,i,j)    = phlag;
            phaseLag_vals(s,j,i)    = -phlag; % antisymmetric
        end
    end
end

% ---------- aggregate ----------
PLI_matrix  = squeeze(mean(PLI_vals, 1));   % nCh x nCh
wPLI_matrix = squeeze(mean(wPLI_vals,1));
mean_phaseLag_matrix = squeeze(angle(mean(exp(1i*phaseLag_vals),1)));

phase_deg = mod(rad2deg(phaseLag_vals) + 360, 360); % n_used x nCh x nCh

phaseLag_matrix_pctls = nan(numel(pctls), nCh, nCh);
for i = 1:nCh
    for j = 1:nCh
        phaseLag_matrix_pctls(:,i,j) = prctile(phase_deg(:,i,j), pctls);
    end
end
% Convert back to radians
phaseLag_matrix_pctls = deg2rad(phaseLag_matrix_pctls);

% percentiles -> [numel(pctls) x nCh x nCh]
PLI_matrix_pctls  = squeeze(prctile(PLI_vals,  pctls, 1));
wPLI_matrix_pctls = squeeze(prctile(wPLI_vals, pctls, 1));

% ---------- package ----------
pli_struct.PLI_matrix        = PLI_matrix;
pli_struct.wPLI_matrix       = wPLI_matrix;
pli_struct.PLI_matrix_pctls  = PLI_matrix_pctls;
pli_struct.wPLI_matrix_pctls = wPLI_matrix_pctls;
pli_struct.pctls             = pctls;
pli_struct.sample_length     = sample_length;
pli_struct.n_samples_requested = n_samples;
pli_struct.n_used            = n_used;
pli_struct.used_starts       = starts;
pli_struct.Hd                = Hd;
pli_struct.mean_overlap      = mean_overlap;
pli_struct.mean_phaseLag_matrix      = mean_phaseLag_matrix;      % nCh x nCh
pli_struct.phaseLag_matrix_pctls     = phaseLag_matrix_pctls;     % numel(pctls) x nCh x nCh
try
    pli_struct.freq_range = [Hd.CutoffFrequency1 Hd.CutoffFrequency2];
catch
    pli_struct.freq_range = [];
end

if verbose, fprintf('Done.\n'); end

end
