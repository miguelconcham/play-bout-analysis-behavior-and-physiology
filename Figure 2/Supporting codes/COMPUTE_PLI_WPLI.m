function pli_struct = COMPUTE_PLI_WPLI(PAG_LFP, Hd, n_samples, sample_length, pctls)
% Compute PLI and wPLI across channel pairs using multiple segments
%
% INPUTS:
%   PAG_LFP       : (nCh × nTime) LFP data
%   Hd            : FIR bandpass filter (designed outside with designfilt)
%   n_samples     : number of segments to extract
%   sample_length : length of each segment (in samples)
%   pctls         : array of percentiles to compute (0–100)
%
% OUTPUT:
%   pli_struct.PLI_matrix         : (nCh × nCh), mean PLI across segments
%   pli_struct.wPLI_matrix        : (nCh × nCh), mean wPLI across segments
%   pli_struct.PLI_matrix_pctls   : (numel(pctls) × nCh × nCh)
%   pli_struct.wPLI_matrix_pctls  : (numel(pctls) × nCh × nCh)
%   pli_struct.freq_range         : filter band
%   pli_struct.Hd                 : filter object

% --------------------------------------------------------
% Step 1: filter data
disp('Filtering data...')
filtered = filtfilt(Hd.Coefficients, 1, double(PAG_LFP'))'; % (ch × time)

% Step 2: analytic signals
disp('Hilbert transform...')
analytic = hilbert(filtered')';   % (ch × time)
nCh = size(analytic,1);
nTime = size(analytic,2);

%Step 3: segment selection (auto overlap if needed) --------------------
n_max_nonoverlap = floor(nTime / sample_length);

if n_samples <= n_max_nonoverlap
    % Non-overlap case
    step = sample_length;
    starts = 1:step:(nTime - sample_length + 1);
    if numel(starts) > n_samples
        starts = starts(1:n_samples);
    end
    mean_overlap = 0;
else
    % Need overlap: compute step so that n_samples segments fit
    step = floor((nTime - sample_length) / (n_samples - 1));
    if step < 1
        step = 1; % extreme case
    end
    starts = 1:step:(1 + step*(n_samples-1));
    starts = starts(1:n_samples); % enforce exact count
    mean_overlap = 1 - (step / sample_length);
end

% Take first n_samples contiguous segments (non-overlapping)
segments = zeros(nCh, sample_length, n_samples, 'like', analytic);
for s = 1:n_samples
    idx = starts(s):(starts(s)+sample_length-1);
    segments(:,:,s) = analytic_full(:, idx);
end

% Step 4: compute PLI/wPLI per segment
PLI_vals  = nan(n_samples, nCh, nCh);
wPLI_vals = nan(n_samples, nCh, nCh);

disp('Computing PLI/wPLI per segment...')
for s = 1:n_samples
    seg = segments(:,:,s); % (nCh × sample_length)
    for i = 1:nCh
        for j = i:nCh
            x = seg(i,:);
            y = seg(j,:);
            
            cs = x .* conj(y);
            im_part = imag(cs);
            
            % PLI
            pli = abs(mean(sign(im_part)));
            % wPLI
            num = abs(mean(im_part));
            den = mean(abs(im_part));
            wpli = num / (den + eps);
            
            PLI_vals(s,i,j)  = pli;
            PLI_vals(s,j,i)  = pli;
            wPLI_vals(s,i,j) = wpli;
            wPLI_vals(s,j,i) = wpli;
        end
    end
end

% Step 5: aggregate across samples
PLI_matrix  = squeeze(mean(PLI_vals,1));   % mean across segments
wPLI_matrix = squeeze(mean(wPLI_vals,1));

PLI_matrix_pctls  = squeeze(prctile(PLI_vals,  pctls, 1));
wPLI_matrix_pctls = squeeze(prctile(wPLI_vals, pctls, 1));

% Step 6: package results
pli_struct.PLI_matrix        = PLI_matrix;
pli_struct.wPLI_matrix       = wPLI_matrix;
pli_struct.PLI_matrix_pctls  = PLI_matrix_pctls;
pli_struct.wPLI_matrix_pctls = wPLI_matrix_pctls;
pli_struct.pctls             = pctls;
pli_struct.freq_range        = [Hd.CutoffFrequency1 Hd.CutoffFrequency2];
pli_struct.Hd                = Hd;
pli_struct.mean_overlap      = mean_overlap;


disp('PLI/wPLI distribution estimation done.')

end
