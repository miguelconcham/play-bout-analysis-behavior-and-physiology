function pli_struct = Compute_pli_wpli(PAG_LFP, Hd, n_samples, sample_length, pctls)
% Compute PLI and wPLI across channel pairs using multiple segments.
%
% OUTPUT:
%   pli_struct.PLI_matrix         : (nCh x nCh), mean PLI across segments
%   pli_struct.wPLI_matrix        : (nCh x nCh), mean wPLI across segments
%   pli_struct.PLI_matrix_pctls   : (numel(pctls) x nCh x nCh)
%   pli_struct.wPLI_matrix_pctls  : (numel(pctls) x nCh x nCh)
%   pli_struct.freq_range         : filter band
%   pli_struct.Hd                 : filter object

%% Filter and compute analytic signals
disp('Filtering data...')
filtered = filtfilt(Hd.Coefficients, 1, double(PAG_LFP')')';
disp('Hilbert transform...')
analytic = hilbert(filtered')';
[nCh, nTime] = size(analytic);

%% Select segment starts
n_max_nonoverlap = floor(nTime / sample_length);
if n_samples <= n_max_nonoverlap
    step = sample_length;
    starts = 1:step:(nTime - sample_length + 1);
    if numel(starts) > n_samples
        starts = starts(1:n_samples);
    end
    mean_overlap = 0;
else
    step = floor((nTime - sample_length) / (n_samples - 1));
    if step < 1
        step = 1;
    end
    starts = 1:step:(1 + step * (n_samples - 1));
    starts = starts(1:n_samples);
    mean_overlap = 1 - (step / sample_length);
end

segments = zeros(nCh, sample_length, n_samples, 'like', analytic);
for s = 1:n_samples
    idx = starts(s):(starts(s) + sample_length - 1);
    segments(:, :, s) = analytic(:, idx);
end

%% Compute PLI/wPLI per segment
PLI_vals  = nan(n_samples, nCh, nCh);
wPLI_vals = nan(n_samples, nCh, nCh);

disp('Computing PLI/wPLI per segment...')
for s = 1:n_samples
    seg = segments(:, :, s);
    for i = 1:nCh
        for j = i:nCh
            x = seg(i, :);
            y = seg(j, :);
            cs = x .* conj(y);
            im_part = imag(cs);
            pli = abs(mean(sign(im_part)));
            wpli = abs(mean(im_part)) / (mean(abs(im_part)) + eps);
            PLI_vals(s, i, j)  = pli;
            PLI_vals(s, j, i)  = pli;
            wPLI_vals(s, i, j) = wpli;
            wPLI_vals(s, j, i) = wpli;
        end
    end
end

%% Aggregate across segments
PLI_matrix  = squeeze(mean(PLI_vals, 1));
wPLI_matrix = squeeze(mean(wPLI_vals, 1));
PLI_matrix_pctls  = squeeze(prctile(PLI_vals, pctls, 1));
wPLI_matrix_pctls = squeeze(prctile(wPLI_vals, pctls, 1));

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
