function [stack, f_range] = normalize_spectrum_like_score(acg, centers, delta_lag, f_grid)
%NORMALIZE_SPECTRUM_LIKE_SCORE  Power divided by the full-spectrum mean.
%
% Each row is one autocorrelogram. Its power is the folded spectrum after
% 100 ms. Every positive Fourier frequency, through the bin below the
% Nyquist frequency, enters the mean. The displayed bins are that power
% divided by the mean, interpolated onto f_grid. This is the power-spectrum
% form of OscillationScore, which divides its 2-5 Hz peak by the mean of
% the whole transform. f_range is the frequencies that entered the mean.

f_grid = double(f_grid(:)).';
acg = double(acg);
f_range = [nan nan];
if isempty(acg)
    stack = nan(0, numel(f_grid));
    return
end
if isvector(acg)
    acg = acg(:).';
end
n_rows = size(acg, 1);
stack = nan(n_rows, numel(f_grid));
if n_rows < 1 || isempty(f_grid)
    return
end
[freq, spec] = folded_acg_spectrum(acg, centers, delta_lag);
used = zeros(0, 1);
for row = 1:n_rows
    [row_freq, row_spec] = sorted_spectrum(freq{row}, spec{row});
    if numel(row_freq) < 2
        continue
    end
    reference = mean(row_spec);
    if ~isfinite(reference) || reference <= 0
        continue
    end
    stack(row, :) = interp1(row_freq, row_spec, f_grid, 'linear', nan) / reference;
    used = [used; row_freq(:)]; %#ok<AGROW>
end
finite = isfinite(used) & used > 0;
if any(finite)
    f_range = [min(used(finite)), max(used(finite))];
end
end


function [row_freq, row_spec] = sorted_spectrum(freq, spec)
row_freq = zeros(0, 1);
row_spec = zeros(0, 1);
freq = double(freq(:));
spec = double(spec(:));
n = min(numel(freq), numel(spec));
if n < 2
    return
end
freq = freq(1:n);
spec = spec(1:n);
keep = isfinite(freq) & isfinite(spec) & freq > 0;
freq = freq(keep);
spec = spec(keep);
if numel(freq) < 2
    return
end
[freq, order] = sort(freq);
spec = spec(order);
[row_freq, unique_at] = unique(freq, 'stable');
row_spec = spec(unique_at);
end
