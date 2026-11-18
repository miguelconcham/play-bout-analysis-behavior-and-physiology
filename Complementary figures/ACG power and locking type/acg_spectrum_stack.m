function stack = acg_spectrum_stack(acg, centers, delta_lag, f_grid)
%ACG_SPECTRUM_STACK  Power spectrum of each autocorrelogram on one frequency grid.
%
% Each row of acg is one neuron. The spectrum is the folded autocorrelogram
% after 100 ms, from folded_acg_spectrum, interpolated onto f_grid. A row
% with fewer than two finite frequencies stays NaN.

f_grid = double(f_grid(:)).';
acg = double(acg);
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
for row = 1:n_rows
    if numel(freq{row}) < 2 || numel(spec{row}) < 2
        continue
    end
    row_freq = double(freq{row}(:));
    row_spec = double(spec{row}(:));
    n = min(numel(row_freq), numel(row_spec));
    row_freq = row_freq(1:n);
    row_spec = row_spec(1:n);
    keep = isfinite(row_freq) & isfinite(row_spec);
    row_freq = row_freq(keep);
    row_spec = row_spec(keep);
    if numel(row_freq) < 2
        continue
    end
    [row_freq, order] = sort(row_freq);
    row_spec = row_spec(order);
    [row_freq, unique_at] = unique(row_freq, 'stable');
    row_spec = row_spec(unique_at);
    if numel(row_freq) < 2
        continue
    end
    stack(row, :) = interp1(row_freq, row_spec, f_grid, 'linear', nan);
end
end
