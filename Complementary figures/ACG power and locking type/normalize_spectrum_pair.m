function [high_out, low_out] = normalize_spectrum_pair(high_spec, low_spec, reference)
%NORMALIZE_SPECTRUM_PAIR  Scale high and low spectra by one reference.
%
% Each neuron uses the total power of its reference spectrum, the all-runs
% spectrum. Both the high-power and low-power rows are divided by that
% same total. When the reference row has no power, the total of the
% finite high and low bins is used instead, still one number for both.

high_out = double(high_spec);
low_out = double(low_spec);
reference = double(reference);
if isempty(high_out) && isempty(low_out)
    return
end
if isvector(high_out)
    high_out = high_out(:).';
end
if isvector(low_out)
    low_out = low_out(:).';
end
if isvector(reference)
    reference = reference(:).';
end
n_rows = max(size(high_out, 1), size(low_out, 1));
n_freq = max(size(high_out, 2), size(low_out, 2));
if size(high_out, 1) ~= n_rows
    high_out = nan(n_rows, n_freq);
end
if size(low_out, 1) ~= n_rows
    low_out = nan(n_rows, n_freq);
end
for i = 1:n_rows
    total = nan;
    if size(reference, 1) >= i && size(reference, 2) > 0
        total = row_total(reference(i, :));
    end
    if ~isfinite(total)
        total = row_total([high_out(i, :), low_out(i, :)]);
    end
    high_out(i, :) = divide_row(high_out(i, :), total);
    low_out(i, :) = divide_row(low_out(i, :), total);
end
end


function total = row_total(row)
finite = isfinite(row);
if ~any(finite)
    total = nan;
    return
end
total = sum(row(finite));
if ~isfinite(total) || total == 0
    total = nan;
end
end


function row = divide_row(row, total)
if ~isfinite(total)
    row(:) = nan;
    return
end
finite = isfinite(row);
row(finite) = row(finite) / total;
end
