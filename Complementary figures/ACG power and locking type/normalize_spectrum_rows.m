function stack = normalize_spectrum_rows(stack)
%NORMALIZE_SPECTRUM_ROWS  Divide each spectrum by its own total power.
%
% Finite bins of a row are divided by their sum. A row with no finite
% power, or with a zero sum, stays NaN.

stack = double(stack);
if isempty(stack)
    return
end
if isvector(stack)
    stack = stack(:).';
end
for i = 1:size(stack, 1)
    stack(i, :) = divide_by_total(stack(i, :), spectrum_total(stack(i, :)));
end
end


function total = spectrum_total(row)
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


function row = divide_by_total(row, total)
if ~isfinite(total)
    row(:) = nan;
    return
end
finite = isfinite(row);
row(finite) = row(finite) / total;
end
