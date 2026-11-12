function Y = normalize_acg_rows(X)
%NORMALIZE_ACG_ROWS  Scale each autocorrelogram so its finite bins sum to 1.
%
% A neuron with more spike pairs then has the same weight as a neuron with
% fewer pairs. Bins that are NaN stay NaN and are left out of the sum.
% A row with no finite positive mass becomes NaN.

if isempty(X)
    Y = X;
    return
end
X = double(X);
filled = X;
filled(isnan(filled)) = 0;
mass = sum(filled, 2);
mass(~isfinite(mass) | mass <= 0) = NaN;
Y = X ./ mass;
Y(~isfinite(Y)) = NaN;
end
