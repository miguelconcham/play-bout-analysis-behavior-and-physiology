function physical = z_to_physical(z_values, z_samples, physical_samples)
%Z_TO_PHYSICAL  Express session-z values in the paired physical unit.
%
% z_samples and physical_samples come from the same rows. The saved speed
% and distance z-scores are an affine transform of the pixel series, so the
% line through those pairs converts a bin edge from z to cm or cm/s.

original_size = size(z_values);
z_values = double(z_values(:));
physical = nan(numel(z_values), 1);
if isempty(z_samples) || isempty(physical_samples)
    physical = reshape(physical, original_size);
    return
end
z_samples = double(z_samples(:));
physical_samples = double(physical_samples(:));
n = min(numel(z_samples), numel(physical_samples));
ok = isfinite(z_samples(1:n)) & isfinite(physical_samples(1:n));
z_ok = z_samples(ok);
if numel(z_ok) < 2 || numel(unique(z_ok)) < 2
    physical = reshape(physical, original_size);
    return
end
coef = polyfit(z_ok, physical_samples(ok), 1);
finite_value = isfinite(z_values);
physical(finite_value) = coef(1) * z_values(finite_value) + coef(2);
physical = reshape(physical, original_size);
end
