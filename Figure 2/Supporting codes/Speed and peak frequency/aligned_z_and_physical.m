function [z_kept, physical_kept] = aligned_z_and_physical(model_data, z_name, physical_name)
%ALIGNED_Z_AND_PHYSICAL  Physical column on the rows used for the z analysis.
%
% Rows with a missing z-score or a missing power are dropped, matching the
% coverage code. A missing physical column returns two empty vectors.

z_kept = [];
physical_kept = [];
if ~model_has_column(model_data, z_name) || ~model_has_column(model_data, physical_name) ...
        || ~model_has_column(model_data, 'Power')
    return
end
z = double(model_data.(z_name)(:));
power = double(model_data.Power(:));
physical = double(model_data.(physical_name)(:));
n = min([numel(z), numel(power), numel(physical)]);
keep = isfinite(z(1:n)) & isfinite(power(1:n));
z_kept = z(1:n);
z_kept = z_kept(keep);
physical_kept = physical(1:n);
physical_kept = physical_kept(keep);
end

function present = model_has_column(model_data, name)
if isstruct(model_data)
    present = isfield(model_data, name);
    return
end
present = any(strcmp(name, model_data.Properties.VariableNames));
end
