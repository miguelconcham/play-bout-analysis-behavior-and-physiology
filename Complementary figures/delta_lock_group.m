function group_name = delta_lock_group(ppc_p, preferred_angle, alpha, non_entrained_p)
%DELTA_LOCK_GROUP  Trough, peak, or non-entrained from the delta table.
%
% group_name is 'trough', 'peak', 'nonentrained', or '' when the cell is
% between the lock threshold and the non-entrained threshold, or when the
% preferred angle is missing. Trough is a locked cell whose preferred
% angle lies outside [-pi/2, pi/2], the same split as Figure 3.

group_name = '';
if ~(isscalar(ppc_p) && isscalar(preferred_angle) && isfinite(ppc_p) && isfinite(preferred_angle))
    return
end
if ppc_p < alpha
    if preferred_angle > pi / 2 || preferred_angle < -pi / 2
        group_name = 'trough';
    else
        group_name = 'peak';
    end
elseif ppc_p >= non_entrained_p
    group_name = 'nonentrained';
end
end
