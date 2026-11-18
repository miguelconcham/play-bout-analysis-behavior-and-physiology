function [group_id, group_index] = assign_ccg_group(lock_name, oscillatory)
%ASSIGN_CCG_GROUP  Map a lock label and an oscillation flag onto one group.
%
% lock_name accepts 'peak', 'trough', and 'unlocked'. 'nonentrained' and
% 'non-entrained' are unlocked. oscillatory is true, false, 1, or 0.
% group_id is '' when the pair is not one of the six groups.

names = group_ccg_names();
group_id = '';
group_index = 0;
lock_name = normalize_lock(lock_name);
if isempty(lock_name) || ~(islogical(oscillatory) || isnumeric(oscillatory)) ...
        || ~isscalar(oscillatory) || ~isfinite(double(oscillatory))
    return
end
is_oscillatory = logical(oscillatory);
hit = find(strcmp(names.lock, lock_name) & names.oscillatory == is_oscillatory, 1);
if isempty(hit)
    return
end
group_id = names.group_id{hit};
group_index = hit;
end


function lock_name = normalize_lock(lock_name)
lock_name = lower(ccg_plain_text(lock_name));
if any(strcmp(lock_name, {'nonentrained', 'non-entrained', 'unlocked'}))
    lock_name = 'unlocked';
elseif ~any(strcmp(lock_name, {'peak', 'trough'}))
    lock_name = '';
end
end
