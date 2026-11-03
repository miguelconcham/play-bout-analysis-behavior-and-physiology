function partner_index = assign_partner_index(bout_start, window_start, window_end)
%ASSIGN_PARTNER_INDEX  Partner block that contains each play-bout onset.
%
%   Windows are partner 1, 2, ... in the order given. A bout belongs to the
%   first window with window_start <= onset < window_end. Onsets outside
%   every window are NaN.

bout_start = double(bout_start(:));
window_start = double(window_start(:));
window_end = double(window_end(:));
partner_index = nan(numel(bout_start), 1);
if isempty(bout_start) || isempty(window_start)
    return
end
if numel(window_start) ~= numel(window_end)
    error('assign_partner_index:Windows', ...
        'Partner windows need one start and one end.');
end

for p = 1:numel(window_start)
    in_window = bout_start >= window_start(p) & bout_start < window_end(p);
    partner_index(in_window & isnan(partner_index)) = p;
end
end
