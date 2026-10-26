function partner_index = partner_index_in_windows(bout_start, bout_end, window_start, window_end)
%PARTNER_INDEX_IN_WINDOWS  Partner index used by the original play-bout script.
%
%   Session blocks are in the order given, after Tickling has been removed.
%   A bout takes the first block that contains it completely:
%   block start <= bout start and bout end < block end.
%   max() on that test is the original rule, so a bout outside every block
%   stays on partner 1.

bout_start = double(bout_start(:));
bout_end = double(bout_end(:));
window_start = double(window_start(:));
window_end = double(window_end(:));
if numel(bout_start) ~= numel(bout_end)
    error('partner_index_in_windows:Bouts', 'Each bout needs a start and an end.');
end
if numel(window_start) ~= numel(window_end) || isempty(window_start)
    error('partner_index_in_windows:Windows', ...
        'Partner blocks need one start and one end.');
end

inside = bout_start.' >= window_start & bout_end.' < window_end;
[~, partner_index] = max(inside, [], 1);
partner_index = double(partner_index(:));
end
