function rows = oscillatory_member_rows(member, score)
%OSCILLATORY_MEMBER_ROWS  Oscillatory neurons, strongest score first.
%
% member is true for each neuron in the full population. score is the
% modulation score in that same order. Rows with a missing score stay in
% the group and are placed after the finite scores.

member = logical(member(:));
score = double(score(:));
if numel(member) ~= numel(score)
    error('oscillatory_member_rows:size', ...
        'The membership flags and scores must have one entry per neuron.');
end
rows = find(member);
finite = isfinite(score(rows));
finite_rows = rows(finite);
[~, order] = sort(score(finite_rows), 'descend');
rows = [finite_rows(order); rows(~finite)];
end
