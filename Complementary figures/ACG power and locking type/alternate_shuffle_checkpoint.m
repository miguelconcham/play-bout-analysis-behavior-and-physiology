function out_file = alternate_shuffle_checkpoint(checkpoint, n_shuffles)
%ALTERNATE_SHUFFLE_CHECKPOINT  Second path when the checkpoint cannot be written.
%
% The copy stays in the same folder as checkpoint. Its name ends in _nN
% for this run's permutation count, so a 2000-permutation estimate is not
% written over a 500-permutation file. When checkpoint already ends in
% that count, the copy is _nN_v2 and the open file is left alone.

if nargin < 2 || isempty(n_shuffles) || ~isscalar(n_shuffles) || n_shuffles < 1
    error('alternate_shuffle_checkpoint:count', ...
        'n_shuffles must be the permutation count for this run.');
end
out_dir = fileparts(checkpoint);
[~, base] = fileparts(checkpoint);
count_tag = sprintf('_n%d', round(n_shuffles));
had_count = ~isempty(regexp(base, '_n\d+$', 'once'));
new_base = regexprep(base, '_n\d+$', count_tag);
if ~had_count
    new_base = [base count_tag];
end
out_file = fullfile(out_dir, [new_base '.mat']);
if strcmpi(out_file, checkpoint)
    out_file = fullfile(out_dir, [new_base '_v2.mat']);
end
end
