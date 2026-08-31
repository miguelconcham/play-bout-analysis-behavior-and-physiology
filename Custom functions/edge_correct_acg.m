function acg = edge_correct_acg(acg, bin_edges, run_length)
%EDGE_CORRECT_ACG  Divide finite-window ACG counts by overlap (T - |lag|).
%
%   acg = edge_correct_acg(counts, bin_edges, run_length)
%
% Pair counts at lag tau can only come from the overlapping part of a run
% of length T, which shrinks as (T - |tau|). Bins with |lag| >= T, or with
% overlap shorter than one bin, are set to NaN so they are not treated as
% true zeros when averaging.

    if isempty(acg)
        return
    end
    bin_edges = bin_edges(:)';
    if numel(bin_edges) < 2
        return
    end
    bin_w = mean(diff(bin_edges));
    centers = bin_edges(1:end-1) + 0.5 * bin_w;
    run_length = run_length(:);
    n_row = size(acg, 1);
    if isscalar(run_length) && n_row > 1
        run_length = repmat(run_length, n_row, 1);
    end
    if numel(run_length) ~= n_row
        error('edge_correct_acg:size', 'run_length must match the number of ACG rows.');
    end
    overlap = run_length - abs(centers);
    overlap(overlap < bin_w | ~isfinite(overlap)) = NaN;
    acg = double(acg) ./ overlap;
end
