function seconds = mixture_remaining_seconds( ...
        elapsed_seconds, done_units, remaining_units)
%MIXTURE_REMAINING_SECONDS  Project the unfinished mixture search.
%
% done_units and remaining_units use the same work scale: centers times
% restarts times the component count of each fit.

if ~isfinite(elapsed_seconds) || ~isfinite(done_units) || done_units <= 0
    seconds = NaN;
    return
end
seconds = elapsed_seconds / done_units * remaining_units;
end
