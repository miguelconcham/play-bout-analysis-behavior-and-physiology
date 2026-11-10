function ok = acg_measure_mask(scores_a, scores_b, measure_name)
%ACG_MEASURE_MASK  Neurons whose paired autocorrelogram property is defined.
%
% The same neurons are kept as in the locking-period paired plots. A
% prominence needs both curves to have been searched. A latency needs the
% peak to be present on both curves. A second-cycle measure also needs
% the second window to have been observed on both curves.

ok = isfinite(scores_a.(measure_name)) & isfinite(scores_b.(measure_name));
switch measure_name
    case 'second_prominence'
        ok = ok & scores_a.second_measurable & scores_b.second_measurable;
    case 'first_prominence'
        ok = ok & scores_a.first_measurable & scores_b.first_measurable;
    case 'first_lag'
        ok = ok & scores_a.first_present & scores_b.first_present;
    case 'second_lag'
        ok = ok & scores_a.second_present & scores_b.second_present ...
            & scores_a.second_measurable & scores_b.second_measurable;
end
ok = ok(:);
end
