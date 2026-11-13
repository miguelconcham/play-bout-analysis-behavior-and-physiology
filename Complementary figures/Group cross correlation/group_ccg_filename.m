function name = group_ccg_filename(animal_code)
%GROUP_CCG_FILENAME  Per-session file written by Estimate_group_cross_correlation.

name = [ccg_plain_text(animal_code), '_group_ccg.mat'];
end
