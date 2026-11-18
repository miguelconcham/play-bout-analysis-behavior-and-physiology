function print_distance_coverage(report)
%PRINT_DISTANCE_COVERAGE  Command-window count of sessions at each distance.
%
% The bin table is the distance analogue of the high-speed session count:
% how many electrode-1 sessions have play, non-play, or both at that
% within-session distance z.

if isempty(report.animals)
    fprintf('No electrode-1 sessions.\n');
    return
end

options = report.options;
n_sess = report.n_sessions;
wind_index = find(isfinite(report.wind_s), 1, 'first');
hop_index = find(isfinite(report.hop_s), 1, 'first');
fprintf('\nInter-animal distance, electrode 1.\n');
if ~isempty(hop_index) && ~isempty(wind_index)
    fprintf('Frame hop %.4f s, window %.3f s.\n', ...
        report.hop_s(hop_index), report.wind_s(wind_index));
    fprintf(['Durations are spans of overlapping spectrogram frames. ' ...
        '%.1f s is about %.1f windows.\n'], ...
        options.min_nonplay_s, options.min_nonplay_s / report.wind_s(wind_index));
end
fprintf(['A session has enough far distance when non-play >= %.1f s and play >= %.1f s at z >= %.1f.\n' ...
    'A session has enough close distance when those same minima hold at z <= %.1f.\n'], ...
    options.min_nonplay_s, options.min_play_s, options.high_z, -options.high_z);

fprintf(['\nSessions at each %.2f-wide distance bin.\n' ...
    '  play: sessions with >= %.1f s of play in the bin.\n' ...
    '  non-play: sessions with >= %.1f s of non-play in the bin.\n' ...
    '  both: sessions that can be compared at that distance.\n'], ...
    options.curve_bin_width, options.min_bin_play_s, options.min_bin_nonplay_s);

n_bins = numel(report.curve_x);
printed = false;
for b = 1:n_bins
    if report.n_play(b) == 0 && report.n_nonplay(b) == 0 && report.n_both(b) == 0
        continue
    end
    printed = true;
    if b == n_bins
        right_bracket = ']';
    else
        right_bracket = ')';
    end
    fprintf('  [%.2f, %.2f%c  play %d, non-play %d, both %d\n', ...
        report.curve_edges(b), report.curve_edges(b + 1), right_bracket, ...
        report.n_play(b), report.n_nonplay(b), report.n_both(b));
end
if ~printed
    fprintf('  no bin meets a duration minimum\n');
end

fprintf('\nSessions with enough non-play and play at distance z >= %.1f: %d / %d\n', ...
    options.high_z, nnz(report.enough_far), n_sess);
print_difference_test(report.diff_far(report.enough_far), report.animals(report.enough_far), ...
    sprintf('distance z >= %.1f, sessions with enough overlap', options.high_z));

fprintf('Sessions with enough non-play and play at distance z <= %.1f: %d / %d\n', ...
    -options.high_z, nnz(report.enough_close), n_sess);
print_difference_test(report.diff_close(report.enough_close), report.animals(report.enough_close), ...
    sprintf('distance z <= %.1f, sessions with enough overlap', -options.high_z));

far_bin_is_far = isfinite(report.far_bin_lo) & report.far_bin_lo >= options.high_z;
fprintf('Sessions whose farthest supported bin starts at z >= %.1f: %d / %d\n', ...
    options.high_z, nnz(far_bin_is_far), n_sess);
print_difference_test(report.diff_far_bin(far_bin_is_far), report.animals(far_bin_is_far), ...
    sprintf('farthest supported distance bin at z >= %.1f', options.high_z));

close_bin_is_close = isfinite(report.close_bin_hi) & report.close_bin_hi <= -options.high_z;
fprintf('Sessions whose closest supported bin stays at z <= %.1f: %d / %d\n', ...
    -options.high_z, nnz(close_bin_is_close), n_sess);
print_difference_test(report.diff_close_bin(close_bin_is_close), report.animals(close_bin_is_close), ...
    sprintf('closest supported distance bin at z <= %.1f', -options.high_z));

fprintf('\nPer session, longest far non-play span first.\n');
[~, order] = sort(report.nonplay_far_s, 'descend');
for k = 1:n_sess
    s = order(k);
    fprintf('%s P%d\n', report.animals{s}, report.partners(s));
    fprintf('  z >= %.1f: non-play %.2f s, play %.2f s, session max non-play z %.2f, session max play z %.2f\n', ...
        options.high_z, report.nonplay_far_s(s), report.play_far_s(s), ...
        report.max_nonplay(s), report.max_play(s));
    fprintf('  difference at z >= %.1f: %s\n', options.high_z, ...
        diff_text(report.diff_far(s), report.diff_far_z(s), report.enough_far(s), ...
        'passes minimum', 'below the minimum duration'));
    if isfinite(report.far_bin_lo(s))
        fprintf('  farthest supported bin [%.2f, %.2f]: non-play %.2f s, play %.2f s, %s\n', ...
            report.far_bin_lo(s), report.far_bin_hi(s), ...
            report.far_bin_nonplay_s(s), report.far_bin_play_s(s), ...
            diff_text(report.diff_far_bin(s), report.diff_far_bin_z(s), ...
            report.far_bin_lo(s) >= options.high_z, ...
            'bin starts at the far cut', 'bin starts below the far cut'));
    else
        fprintf('  no %.2f-wide bin has >= %.1f s play and >= %.1f s non-play\n', ...
            options.ceiling_bin_width, options.min_ceiling_play_s, ...
            options.min_ceiling_nonplay_s);
    end
    fprintf('  z <= %.1f: non-play %.2f s, play %.2f s, session min non-play z %.2f, session min play z %.2f\n', ...
        -options.high_z, report.nonplay_close_s(s), report.play_close_s(s), ...
        report.min_nonplay(s), report.min_play(s));
    fprintf('  difference at z <= %.1f: %s\n', -options.high_z, ...
        diff_text(report.diff_close(s), report.diff_close_z(s), report.enough_close(s), ...
        'passes minimum', 'below the minimum duration'));
    if isfinite(report.close_bin_lo(s))
        fprintf('  closest supported bin [%.2f, %.2f]: non-play %.2f s, play %.2f s, %s\n', ...
            report.close_bin_lo(s), report.close_bin_hi(s), ...
            report.close_bin_nonplay_s(s), report.close_bin_play_s(s), ...
            diff_text(report.diff_close_bin(s), report.diff_close_bin_z(s), ...
            report.close_bin_hi(s) <= -options.high_z, ...
            'bin stays at the close cut', 'bin extends past the close cut'));
    end
end
end

function text = diff_text(diff_log, diff_z, passes, pass_label, fail_label)
if ~isfinite(diff_log)
    text = 'no overlap';
    return
end
amp_pct = 100 * (10 ^ diff_log - 1);
if passes
    flag = pass_label;
else
    flag = fail_label;
end
text = sprintf('%+.4f log10 (%+.1f%% amplitude, %+.3f session z), %s', ...
    diff_log, amp_pct, diff_z, flag);
end

function print_difference_test(vals, animals, label)
vals = vals(:);
animals = animals(:);
keep = isfinite(vals);
vals = vals(keep);
animals = animals(keep);
fprintf('%s: n = %d\n', label, numel(vals));
if numel(vals) < 2
    return
end
[~, p_val, ci, stats] = ttest(vals);
fprintf('  one-sample t: mean %.4f log10, 95%% CI [%.4f, %.4f], t(%d) = %.2f, p = %.4f\n', ...
    mean(vals), ci(1), ci(2), stats.df, stats.tstat, p_val);
if numel(unique(animals)) < 2
    fprintf('  LME skipped: only one animal\n');
    return
end
if exist('fitlme', 'file') ~= 2
    fprintf('  LME skipped: fitlme is not available\n');
    return
end
lme = fitlme(table(vals, categorical(animals), ...
    'VariableNames', {'Difference', 'Animal'}), ...
    'Difference ~ 1 + (1|Animal)');
coef = fixedEffects(lme);
ci_lme = coefCI(lme);
fprintf('  LME Difference ~ 1 + (1|Animal): mean %.4f, 95%% CI [%.4f, %.4f], p = %.4f\n', ...
    coef, ci_lme(1), ci_lme(2), lme.Coefficients.pValue(1));
end
