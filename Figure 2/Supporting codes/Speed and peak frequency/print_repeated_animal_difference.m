function print_repeated_animal_difference(differences, animals, label)
%PRINT_REPEATED_ANIMAL_DIFFERENCE  Is play minus non-play above 0 in general?
%
% Sessions stay nested in animals. The equal-animal test gives every animal
% one mean. The mixed model keeps every session and adds (1|Animal).

summary = animal_level_means(differences, animals);
fprintf('\n%s\n', label);
fprintf('  Sessions with both play and non-play: %d. Animals: %d.\n', ...
    summary.n_sessions, summary.n_animals);
if summary.n_sessions == 0
    fprintf('  No session has both play and non-play here.\n');
    return
end
for k = 1:summary.n_animals
    fprintf('  %s  sessions %d  mean diff %+.4f\n', ...
        summary.animals{k}, summary.counts(k), summary.means(k));
end
fprintf('  Equal-animal mean %+.4f log10\n', summary.equal_animal_mean);
if summary.n_animals < 2
    fprintf('  A general test needs at least two animals.\n');
    return
end
if exist('ttest', 'file') == 2
    [~, p_val, ci, stats] = ttest(summary.means);
    fprintf(['  Equal-animal t, one mean per animal: 95%% CI [%.4f, %.4f], ' ...
        't(%d) = %.2f, p = %.4f\n'], ci(1), ci(2), stats.df, stats.tstat, p_val);
    print_above_zero(ci, 'equal-animal');
else
    fprintf('  Equal-animal t-test was not available.\n');
end

if exist('fitlme', 'file') ~= 2
    fprintf(['  LME skipped: fitlme is not available. The equal-animal ' ...
        'mean gives each animal one observation.\n']);
    return
end
if summary.n_sessions < 3
    fprintf('  LME skipped: fewer than three sessions.\n');
    return
end
keep = isfinite(differences(:));
session_diff = differences(keep);
session_animals = animals(:);
session_animals = session_animals(keep);
try
    lme = fitlme(table(session_diff, categorical(session_animals), ...
        'VariableNames', {'Difference', 'Animal'}), ...
        'Difference ~ 1 + (1|Animal)');
    coef = fixedEffects(lme);
    ci_lme = coefCI(lme);
    fprintf(['  LME Difference ~ 1 + (1|Animal), every session kept: ' ...
        'mean %.4f, 95%% CI [%.4f, %.4f], p = %.4f\n'], ...
        coef, ci_lme(1), ci_lme(2), lme.Coefficients.pValue(1));
    print_above_zero(ci_lme, 'LME');
catch
    fprintf(['  LME could not separate the animal effect from the ' ...
        'session scatter.\n']);
end
end

function print_above_zero(ci, name)
if ci(1) > 0
    fprintf(['  The %s 95%% CI is above 0, so play power is higher ' ...
        'in general.\n'], name);
elseif ci(2) < 0
    fprintf(['  The %s 95%% CI is below 0, so non-play power is higher ' ...
        'in general.\n'], name);
else
    fprintf('  The %s 95%% CI includes 0.\n', name);
end
end
