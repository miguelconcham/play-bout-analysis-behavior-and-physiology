function results = grand_average_erp_lme(session_erp, animal_id, time_s, alpha)
%GRAND_AVERAGE_ERP_LME  Grand-mean ERP from one mixed model per time point.
%
% Each row of session_erp is one recording's mean waveform. The model is
%
%   Voltage ~ 1 + (1|Animal)
%
% The intercept is the grand mean. Its 95% interval is results.ci. Recordings
% from the same animal share a random intercept, so repeated days are not
% treated as independent animals. At least two animals and two recordings
% are required.

if nargin < 4 || isempty(alpha)
    alpha = 0.05;
end
if ~isscalar(alpha) || ~isfinite(alpha) || alpha <= 0 || alpha >= 1
    error('grand_average_erp_lme:alpha', 'alpha must lie in (0, 1).');
end
if exist('fitlme', 'file') ~= 2
    error('grand_average_erp_lme:fitlme', ...
        'fitlme was not found. The Statistics and Machine Learning Toolbox is required.');
end

time_s = double(time_s(:));
session_erp = double(session_erp);
animal_id = animal_labels(animal_id);
if isvector(session_erp) && numel(session_erp) == numel(time_s)
    session_erp = session_erp(:).';
end
if size(session_erp, 2) ~= numel(time_s)
    error('grand_average_erp_lme:time', ...
        'session_erp must have one column per time point.');
end
if numel(animal_id) ~= size(session_erp, 1)
    error('grand_average_erp_lme:animal', ...
        'animal_id must have one entry per session ERP.');
end

finite_row = all(isfinite(session_erp), 2);
session_erp = session_erp(finite_row, :);
animal_id = animal_id(finite_row);
if size(session_erp, 1) < 2
    error('grand_average_erp_lme:sessions', ...
        'At least two session ERPs are required.');
end
if numel(unique(animal_id)) < 2
    error('grand_average_erp_lme:animals', ...
        'At least two animals are required for Voltage ~ 1 + (1|Animal).');
end

n_time = numel(time_s);
estimate = nan(n_time, 1);
se = nan(n_time, 1);
ci = nan(n_time, 2);
pvals = nan(n_time, 1);
n_failed = 0;
first_error = '';

for time_n = 1:n_time
    if mod(time_n, 1000) == 0
        fprintf('Mixed model %d/%d\n', time_n, n_time);
    end
    y = session_erp(:, time_n);
    keep = isfinite(y);
    if nnz(keep) < 2 || numel(unique(animal_id(keep))) < 2
        continue
    end
    tbl = table(y(keep), categorical(animal_id(keep)), ...
        'VariableNames', {'Voltage', 'Animal'});
    try
        lme = fitlme(tbl, 'Voltage ~ 1 + (1|Animal)');
        row = intercept_row(lme.Coefficients);
        if isempty(row)
            error('grand_average_erp_lme:intercept', ...
                'The fit did not return an intercept.');
        end
        ci_t = coefCI(lme, 'Alpha', alpha);
        estimate(time_n) = lme.Coefficients.Estimate(row);
        se(time_n) = lme.Coefficients.SE(row);
        pvals(time_n) = lme.Coefficients.pValue(row);
        ci(time_n, :) = ci_t(row, :);
    catch fit_error
        n_failed = n_failed + 1;
        if n_failed == 1
            first_error = fit_error.message;
        end
    end
end

if n_failed > 0
    warning('grand_average_erp_lme:fit', ...
        '%d of %d time points were not fit. First error: %s', ...
        n_failed, n_time, first_error);
end
if ~any(isfinite(estimate))
    error('grand_average_erp_lme:fit', ...
        'The mixed model did not fit at any time point. First error: %s', ...
        first_error);
end

results = struct();
results.time_s = time_s;
results.mean = estimate;
results.se = se;
results.ci = ci;
results.pvals = pvals;
results.alpha = alpha;
results.formula = 'Voltage ~ 1 + (1|Animal)';
results.session_erp = session_erp;
results.animal_id = animal_id;
results.n_sessions = size(session_erp, 1);
results.n_animals = numel(unique(animal_id));
results.n_failed = n_failed;
end

function animal_id = animal_labels(animal_id)
if isstring(animal_id)
    animal_id = cellstr(animal_id(:));
elseif ischar(animal_id)
    animal_id = cellstr(animal_id);
elseif iscell(animal_id)
    animal_id = animal_id(:);
else
    error('grand_average_erp_lme:animal', ...
        'animal_id must be a cell array of recording-animal labels.');
end
for label_n = 1:numel(animal_id)
    animal_id{label_n} = char(animal_id{label_n});
end
end

function row = intercept_row(coef)
names = coef.Name;
if isstring(names) || iscategorical(names)
    names = cellstr(names);
elseif ~iscell(names)
    names = cellstr(string(names));
end
row = find(strcmp(names, '(Intercept)'), 1);
end
