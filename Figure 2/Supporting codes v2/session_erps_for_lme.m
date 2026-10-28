function assembled = session_erps_for_lme(playbout_erp, options)
%SESSION_ERPS_FOR_LME  Clean session ERPs that enter the grand-average model.
%
% Electrode options.electrode_index is kept (the first mid-LPAG site when
% that value is 1). Saturated trials are removed, each remaining trial is
% baseline-corrected, and one mean ERP is kept per recording. Recordings
% with no remaining trial are left out of assembled.session_erp. Trial
% counts still include those recordings.

if nargin < 2 || isempty(options)
    options = playbout_erp_defaults();
end
defaults = playbout_erp_defaults();
options = fill_missing_options(options, defaults);

if isempty(playbout_erp)
    error('session_erps_for_lme:empty', 'playbout_erp has no recordings.');
end
if ~isscalar(options.electrode_index) || ~isfinite(options.electrode_index)
    error('session_erps_for_lme:electrode', ...
        'electrode_index must be a finite scalar.');
end

selected = playbout_erp([playbout_erp.electrode_index] == options.electrode_index);
if isempty(selected)
    error('session_erps_for_lme:electrode', ...
        'No recording has electrode_index %g.', options.electrode_index);
end

time_s = selected(1).time_s(:).';
if isempty(time_s)
    error('session_erps_for_lme:time', 'The first recording has no time axis.');
end
n_sessions = numel(selected);
session_erp = nan(n_sessions, numel(time_s));
animal_id = cell(n_sessions, 1);
session_id = cell(n_sessions, 1);
n_complete = zeros(n_sessions, 1);
n_saturated = zeros(n_sessions, 1);
n_kept = zeros(n_sessions, 1);
lfp_abs_max = 0;

for session_n = 1:n_sessions
    session = selected(session_n);
    session_time = session.time_s(:).';
    if numel(session_time) ~= numel(time_s) || any(session_time ~= time_s)
        error('session_erps_for_lme:time', ...
            '%s does not share the time axis of the first recording.', ...
            session_label(session));
    end
    if size(session.lfp_trials, 2) ~= numel(time_s)
        error('session_erps_for_lme:trials', ...
            '%s trials do not match its time axis.', session_label(session));
    end

    finite_values = session.lfp_trials(isfinite(session.lfp_trials));
    if ~isempty(finite_values)
        lfp_abs_max = max(lfp_abs_max, max(abs(finite_values)));
    end

    [session_erp(session_n, :), info] = clean_session_erp( ...
        session.lfp_trials, time_s, options.saturation_limit, options.baseline_window);
    animal_id{session_n} = char(session.animal_id);
    session_id{session_n} = char(session.session_id);
    n_complete(session_n) = info.n_complete;
    n_saturated(session_n) = info.n_saturated;
    n_kept(session_n) = info.n_kept;
end

used = n_kept > 0 & all(isfinite(session_erp), 2);
if ~any(used)
    error('session_erps_for_lme:trials', ...
        'Every recording lost all of its trials to saturation or incomplete windows.');
end

assembled = struct();
assembled.time_s = time_s;
assembled.session_erp = session_erp(used, :);
assembled.animal_id = animal_id(used);
assembled.session_id = session_id(used);
assembled.electrode_index = options.electrode_index;
assembled.n_complete = n_complete(used);
assembled.n_saturated = n_saturated(used);
assembled.n_kept = n_kept(used);
assembled.n_trials_complete = sum(n_complete);
assembled.n_trials_saturated = sum(n_saturated);
assembled.n_trials_kept = sum(n_kept);
assembled.n_recordings_used = sum(used);
assembled.n_recordings_selected = n_sessions;
assembled.saturation_limit = options.saturation_limit;
assembled.baseline_window = options.baseline_window;
assembled.lfp_abs_max = lfp_abs_max;
end

function options = fill_missing_options(options, defaults)
fields = fieldnames(defaults);
for field_n = 1:numel(fields)
    name = fields{field_n};
    if ~isfield(options, name)
        options.(name) = defaults.(name);
    end
end
end

function label = session_label(session)
if isfield(session, 'session_id') && ~isempty(session.session_id)
    label = char(session.session_id);
else
    label = 'recording';
end
end
