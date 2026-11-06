function options = irasa_play_session_options(options)
%IRASA_PLAY_SESSION_OPTIONS  Window settings for the play-session IRASA.
%
% window_seconds is the length of each LFP epoch passed to IRASA. The
% default is 10 s, so the raw spectrum has bins of 0.1 Hz. Pass another
% positive duration to change it.
% window_step_seconds is the start-to-start step. Leave it empty to use
% the same value as window_seconds, which makes the windows contiguous
% and non-overlapping.

if nargin < 1 || isempty(options)
    options = struct();
end
if ~isstruct(options) || numel(options) ~= 1
    error('irasa_play_session_options:input', ...
        'Options must be a scalar struct.');
end

options = fill_option(options, 'window_seconds', 10);
if ~isfield(options, 'window_step_seconds') || isempty(options.window_step_seconds)
    options.window_step_seconds = options.window_seconds;
end

if ~isscalar(options.window_seconds) || ~isfinite(options.window_seconds) || ...
        options.window_seconds <= 0
    error('irasa_play_session_options:window', ...
        'window_seconds must be a positive finite scalar.');
end
if ~isscalar(options.window_step_seconds) || ...
        ~isfinite(options.window_step_seconds) || options.window_step_seconds <= 0
    error('irasa_play_session_options:step', ...
        'window_step_seconds must be a positive finite scalar.');
end
end

function options = fill_option(options, name, value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = value;
end
end
