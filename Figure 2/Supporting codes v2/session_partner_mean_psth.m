function obs = session_partner_mean_psth(power, partner, session_id, bout_length, time, zscore_limit, length_limit)
%SESSION_PARTNER_MEAN_PSTH  Mean play-bout PSTH for each recording day and partner.
%
%   power is nBout x nTime. partner and bout_length are one value per bout.
%   session_id names the recording day of each bout. A sample enters the mean
%   when its time is before the bout ends, its absolute value is within
%   zscore_limit, and the bout is at least length_limit seconds. Only partners
%   1 and 2 are kept.

if nargin < 7 || isempty(length_limit)
    length_limit = 0;
end
if nargin < 6 || isempty(zscore_limit)
    zscore_limit = 4;
end

power = double(power);
partner = double(partner(:));
bout_length = double(bout_length(:));
time = double(time(:)');
session_id = session_labels(session_id);

n_bouts = size(power, 1);
if numel(partner) ~= n_bouts || numel(bout_length) ~= n_bouts || numel(session_id) ~= n_bouts
    error('session_partner_mean_psth:Size', ...
        'power, partner, bout_length, and session_id must have one row per bout.');
end
if size(power, 2) ~= numel(time)
    error('session_partner_mean_psth:Time', ...
        'time must have one entry per power column.');
end

keep = (partner == 1 | partner == 2) & bout_length >= length_limit & isfinite(bout_length);
keep = keep & ~cellfun(@isempty, session_id);

obs.time = time;
obs.session_id = cell(0, 1);
obs.partner = zeros(0, 1);
obs.power = zeros(0, numel(time));
obs.n_bouts = zeros(0, numel(time));
if ~any(keep)
    return
end

sessions = unique(session_id(keep), 'stable');
for s = 1:numel(sessions)
    for p = [1 2]
        trial = keep & strcmp(session_id, sessions{s}) & partner == p;
        if ~any(trial)
            continue
        end
        block = power(trial, :);
        in_bout = time < bout_length(trial);
        block(~in_bout) = NaN;
        block(abs(block) > zscore_limit) = NaN;
        n_at_time = sum(isfinite(block), 1);
        if ~any(n_at_time)
            continue
        end
        obs.session_id(end + 1, 1) = sessions(s); %#ok<AGROW>
        obs.partner(end + 1, 1) = p; %#ok<AGROW>
        obs.power(end + 1, :) = mean(block, 1, 'omitnan'); %#ok<AGROW>
        obs.n_bouts(end + 1, :) = n_at_time; %#ok<AGROW>
    end
end
end

function session_id = session_labels(session_id)
if iscell(session_id)
    session_id = session_id(:);
    for k = 1:numel(session_id)
        value = session_id{k};
        if ischar(value)
            continue
        elseif is_text_scalar(value)
            session_id{k} = char(value);
        else
            error('session_partner_mean_psth:Session', ...
                'Each session id must be the recording-day name.');
        end
    end
elseif ischar(session_id)
    session_id = cellstr(session_id);
elseif is_text_scalar(session_id)
    session_id = cellstr(session_id);
else
    error('session_partner_mean_psth:Session', ...
        'session_id must name the recording day of each bout.');
end
session_id = session_id(:);
end

function tf = is_text_scalar(value)
tf = false;
if exist('isstring', 'builtin') == 5 || exist('isstring', 'file') == 2
    tf = isstring(value);
end
end
