function report = distance_coverage_report(sessions, options)
%DISTANCE_COVERAGE_REPORT  Play vs non-play coverage along inter-animal distance.
%
% sessions is the struct array from distance_coverage_sessions. options
% defaults to distance_coverage_options. Durations are frame counts times
% the spectrogram hop, the same span used for the high-speed coverage.
%
% A distance bin is [left, right) except the last bin, which includes its
% right edge. A session counts in a bin when the play frames in that bin
% and the non-play frames in that bin each span the bin minimum.
% The far tail is distance z >= high_z. The close tail is z <= -high_z.

if nargin < 2 || isempty(options)
    options = distance_coverage_options();
end
validate_options(options);

if isempty(sessions)
    sessions = repmat(struct( ...
        'animal', '', 'partner', NaN, 'distance_z', [], 'power', [], ...
        'is_play', [], 'hop_s', NaN, 'wind_s', NaN), 0, 1);
end

n_sess = numel(sessions);
animals = cell(n_sess, 1);
partners = nan(n_sess, 1);
hop_s = nan(n_sess, 1);
wind_s = nan(n_sess, 1);
nonplay_far_s = zeros(n_sess, 1);
play_far_s = zeros(n_sess, 1);
enough_far = false(n_sess, 1);
diff_far = nan(n_sess, 1);
diff_far_z = nan(n_sess, 1);
nonplay_close_s = zeros(n_sess, 1);
play_close_s = zeros(n_sess, 1);
enough_close = false(n_sess, 1);
diff_close = nan(n_sess, 1);
diff_close_z = nan(n_sess, 1);
max_nonplay = nan(n_sess, 1);
max_play = nan(n_sess, 1);
min_nonplay = nan(n_sess, 1);
min_play = nan(n_sess, 1);
far_bin_lo = nan(n_sess, 1);
far_bin_hi = nan(n_sess, 1);
far_bin_play_s = nan(n_sess, 1);
far_bin_nonplay_s = nan(n_sess, 1);
diff_far_bin = nan(n_sess, 1);
diff_far_bin_z = nan(n_sess, 1);
close_bin_lo = nan(n_sess, 1);
close_bin_hi = nan(n_sess, 1);
close_bin_play_s = nan(n_sess, 1);
close_bin_nonplay_s = nan(n_sess, 1);
diff_close_bin = nan(n_sess, 1);
diff_close_bin_z = nan(n_sess, 1);
X_cell = cell(n_sess, 1);
Y_cell = cell(n_sess, 1);
play_cell = cell(n_sess, 1);
observed_min = Inf;
observed_max = -Inf;

for s = 1:n_sess
    [X, Y, is_play, hop_s(s), wind_s(s), animals{s}, partners(s)] = ...
        session_vectors(sessions(s), s);
    if isempty(X)
        continue
    end
    observed_min = min(observed_min, min(X));
    observed_max = max(observed_max, max(X));
    X_cell{s} = X;
    Y_cell{s} = Y;
    play_cell{s} = is_play;

    session_sd = std(Y);
    far = X >= options.high_z;
    close = X <= -options.high_z;
    nonplay_far_s(s) = nnz(far & ~is_play) * hop_s(s);
    play_far_s(s) = nnz(far & is_play) * hop_s(s);
    enough_far(s) = nonplay_far_s(s) >= options.min_nonplay_s ...
        && play_far_s(s) >= options.min_play_s;
    [diff_far(s), diff_far_z(s)] = group_difference(Y, far, is_play, session_sd);

    nonplay_close_s(s) = nnz(close & ~is_play) * hop_s(s);
    play_close_s(s) = nnz(close & is_play) * hop_s(s);
    enough_close(s) = nonplay_close_s(s) >= options.min_nonplay_s ...
        && play_close_s(s) >= options.min_play_s;
    [diff_close(s), diff_close_z(s)] = group_difference(Y, close, is_play, session_sd);

    if any(~is_play)
        max_nonplay(s) = max(X(~is_play));
        min_nonplay(s) = min(X(~is_play));
    end
    if any(is_play)
        max_play(s) = max(X(is_play));
        min_play(s) = min(X(is_play));
    end

    [far_bin_lo(s), far_bin_hi(s), far_bin_play_s(s), far_bin_nonplay_s(s), ...
        diff_far_bin(s)] = supported_bin(X, Y, is_play, hop_s(s), ...
        options.ceiling_bin_width, options.min_ceiling_play_s, ...
        options.min_ceiling_nonplay_s, true);
    [close_bin_lo(s), close_bin_hi(s), close_bin_play_s(s), ...
        close_bin_nonplay_s(s), diff_close_bin(s)] = supported_bin( ...
        X, Y, is_play, hop_s(s), options.ceiling_bin_width, ...
        options.min_ceiling_play_s, options.min_ceiling_nonplay_s, false);
    if session_sd > 0
        if isfinite(diff_far_bin(s))
            diff_far_bin_z(s) = diff_far_bin(s) / session_sd;
        end
        if isfinite(diff_close_bin(s))
            diff_close_bin_z(s) = diff_close_bin(s) / session_sd;
        end
    end
end

if isfinite(observed_min)
    curve_edges = coverage_curve_edges(observed_min, observed_max, ...
        options.curve_bin_width, options.high_z, options.curve_z_limit);
else
    curve_edges = coverage_curve_edges([], [], ...
        options.curve_bin_width, options.high_z, options.curve_z_limit);
end

use_all = true(n_sess, 1);
[curve_x, diff_all, n_both, n_play, n_nonplay] = coverage_curves( ...
    X_cell, Y_cell, play_cell, hop_s, curve_edges, ...
    options.min_bin_play_s, options.min_bin_nonplay_s, use_all);
[~, diff_far_curve, n_both_far, n_play_far, n_nonplay_far] = coverage_curves( ...
    X_cell, Y_cell, play_cell, hop_s, curve_edges, ...
    options.min_bin_play_s, options.min_bin_nonplay_s, enough_far);
[~, diff_close_curve, n_both_close, n_play_close, n_nonplay_close] = coverage_curves( ...
    X_cell, Y_cell, play_cell, hop_s, curve_edges, ...
    options.min_bin_play_s, options.min_bin_nonplay_s, enough_close);

report = struct( ...
    'options', options, ...
    'n_sessions', n_sess, ...
    'animals', {animals}, ...
    'partners', partners, ...
    'hop_s', hop_s, ...
    'wind_s', wind_s, ...
    'nonplay_far_s', nonplay_far_s, ...
    'play_far_s', play_far_s, ...
    'enough_far', enough_far, ...
    'diff_far', diff_far, ...
    'diff_far_z', diff_far_z, ...
    'nonplay_close_s', nonplay_close_s, ...
    'play_close_s', play_close_s, ...
    'enough_close', enough_close, ...
    'diff_close', diff_close, ...
    'diff_close_z', diff_close_z, ...
    'max_nonplay', max_nonplay, ...
    'max_play', max_play, ...
    'min_nonplay', min_nonplay, ...
    'min_play', min_play, ...
    'far_bin_lo', far_bin_lo, ...
    'far_bin_hi', far_bin_hi, ...
    'far_bin_play_s', far_bin_play_s, ...
    'far_bin_nonplay_s', far_bin_nonplay_s, ...
    'diff_far_bin', diff_far_bin, ...
    'diff_far_bin_z', diff_far_bin_z, ...
    'close_bin_lo', close_bin_lo, ...
    'close_bin_hi', close_bin_hi, ...
    'close_bin_play_s', close_bin_play_s, ...
    'close_bin_nonplay_s', close_bin_nonplay_s, ...
    'diff_close_bin', diff_close_bin, ...
    'diff_close_bin_z', diff_close_bin_z, ...
    'curve_edges', curve_edges, ...
    'curve_x', curve_x, ...
    'diff_all', diff_all, ...
    'n_both', n_both, ...
    'n_play', n_play, ...
    'n_nonplay', n_nonplay, ...
    'diff_far_curve', diff_far_curve, ...
    'n_both_far', n_both_far, ...
    'n_play_far', n_play_far, ...
    'n_nonplay_far', n_nonplay_far, ...
    'diff_close_curve', diff_close_curve, ...
    'n_both_close', n_both_close, ...
    'n_play_close', n_play_close, ...
    'n_nonplay_close', n_nonplay_close);
end

function [X, Y, is_play, hop_s, wind_s, animal, partner] = session_vectors(session, index)
X = double(session.distance_z(:));
Y = double(session.power(:));
is_play = session.is_play(:);
animal = session.animal;
partner = session.partner;
hop_s = session.hop_s;
if isfield(session, 'wind_s') && ~isempty(session.wind_s)
    wind_s = session.wind_s;
else
    wind_s = 0.25;
end
if ~islogical(is_play)
    if isnumeric(is_play)
        is_play = is_play ~= 0 & isfinite(is_play);
    else
        error('distance_coverage_report:PlayMask', ...
            'Session %d is_play must be logical or numeric.', index);
    end
end
if numel(X) ~= numel(Y) || numel(is_play) ~= numel(X)
    error('distance_coverage_report:LengthMismatch', ...
        'Session %d distance, power, and play mask differ in length.', index);
end
if ~isscalar(hop_s) || ~isfinite(hop_s) || hop_s <= 0
    error('distance_coverage_report:Hop', ...
        'Session %d needs a positive frame hop in seconds.', index);
end
keep = isfinite(X) & isfinite(Y);
X = X(keep);
Y = Y(keep);
is_play = is_play(keep);
if exist('isstring', 'builtin') == 5 && isstring(animal)
    animal = char(animal);
end
animal = char(animal);
end

function [diff_val, diff_z] = group_difference(Y, in_tail, is_play, session_sd)
diff_val = NaN;
diff_z = NaN;
if any(in_tail & is_play) && any(in_tail & ~is_play)
    diff_val = mean(Y(in_tail & is_play)) - mean(Y(in_tail & ~is_play));
    if session_sd > 0
        diff_z = diff_val / session_sd;
    end
end
end

function [lo, hi, play_s, nonplay_s, diff_val] = supported_bin( ...
    X, Y, is_play, hop_s, bin_width, min_play_s, min_nonplay_s, from_high)

lo = NaN;
hi = NaN;
play_s = NaN;
nonplay_s = NaN;
diff_val = NaN;
if isempty(X)
    return
end
edges = span_edges(min(X), max(X), bin_width);
bin_order = 1:(numel(edges) - 1);
if from_high
    bin_order = fliplr(bin_order);
end
for b = bin_order
    in_bin = samples_in_bin(X, edges, b);
    this_play_s = nnz(in_bin & is_play) * hop_s;
    this_nonplay_s = nnz(in_bin & ~is_play) * hop_s;
    if this_play_s >= min_play_s && this_nonplay_s >= min_nonplay_s
        lo = edges(b);
        hi = edges(b + 1);
        play_s = this_play_s;
        nonplay_s = this_nonplay_s;
        diff_val = mean(Y(in_bin & is_play)) - mean(Y(in_bin & ~is_play));
        return
    end
end
end

function edges = span_edges(low_value, high_value, bin_width)
edge_lo = floor(low_value / bin_width) * bin_width;
edge_hi = ceil(high_value / bin_width) * bin_width;
if edge_hi <= edge_lo
    edge_hi = edge_lo + bin_width;
end
n_edges = round((edge_hi - edge_lo) / bin_width) + 1;
edges = linspace(edge_lo, edge_hi, n_edges);
end

function edges = coverage_curve_edges(low_value, high_value, bin_width, high_z, z_limit)
if isempty(low_value) || ~isfinite(low_value) || ~isfinite(high_value)
    lo = -high_z;
    hi = high_z + 1;
else
    lo = min(-high_z, low_value);
    hi = max(high_z + 1, high_value);
end
lo = max(-z_limit, floor(lo / bin_width) * bin_width);
hi = min(z_limit, ceil(hi / bin_width) * bin_width);
if hi <= lo
    hi = lo + bin_width;
end
n_edges = round((hi - lo) / bin_width) + 1;
edges = linspace(lo, hi, max(n_edges, 2));
end

function [centers, diff_mat, n_both, n_play, n_nonplay] = coverage_curves( ...
    X_cell, Y_cell, play_cell, hop_s, edges, min_play_s, min_nonplay_s, use_session)

centers = edges(1:end - 1) + diff(edges) / 2;
n_sess = numel(X_cell);
n_bins = numel(centers);
diff_mat = nan(n_sess, n_bins);
n_both = zeros(1, n_bins);
n_play = zeros(1, n_bins);
n_nonplay = zeros(1, n_bins);
for s = 1:n_sess
    if ~use_session(s) || isempty(X_cell{s})
        continue
    end
    X = X_cell{s};
    Y = Y_cell{s};
    is_play = play_cell{s};
    for b = 1:n_bins
        in_bin = samples_in_bin(X, edges, b);
        play_s = nnz(in_bin & is_play) * hop_s(s);
        nonplay_s = nnz(in_bin & ~is_play) * hop_s(s);
        play_ok = play_s >= min_play_s;
        nonplay_ok = nonplay_s >= min_nonplay_s;
        if play_ok
            n_play(b) = n_play(b) + 1;
        end
        if nonplay_ok
            n_nonplay(b) = n_nonplay(b) + 1;
        end
        if play_ok && nonplay_ok
            diff_mat(s, b) = mean(Y(in_bin & is_play)) - mean(Y(in_bin & ~is_play));
            n_both(b) = n_both(b) + 1;
        end
    end
end
end

function in_bin = samples_in_bin(X, edges, bin_index)
if bin_index == numel(edges) - 1
    in_bin = X >= edges(bin_index) & X <= edges(bin_index + 1);
else
    in_bin = X >= edges(bin_index) & X < edges(bin_index + 1);
end
end

function validate_options(options)
needed = {'high_z', 'ceiling_bin_width', 'min_nonplay_s', 'min_play_s', ...
    'min_ceiling_nonplay_s', 'min_ceiling_play_s', 'curve_bin_width', ...
    'min_bin_nonplay_s', 'min_bin_play_s', 'curve_z_limit'};
for k = 1:numel(needed)
    if ~isfield(options, needed{k})
        error('distance_coverage_report:Options', ...
            'Missing option %s.', needed{k});
    end
end
if ~isfinite(options.high_z) || options.high_z <= 0
    error('distance_coverage_report:Options', 'high_z must be positive.');
end
widths = [options.ceiling_bin_width, options.curve_bin_width];
if any(~isfinite(widths) | widths <= 0)
    error('distance_coverage_report:Options', 'Bin widths must be positive.');
end
if ~isfinite(options.curve_z_limit) || options.curve_z_limit <= options.high_z
    error('distance_coverage_report:Options', ...
        'curve_z_limit must be greater than high_z.');
end
mins = [options.min_nonplay_s, options.min_play_s, ...
    options.min_ceiling_nonplay_s, options.min_ceiling_play_s, ...
    options.min_bin_nonplay_s, options.min_bin_play_s];
if any(~isfinite(mins) | mins < 0)
    error('distance_coverage_report:Options', ...
        'Duration minima must be non-negative.');
end
end
