function packed = neuron_acg_with_nulls(spikes, t0, t1, starts, ends, tail_index, bin_edges, options)
%NEURON_ACG_WITH_NULLS  Whole-session ACG, window ACGs, and tail ACGs.
%
% The observed whole-session autocorrelogram uses every spike in [t0, t1].
% Each null keeps that spike count and draws independent uniform times on
% [t0, t1). The same draw is cut into the fixed power windows to make the
% low-tail and high-tail nulls. A tail is kept only when those windows
% contain at least options.min_spikes spikes. Per-window autocorrelograms
% are observed only, and only in windows that themselves have that many
% spikes. Saved curves are sum-to-1 after the refractory bins are blanked.
% Null row k is shuffle k.

bin_edges = double(bin_edges(:))';
bin_w = bin_edges(2) - bin_edges(1);
centers = bin_edges(1:end - 1) + bin_w / 2;
n_bins = numel(centers);
starts = double(starts(:));
ends = double(ends(:));
spikes = double(spikes(:));
spikes = spikes(isfinite(spikes) & spikes >= t0 & spikes <= t1);
spikes = sort(spikes);
n_spikes = numel(spikes);
duration = double(t1) - double(t0);
n_windows = numel(starts);
counts = window_spike_counts(spikes, t0, options.window_s, n_windows);
names = session_acg_names('curves');
use = gate_tails(tail_index, counts, options.min_spikes, n_windows);

packed = blank_packed(names, n_bins, n_windows);
packed.n_spikes = n_spikes;
packed.duration = duration;
packed.window_spike_count = counts;
packed.tail_used = use;
if duration > 0 && n_spikes > 0
    packed.firing_rate = n_spikes / duration;
end
if n_spikes < 2 || ~(duration > 0)
    return
end

previous_rng = rng;
rng(options.seed, 'twister');
rng_cleanup = onCleanup(@() rng(previous_rng)); %#ok<NASGU>

session_rate = interval_spike_acg(spikes, t0, t1, bin_edges, true);
[session_row, session_mass] = sum_to_one_acg(session_rate, centers, options.acg_refractory_s);
packed.curves.session = observed_curve(session_row, n_spikes, duration, session_mass, options.n_shuffles, n_bins);

for c = 2:numel(names)
    name = names{c};
    if ~use.(name)
        continue
    end
    tail_duration = double(sum(tail_index.(name))) * options.window_s;
    tail_spikes = sum(counts(tail_index.(name)));
    tail_rate = pool_interval_acgs(spikes, starts(tail_index.(name)), ends(tail_index.(name)), bin_edges);
    [tail_row, tail_mass] = sum_to_one_acg(tail_rate, centers, options.acg_refractory_s);
    packed.curves.(name) = observed_curve(tail_row, tail_spikes, tail_duration, tail_mass, options.n_shuffles, n_bins);
end

for w = 1:n_windows
    if counts(w) < options.min_spikes
        continue
    end
    window_rate = interval_spike_acg(spikes, starts(w), ends(w), bin_edges, false);
    packed.window_acg(w, :) = single(sum_to_one_acg(window_rate, centers, options.acg_refractory_s));
end

label = '';
if isfield(options, 'label') && ~isempty(options.label)
    label = options.label;
end
log_file = '';
if isfield(options, 'log_file')
    log_file = options.log_file;
end
clock_start = tic;
for s = 1:options.n_shuffles
    drawn = sort(t0 + duration * rand(n_spikes, 1));
    drawn_rate = interval_spike_acg(drawn, t0, t1, bin_edges, true);
    packed.curves.session.null(s, :) = single(sum_to_one_acg( ...
        drawn_rate, centers, options.acg_refractory_s));
    for c = 2:numel(names)
        name = names{c};
        if ~use.(name)
            continue
        end
        drawn_tail = pool_interval_acgs(drawn, starts(tail_index.(name)), ...
            ends(tail_index.(name)), bin_edges);
        packed.curves.(name).null(s, :) = single(sum_to_one_acg( ...
            drawn_tail, centers, options.acg_refractory_s));
    end
    if mod(s, 25) == 0 && ~isempty(label)
        session_acg_log(log_file, sprintf('%s shuffle %d/%d (%.0f s)', ...
            label, s, options.n_shuffles, toc(clock_start)));
    end
end

for c = 1:numel(names)
    name = names{c};
    packed.scores.(name) = score_acg_oscillation( ...
        packed.curves.(name).observed, packed.curves.(name).null, ...
        centers, options.smooth_n, options.oscillation_alpha);
end
end


function use = gate_tails(tail_index, counts, min_spikes, n_windows)
names = session_acg_names('curves');
names = names(2:end);
use = struct();
for i = 1:numel(names)
    name = names{i};
    use.(name) = false;
    if ~isfield(tail_index, name) || isempty(tail_index.(name))
        continue
    end
    index = tail_index.(name)(:);
    if numel(index) ~= n_windows || ~any(index)
        continue
    end
    if sum(counts(index)) >= min_spikes
        use.(name) = true;
    end
end
end


function packed = blank_packed(names, n_bins, n_windows)
packed.n_spikes = 0;
packed.duration = NaN;
packed.firing_rate = NaN;
packed.window_spike_count = zeros(n_windows, 1);
packed.window_acg = nan(n_windows, n_bins, 'single');
packed.tail_used = struct();
packed.curves = struct();
packed.scores = struct();
for i = 1:numel(names)
    packed.curves.(names{i}) = missing_curve(n_bins);
    packed.scores.(names{i}) = score_acg_oscillation(nan(1, n_bins), [], [], 1, 1);
    packed.tail_used.(names{i}) = false;
end
packed.tail_used.session = true;
end


function curve = missing_curve(n_bins)
curve.observed = nan(1, n_bins, 'single');
curve.null = nan(0, n_bins, 'single');
curve.n_spikes = 0;
curve.duration = NaN;
curve.mass = NaN;
end


function curve = observed_curve(row, n_spikes, duration, mass, n_shuffles, n_bins)
curve.observed = single(row);
curve.null = nan(n_shuffles, n_bins, 'single');
curve.n_spikes = n_spikes;
curve.duration = duration;
curve.mass = mass;
end
