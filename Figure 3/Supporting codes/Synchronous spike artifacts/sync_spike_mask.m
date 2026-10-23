function mask = sync_spike_mask(spikes, partner, bout_starts, hist_range, tolerance)
%SYNC_SPIKE_MASK  Spikes that coincide with a partner inside play-bout windows.
%
% spikes and partner are sorted spike times in seconds. A spike is marked
% when any partner spike falls within +/- tolerance and both spikes lie in
% [bout_start + hist_range(1), bout_start + hist_range(2)]. The same
% definition is used by GENERATE_PHASE_COINCIDENCE_STRUCTURE.

spikes = double(spikes(:));
partner = double(partner(:));
bout_starts = double(bout_starts(:));
mask = false(numel(spikes), 1);
if isempty(spikes) || isempty(partner) || isempty(bout_starts)
    return
end
if numel(hist_range) ~= 2
    error('sync_spike_mask:window', 'hist_range must be [start end] in seconds.');
end

for bout_index = 1:numel(bout_starts)
    t0 = bout_starts(bout_index);
    window_start = t0 + hist_range(1);
    window_end = t0 + hist_range(2);
    in_window = spikes >= window_start & spikes <= window_end;
    partner_window = partner(partner >= window_start & partner <= window_end);
    if ~any(in_window) || isempty(partner_window)
        continue
    end
    window_index = find(in_window);
    hit = mark_synchronous(spikes(window_index), partner_window, tolerance);
    mask(window_index(hit)) = true;
end
end

function hit = mark_synchronous(spikes, partner, tolerance)
n_spikes = numel(spikes);
n_partner = numel(partner);
hit = false(n_spikes, 1);
partner_start = 1;
for spike_index = 1:n_spikes
    while partner_start <= n_partner && partner(partner_start) < spikes(spike_index) - tolerance
        partner_start = partner_start + 1;
    end
    partner_index = partner_start;
    while partner_index <= n_partner && partner(partner_index) <= spikes(spike_index) + tolerance
        hit(spike_index) = true;
        partner_index = partner_index + 1;
    end
end
end
