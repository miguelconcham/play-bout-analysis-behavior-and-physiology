function summarize_sync_coverage(coverages, min_spikes)
%SUMMARIZE_SYNC_COVERAGE  Print how much of each area was observed.
%
% coverages is the coverage struct from each session. The target is met
% when every area, in every session, reached min_neuron_fraction of its
% neurons with at least min_spikes spikes. The totals add those sessions.

if nargin < 2 || isempty(min_spikes)
    min_spikes = 1000;
end
if isempty(coverages)
    return
end

n_area = sum([coverages.n_area]);
n_with_spikes = sum([coverages.n_with_spikes]);
n_covered = sum([coverages.n_covered]);
n_below = sum([coverages.n_below_min_spikes]);
n_without_partner = sum([coverages.n_spiking_without_partner]);
target = coverages(1).target;
fprintf(['Observed %d/%d neurons with >= %d spikes (%s) and %d/%d good units ' ...
    'in SC, DLPAG, LPAG, VLPAG, and DR. Target %.0f%% of the spiking neurons in each area.\n'], ...
    n_covered, n_with_spikes, min_spikes, percent_text(n_covered, n_with_spikes), ...
    n_covered, n_area, 100 * target);

names = {};
area_spikes = [];
area_covered = [];
area_units = [];
for session_index = 1:numel(coverages)
    areas = coverages(session_index).areas;
    for area_index = 1:numel(areas)
        area_name = areas(area_index).area;
        slot = find(strcmp(names, area_name), 1);
        if isempty(slot)
            names{end + 1} = area_name; %#ok<AGROW>
            area_spikes(end + 1) = 0; %#ok<AGROW>
            area_covered(end + 1) = 0; %#ok<AGROW>
            area_units(end + 1) = 0; %#ok<AGROW>
            slot = numel(names);
        end
        area_spikes(slot) = area_spikes(slot) + areas(area_index).n_with_spikes;
        area_covered(slot) = area_covered(slot) + areas(area_index).n_covered;
        area_units(slot) = area_units(slot) + areas(area_index).n_area;
    end
end
for slot = 1:numel(names)
    fprintf('  %s  %d/%d with enough spikes (%s), %d good units\n', ...
        names{slot}, area_covered(slot), area_spikes(slot), ...
        percent_text(area_covered(slot), area_spikes(slot)), area_units(slot));
end
if n_below > 0 || n_without_partner > 0
    fprintf(['  Not observed: %d neurons below %d spikes, %d with enough spikes ' ...
        'but no synchronous partner that can support a null.\n'], ...
        n_below, min_spikes, n_without_partner);
end
end

function label = percent_text(numerator, denominator)
if denominator < 1
    label = 'n/a';
else
    label = sprintf('%.0f%%', 100 * numerator / denominator);
end
end
