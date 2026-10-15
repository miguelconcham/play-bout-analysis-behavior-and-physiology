%% Plot_delta_peak_trough_evoked
% Load delta peak/trough LFP stacks written by Estimate_delta_peak_trough_csd
% and plot trial-averaged evoked responses: peak | trough as two subplots,
% one figure per session (animal, or animal/probeN for NPX2).
%
% Memory: sessions are processed one at a time. For each session the peak
% stack is loaded → averaged → cleared, then the trough stack, then plot.
% Large arrays are not kept across sessions.
%
% Edit root_folder if needed, then run this script.

%% User settings
root_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\delta peak trough CSD';
% Fallback if you still use the .pri host:
% root_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\delta peak trough CSD';

save_figures   = false;  % true → PNG under root_folder/figures, then close figure
close_after    = false;  % true → close each figure after plotting (saves RAM/UI)
clim_uV        = [];     % e.g. [-200 200]; empty = autoscale per panel
sort_by_depth  = true;

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = fullfile(root_folder, '..', '..', '..', 'CSD analysis', 'Plot_delta_peak_trough_evoked.m');
end
csd_code_dir = fileparts(this_file);
addpath(csd_code_dir);

peak_root   = fullfile(root_folder, 'peak');
trough_root = fullfile(root_folder, 'trough');
if ~exist(peak_root, 'dir')
    error('Plot_delta_peak_trough_evoked:Root', 'Missing peak/ under %s', root_folder);
end

sessions = local_list_sessions(peak_root);
if isempty(sessions)
    error('Plot_delta_peak_trough_evoked:Empty', 'No session folders with lfp.npy under %s', peak_root);
end
disp(['Found ', num2str(numel(sessions)), ' session(s) under peak/']);
disp('Loading one session at a time (peak then trough), then clearing.');

for s = 1:numel(sessions)
    rel = sessions{s};
    peak_dir   = fullfile(peak_root, rel);
    trough_dir = fullfile(trough_root, rel);

    if ~exist(fullfile(peak_dir, 'lfp.npy'), 'file')
        warning('Skipping %s — no peak lfp.npy', rel);
        continue
    end
    if ~exist(fullfile(trough_dir, 'lfp.npy'), 'file')
        warning('Skipping %s — no matching trough folder', rel);
        continue
    end

    disp(['[', num2str(s), '/', num2str(numel(sessions)), '] ', rel, ' — peak']);
    [ev_peak, t_ms_p, y_p, n_p] = local_load_evoked(peak_dir, sort_by_depth);

    disp(['[', num2str(s), '/', num2str(numel(sessions)), '] ', rel, ' — trough']);
    [ev_tr, t_ms_t, y_t, n_t] = local_load_evoked(trough_dir, sort_by_depth);

    fig = figure('Name', rel, 'Color', 'w', 'Position', [80 80 1200 500]);
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    local_plot_panel(ev_peak, t_ms_p, y_p, clim_uV);
    title(sprintf('Peak  (n=%d)', n_p));
    ylabel('Depth (\mum)');
    xlabel('Time (ms)');

    nexttile;
    local_plot_panel(ev_tr, t_ms_t, y_t, clim_uV);
    title(sprintf('Trough  (n=%d)', n_t));
    xlabel('Time (ms)');

    sgtitle(rel, 'Interpreter', 'none');
    drawnow;

    if save_figures
        out_dir = fullfile(root_folder, 'figures');
        if ~exist(out_dir, 'dir'); mkdir(out_dir); end
        out_name = strrep(rel, filesep, '_');
        saveas(fig, fullfile(out_dir, [out_name, '_evoked.png']));
    end

    % Drop large arrays before the next session
    clear ev_peak ev_tr t_ms_p t_ms_t y_p y_t n_p n_t

    if close_after || save_figures
        close(fig);
        clear fig
    end
end

disp('Done plotting peak/trough evoked LFP.');

%% -------------------------------------------------------------------------
function sessions = local_list_sessions(peak_root)
sessions = {};
req = 'lfp.npy';
stack = {peak_root};
while ~isempty(stack)
    cur = stack{end};
    stack(end) = [];
    d = dir(cur);
    d = d(~ismember({d.name}, {'.', '..'}));
    subdirs = d([d.isdir]);
    files = d(~[d.isdir]);
    if any(strcmp({files.name}, req))
        rel = strrep(cur, [peak_root, filesep], '');
        if strcmp(rel, cur)  % peak_root itself
            rel = '.';
        end
        if ~strcmp(rel, '.')
            sessions{end + 1} = rel; %#ok<AGROW>
        end
    end
    for k = 1:numel(subdirs)
        stack{end + 1} = fullfile(cur, subdirs(k).name); %#ok<AGROW>
    end
end
sessions = sort(sessions);
end

function [evoked, time_ms, y, n_events] = local_load_evoked(folder, sort_by_depth)
% Load one lfp.npy, average over events immediately, clear the stack.

lfp = read_csd_npy(fullfile(folder, 'lfp.npy'));  % keep native (usually single)
x   = double(read_csd_npy(fullfile(folder, 'x_coords.npy'))); x = x(:);
y   = double(read_csd_npy(fullfile(folder, 'y_coords.npy'))); y = y(:);
T0  = double(read_csd_npy(fullfile(folder, 'T0.npy'))); T0 = T0(1);
if exist(fullfile(folder, 'fs.npy'), 'file')
    fs = double(read_csd_npy(fullfile(folder, 'fs.npy'))); fs = fs(1);
else
    fs = 2500;
end

if ndims(lfp) ~= 3
    error('Expected lfp.npy as (n_events, n_chan, n_time), got %s', mat2str(size(lfp)));
end

n_events = size(lfp, 1);
evoked = double(squeeze(mean(lfp, 1)));  % (n_chan, n_time) only
clear lfp   % free the full event stack before returning

if sort_by_depth
    [y, ord] = sort(y);
    evoked = evoked(ord, :);
end

n_time = size(evoked, 2);
time_ms = ((0:n_time - 1) - T0) / fs * 1000;
end

function local_plot_panel(evoked, time_ms, y, clim_uV)
imagesc(time_ms, y, evoked);
axis xy tight;
colormap(gca, 'parula');
cb = colorbar;
cb.Label.String = 'LFP (a.u.)';
if ~isempty(clim_uV)
    caxis(clim_uV);
end
hold on;
plot([0 0], ylim, 'k-');
hold off;
end
