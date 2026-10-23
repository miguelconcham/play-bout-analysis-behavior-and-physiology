%% Plot_delta_peak_trough_evoked
% Inspect delta-aligned LFP before running Python CSD.
%
% For each session under:
%   <root_folder>/peak/<session>/
%   <root_folder>/trough/<session>/
%
% Figure layout (2 x 2):
%   Top:    trial-averaged LFP (channels x time) for peak | trough
%   Bottom: LFP amplitude near t = 0 vs depth, event-wise mean ± 95% CI
%
% Memory: one session at a time. Peak stack loaded → stats → cleared,
% then trough, then plot / save / pause / clear.
%
% Edit root_folder if needed, then run this script.

%% User settings
root_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\delta peak trough CSD';

save_figures     = true;   % write .fig (and optional PNG) under root_folder/figures
save_png         = false;  % also write PNG next to the .fig
pause_sec        = 2;      % pause after each figure so you can inspect (0 = no pause)
close_after      = false;  % true → close figure after save/pause (frees UI memory)
clim_uV          = [];     % e.g. [-200 200] for imagesc; [] = autoscale
sort_by_depth    = true;
amp_half_win_ms  = 10;     % average |t| <= this (ms) around T0 for amplitude vs depth
ci_level         = 0.95;   % confidence level for event-wise CI (normal approx via SEM)

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file)
    % Editable fallback if mfilename is empty (unsaved / pasted run)
    this_file = fullfile(root_folder, '..', '..', '..', 'Figure 2', 'Supporting codes v2', 'CSD analysis', 'Plot_delta_peak_trough_evoked.m');
end
csd_code_dir = fileparts(this_file);
addpath(csd_code_dir);

peak_root   = fullfile(root_folder, 'peak');
trough_root = fullfile(root_folder, 'trough');
if ~exist(peak_root, 'dir')
    error('Plot_delta_peak_trough_evoked:Root', 'Missing peak/ under %s', root_folder);
end

% Discover session relative paths that have lfp.npy under peak/
sessions = local_list_sessions(peak_root);
if isempty(sessions)
    error('Plot_delta_peak_trough_evoked:Empty', 'No session folders with lfp.npy under %s', peak_root);
end
disp(['Found ', num2str(numel(sessions)), ' session(s) under peak/']);
disp('Processing one session at a time (load → plot → clear).');

out_dir = fullfile(root_folder, 'figures');
if save_figures && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

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

    % ---- Peak: load stack, average map + event-wise amp@0, clear stack ----
    disp(['[', num2str(s), '/', num2str(numel(sessions)), '] ', rel, ' — peak']);
    peak_stats = local_load_session_stats(peak_dir, sort_by_depth, amp_half_win_ms, ci_level);

    % ---- Trough: same ----
    disp(['[', num2str(s), '/', num2str(numel(sessions)), '] ', rel, ' — trough']);
    trough_stats = local_load_session_stats(trough_dir, sort_by_depth, amp_half_win_ms, ci_level);

    % ---- Figure: top = evoked maps; bottom = amp vs depth with CI ----
    fig = figure('Name', rel, 'Color', 'w', 'Position', [60 60 1200 820]);
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Top-left: peak evoked LFP
    nexttile;
    local_plot_evoked(peak_stats.evoked, peak_stats.time_ms, peak_stats.y, clim_uV);
    title(sprintf('Peak evoked LFP (n=%d events)', peak_stats.n_events));
    ylabel('Depth (\mum)');
    xlabel('Time (ms)');

    % Top-right: trough evoked LFP
    nexttile;
    local_plot_evoked(trough_stats.evoked, trough_stats.time_ms, trough_stats.y, clim_uV);
    title(sprintf('Trough evoked LFP (n=%d events)', trough_stats.n_events));
    xlabel('Time (ms)');

    % Bottom-left: peak amplitude near t=0 vs depth
    nexttile;
    local_plot_amp_vs_depth(peak_stats);
    title(sprintf('Peak amp near t=0 (±%g ms), mean ± %.0f%% CI', ...
        amp_half_win_ms, 100 * ci_level));
    ylabel('Depth (\mum)');
    xlabel('LFP amplitude (a.u.)');

    % Bottom-right: trough amplitude near t=0 vs depth
    nexttile;
    local_plot_amp_vs_depth(trough_stats);
    title(sprintf('Trough amp near t=0 (±%g ms), mean ± %.0f%% CI', ...
        amp_half_win_ms, 100 * ci_level));
    xlabel('LFP amplitude (a.u.)');

    sgtitle(rel, 'Interpreter', 'none');
    drawnow;

    % ---- Save .fig for later MATLAB sessions ----
    if save_figures
        out_name = strrep(rel, filesep, '_');
        fig_path = fullfile(out_dir, [out_name, '_evoked.fig']);
        savefig(fig, fig_path);
        disp(['  saved ', fig_path]);
        if save_png
            saveas(fig, fullfile(out_dir, [out_name, '_evoked.png']));
        end
    end

    % ---- Pause so you can inspect before the next session ----
    if pause_sec > 0
        disp(['  pausing ', num2str(pause_sec), ' s — inspect figure "', rel, '"']);
        pause(pause_sec);
    end

    % Free session data before next iteration
    clear peak_stats trough_stats

    if close_after
        close(fig);
        clear fig
    end
end

disp('Done plotting peak/trough evoked LFP.');
if save_figures
    disp(['Figures (.fig) in ', out_dir]);
end

%% =========================================================================
%  Local helpers
%% =========================================================================

function sessions = local_list_sessions(peak_root)
% Recursively find folders under peak_root that contain lfp.npy.
% Returns relative paths (e.g. 'AnimalName' or 'AnimalName\probe2').
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
        if strcmp(rel, cur)
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

function stats = local_load_session_stats(folder, sort_by_depth, amp_half_win_ms, ci_level)
% Load one lfp.npy stack, compute:
%   - evoked: mean over events (n_chan x n_time)
%   - amp_mean / amp_lo / amp_hi: event-wise amplitude near t=0 vs channel
% Then clear the full event stack.

lfp = read_csd_npy(fullfile(folder, 'lfp.npy'));  % (n_events, n_chan, n_time), usually single
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
n_chan   = size(lfp, 2);
n_time   = size(lfp, 3);

% Time axis (ms) with 0 at the aligned peak/trough
time_ms = ((0:n_time - 1) - T0) / fs * 1000;

% Samples used for "amplitude around time 0"
amp_mask = abs(time_ms) <= amp_half_win_ms;
if ~any(amp_mask)
    % fallback: single sample at T0 (1-based index T0+1)
    amp_mask = false(1, n_time);
    amp_mask(min(max(T0 + 1, 1), n_time)) = true;
end

% Event-wise amplitude per channel: mean LFP in the t≈0 window
% amp_events: (n_events x n_chan)
amp_events = mean(lfp(:, :, amp_mask), 3);

% Trial-averaged map for imagesc
evoked = double(squeeze(mean(lfp, 1)));  % (n_chan, n_time)

% Done with the huge stack
clear lfp

% Event-wise mean and CI across events (per channel)
amp_mean = mean(amp_events, 1)';                 % (n_chan x 1)
amp_sem  = std(amp_events, 0, 1)' / sqrt(max(n_events, 1));
z = local_norminv((1 + ci_level) / 2);           % ~1.96 for 95%
amp_lo = amp_mean - z * amp_sem;
amp_hi = amp_mean + z * amp_sem;

if sort_by_depth
    [y, ord] = sort(y);
    evoked = evoked(ord, :);
    amp_mean = amp_mean(ord);
    amp_lo = amp_lo(ord);
    amp_hi = amp_hi(ord);
    x = x(ord);
end

stats = struct();
stats.evoked   = evoked;
stats.time_ms  = time_ms;
stats.y        = y;
stats.x        = x;
stats.n_events = n_events;
stats.amp_mean = amp_mean;
stats.amp_lo   = amp_lo;
stats.amp_hi   = amp_hi;
stats.fs       = fs;
stats.T0       = T0;
end

function z = local_norminv(p)
% Inverse CDF of standard normal (avoids needing Statistics Toolbox).
% Rational approximation (Abramowitz & Stegun 26.2.23), good for CI use.
if p <= 0 || p >= 1
    error('local_norminv:p', 'p must be in (0,1)');
end
if p == 0.5
    z = 0;
    return
end
% Use erfinv if available (MATLAB core)
z = sqrt(2) * erfinv(2 * p - 1);
end

function local_plot_evoked(evoked, time_ms, y, clim_uV)
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

function local_plot_amp_vs_depth(stats)
% Horizontal mean ± CI at each depth (y on vertical axis to match imagesc).
y = stats.y(:);
m = stats.amp_mean(:);
lo = stats.amp_lo(:);
hi = stats.amp_hi(:);

% CI as a patch (x = amplitude, y = depth)
yy = [y; flipud(y)];
xx = [lo; flipud(hi)];
fill(xx, yy, [0.75 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
plot(m, y, 'b-', 'LineWidth', 1.5);
plot([0 0], [min(y) max(y)], 'k:');
hold off;
axis tight;
set(gca, 'YDir', 'normal');  % same orientation as imagesc(..., axis xy)
end
