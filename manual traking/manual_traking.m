%% manual_traking
% Click the focal animal, then the partner, on frames taken from a video.
% Run this file one section at a time. Each section uses the variables
% created above it. Do not use "Run All": the clicking sections wait for
% one mouse click per sample.
%
% Return in the figure window stops a clicking section and keeps every
% click collected so far.
%
% Workspace / saved struct: traking_structure
%   .frames2stract   video frame index of each click (row)
%   .time            seconds, frames2stract / frame_rate
%   .animal_pos      [n x 2] focal-animal clicks, [x y] in image pixels
%   .partner_pos     [n x 2] partner clicks
%   .outher_border   outer-arena polygon, [x y] vertices
%   .inner_border    inner-arena polygon
%   .inner_box       inner-box polygon
%
% Those field names, including the spellings frames2stract and
% outher_border, are what Figure 1 / Figure 2 load from Traking backups.
% .time follows frame_rate. Several analysis scripts ignore .time and
% divide frames2stract by 30, so frame_rate should be the real video rate
% (30 for the NPX play videos, 15 for the monitoring clip below).

%% Parameters
% Example on the lab file server (left commented; set video_file below):
% video_file = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\20231012\_NMMTT_230929 23-10-12 14-32-19.avi';

video_file = 'control_HC_monitoring_0002 25-06-09 16-09-46.avi';
track_window_s = [310 1290];   % [start end] in seconds, the old ranagetrak
frame_rate = 15;                % frames per second of video_file
samples_per_second = 3;        % clicks per second inside track_window_s
brightness_offset = 50;        % added to the frame so a dim scene stays visible
checkpoint_every = 200;        % samples between automatic saves

% 1 tracks only samples that are still NaN (a new list, or a resumed one).
% A value above 1 re-clicks from that sample through the end.
animal_start_sample = 1;
partner_start_sample = 1;

% true: section "Build the sample list" reloads resume_file instead of
% making a new frame list. Samples still NaN are the ones left to click.
resume_from_file = false;
resume_file = 'traking_structure PD2.mat';

animal_save_name = 'traking_structure PD2';
partner_save_name = 'traking_structure_P';
final_save_name = 'traking_structure_C2_P';
output_folder = pwd;

%% Build the sample list
% Samples are spaced by frame_rate / samples_per_second. With the defaults
% that step is 5 video frames (3 clicks per second at 15 fps).

if exist(video_file, 'file') ~= 2
    error('manual_traking:MissingVideo', 'Video file not found: %s', video_file);
end
if numel(track_window_s) ~= 2 || track_window_s(2) <= track_window_s(1)
    error('manual_traking:Window', 'track_window_s must be [start_s end_s] with end > start.');
end
if frame_rate <= 0 || samples_per_second <= 0
    error('manual_traking:Rate', 'frame_rate and samples_per_second must be positive.');
end

vidObj = VideoReader(video_file);
if abs(double(vidObj.FrameRate) - frame_rate) > 0.05
    warning('manual_traking:FrameRate', ...
        'Video frame rate is %.3f Hz; frame_rate is %.3f Hz.', ...
        vidObj.FrameRate, frame_rate);
end

if resume_from_file
    loaded = load(resume_file, 'traking_structure');
    traking_structure = loaded.traking_structure;
    frames2stract = traking_structure.frames2stract(:).';
    traking_values_animal = traking_structure.animal_pos;
    n_samples = numel(frames2stract);
    if size(traking_values_animal, 1) ~= n_samples
        error('manual_traking:LengthMismatch', ...
            'animal_pos has %d rows but frames2stract has %d.', ...
            size(traking_values_animal, 1), n_samples);
    end
    fprintf('Resumed %d samples from %s.\n', n_samples, resume_file);
else
    duration_s = track_window_s(2) - track_window_s(1);
    n_samples = ceil(duration_s * samples_per_second);
    start_frame = floor(track_window_s(1) * frame_rate);
    frame_step = frame_rate / samples_per_second;
    frames2stract = round(start_frame + (0:n_samples - 1) * frame_step);

    traking_values_animal = nan(n_samples, 2);
    traking_structure = struct();
    traking_structure.frames2stract = frames2stract;
    traking_structure.animal_pos = traking_values_animal;
end

n_video_frames = video_frame_count(vidObj);
if ~isempty(n_video_frames) && n_video_frames > 0
    if any(frames2stract < 1) || any(frames2stract > n_video_frames)
        error('manual_traking:FrameRange', ...
            ['Requested frames %d-%d, but the video has %d frames. ', ...
             'Adjust track_window_s or frame_rate.'], ...
            min(frames2stract), max(frames2stract), n_video_frames);
    end
end

traking_time = frames2stract / frame_rate;
traking_structure.frames2stract = frames2stract;
traking_structure.time = traking_time;
fprintf('Tracking %d samples, video frames %d-%d (%.1f-%.1f s).\n', ...
    n_samples, frames2stract(1), frames2stract(end), traking_time(1), traking_time(end));

%% Track the focal animal
% One left-click per sample, on the animal. A green x marks the stored
% point for this sample (empty until the first click). The red dot is the
% click just accepted. Every checkpoint_every clicks the figure is
% reopened so it does not keep every previous frame in memory.

samples_to_click = samples_for_clicking(traking_values_animal, animal_start_sample);
if isempty(samples_to_click)
    fprintf(['Focal-animal track has no samples to click. ', ...
        'Set animal_start_sample above 1 to replace from that sample.\n']);
else
    fprintf('Focal animal: %d samples, starting at %d.\n', ...
        numel(samples_to_click), samples_to_click(1));
    vidObj = VideoReader(video_file);
    fig = open_tracking_figure([0.25 0.25 0.5 0.75]);
    animal_clock = tic;
    n_clicked = 0;

    for sample_index = samples_to_click
        n_clicked = n_clicked + 1;
        show_video_frame(fig, vidObj, frames2stract(sample_index), brightness_offset, ...
            sample_title(frames2stract(sample_index), frame_rate, sample_index, n_samples));
        plot(traking_values_animal(sample_index, 1), traking_values_animal(sample_index, 2), ...
            'xg', 'MarkerSize', 10);

        new_traking_point = ginput(1);
        if isempty(new_traking_point)
            fprintf('Stopped animal tracking at sample %d.\n', sample_index);
            break;
        end

        traking_values_animal(sample_index, :) = new_traking_point;
        plot(new_traking_point(1), new_traking_point(2), '.r', 'MarkerSize', 10);
        traking_structure = store_animal_track(traking_structure, frames2stract, traking_values_animal, traking_time);

        if checkpoint_every > 0 && mod(n_clicked, checkpoint_every) == 0
            save_tracking_checkpoint(traking_structure, output_folder, animal_save_name, 'traking_structure_');
            close(fig);
            fig = open_tracking_figure([0.25 0.25 0.5 0.75]);
        end
    end

    traking_structure = store_animal_track(traking_structure, frames2stract, traking_values_animal, traking_time);
    save_tracking_checkpoint(traking_structure, output_folder, animal_save_name, 'traking_structure_');
    if isgraphics(fig)
        close(fig);
    end
    toc(animal_clock);
end

%% Track the partner
% Same frames as the focal animal. A green x is the focal animal. A green
% dot is the partner: the stored click when one exists, otherwise the
% focal-animal point as a placeholder. The red dot is the click just
% accepted. Unclicked samples stay NaN in partner_pos so a checkpoint is
% not filled with a copy of the focal animal.
%
% With partner_start_sample at 1, only NaN samples are clicked, so a
% checkpoint resumes without overwriting finished frames. Set it above 1
% to re-click from that sample through the end. Clear
% traking_structure.partner_pos before running this section to discard a
% previous partner track and start over.

n_samples = numel(frames2stract);
if size(traking_values_animal, 1) ~= n_samples || size(traking_structure.animal_pos, 1) ~= n_samples
    error('manual_traking:LengthMismatch', ...
        'Focal-animal clicks and frames2stract must have the same length before partner tracking.');
end

traking_values_partner = nan(n_samples, 2);
if isfield(traking_structure, 'partner_pos') && size(traking_structure.partner_pos, 1) == n_samples
    traking_values_partner = traking_structure.partner_pos;
elseif isfield(traking_structure, 'partner_pos')
    error('manual_traking:LengthMismatch', ...
        'partner_pos has %d rows but this session has %d samples.', ...
        size(traking_structure.partner_pos, 1), n_samples);
end

samples_to_click = samples_for_clicking(traking_values_partner, partner_start_sample);
if isempty(samples_to_click)
    fprintf(['Partner track has no samples to click. ', ...
        'Set partner_start_sample above 1 to replace from that sample.\n']);
else
    fprintf('Partner: %d samples, starting at %d.\n', ...
        numel(samples_to_click), samples_to_click(1));
    vidObj = VideoReader(video_file);
    fig = open_tracking_figure([0.25 0.25 0.5 0.75]);
    partner_clock = tic;
    n_clicked = 0;

    for sample_index = samples_to_click
        n_clicked = n_clicked + 1;
        show_video_frame(fig, vidObj, frames2stract(sample_index), brightness_offset, ...
            sample_title(frames2stract(sample_index), frame_rate, sample_index, n_samples));
        plot(traking_structure.animal_pos(sample_index, 1), traking_structure.animal_pos(sample_index, 2), ...
            'xg', 'MarkerSize', 10);
        partner_marker = traking_values_partner(sample_index, :);
        if any(~isfinite(partner_marker))
            partner_marker = traking_structure.animal_pos(sample_index, :);
        end
        plot(partner_marker(1), partner_marker(2), 'g.', 'MarkerSize', 20);

        new_traking_point = ginput(1);
        if isempty(new_traking_point)
            fprintf('Stopped partner tracking at sample %d.\n', sample_index);
            break;
        end

        traking_values_partner(sample_index, :) = new_traking_point;
        traking_structure.partner_pos = traking_values_partner;
        plot(new_traking_point(1), new_traking_point(2), '.r', 'MarkerSize', 10);

        if checkpoint_every > 0 && mod(n_clicked, checkpoint_every) == 0
            save_tracking_checkpoint(traking_structure, output_folder, partner_save_name, 'traking_structure_P_');
            close(fig);
            fig = open_tracking_figure([0.25 0.25 0.5 0.75]);
        end
    end

    traking_structure.partner_pos = traking_values_partner;
    save_tracking_checkpoint(traking_structure, output_folder, partner_save_name, 'traking_structure_P_');
    if isgraphics(fig)
        close(fig);
    end
    toc(partner_clock);
end

%% Correct selected focal-animal samples
% Optional. Replace the example list below before running this section.
% Indices are positions in frames2stract (sample 1 is the first click),
% not raw video frame numbers.
%
% A click with x greater than reject_click_above_x is discarded and the
% stored focal-animal point is left as it was. Use that right-hand region
% to skip a sample. Raise the threshold if the animal can stand to the
% right of column 721.

frames2correct = 3:20;
frames2correct = frames2correct(:).';
reject_click_above_x = 721;

if any(frames2correct < 1) || any(frames2correct > numel(frames2stract))
    error('manual_traking:CorrectIndex', ...
        'frames2correct must be sample indices between 1 and %d.', numel(frames2stract));
end

vidObj = VideoReader(video_file);
fig = open_tracking_figure([0 0 1 1]);

for sample_index = frames2correct
    show_video_frame(fig, vidObj, frames2stract(sample_index), brightness_offset, ...
        sample_title(frames2stract(sample_index), frame_rate, sample_index, numel(frames2stract)));
    plot(traking_structure.animal_pos(sample_index, 1), traking_structure.animal_pos(sample_index, 2), ...
        'rx', 'MarkerSize', 20);
    if isfield(traking_structure, 'partner_pos')
        plot(traking_structure.partner_pos(sample_index, 1), traking_structure.partner_pos(sample_index, 2), ...
            'g.', 'MarkerSize', 20);
    end

    new_traking_point = ginput(1);
    if isempty(new_traking_point)
        fprintf('Stopped corrections at sample %d.\n', sample_index);
        break;
    end

    % Write this sample immediately so a stopped section does not blank
    % the samples that were not reached.
    if new_traking_point(1) <= reject_click_above_x
        traking_structure.animal_pos(sample_index, :) = new_traking_point;
    end
    traking_values_animal(sample_index, :) = traking_structure.animal_pos(sample_index, :);
    plot(new_traking_point(1), new_traking_point(2), '.r', 'MarkerSize', 16);
end

if isgraphics(fig)
    close(fig);
end
save_tracking_file(traking_structure, fullfile(output_folder, partner_save_name));

%% Draw the arena borders
% Optional. On the frame shown, draw three polygons and double-click each
% one to finish it, in this order: outer border, inner border, inner box.
% Vertices are stored under the field names the rest of the repository
% already reads (outher_border keeps that spelling).

border_brightness_offset = 25;
border_frame = round(max(1, max(traking_structure.frames2stract) - 100));
vidObj = VideoReader(video_file);
n_video_frames = video_frame_count(vidObj);
if ~isempty(n_video_frames) && n_video_frames > 0
    border_frame = min(max(border_frame, 1), n_video_frames);
end

figure;
current_image = read(vidObj, border_frame);
imagesc(current_image + border_brightness_offset);
traking_structure.outher_border = polygon_vertices('Outer border: click vertices, then double-click to finish');
traking_structure.inner_border = polygon_vertices('Inner border: click vertices, then double-click to finish');
traking_structure.inner_box = polygon_vertices('Inner box: click vertices, then double-click to finish');
save_tracking_file(traking_structure, fullfile(output_folder, partner_save_name));

%% Preview smoothed tracks
% Optional check before the final save. Clicks are sorted, duplicate frame
% indices are dropped, and the remaining points are interpolated onto every
% video frame from the first sample to the last. The loess window is in
% those interpolated frames (3 frames is a light smooth). The preview does
% not write back into traking_structure; the saved clicks stay unsmoothed.
% Both animal_pos and partner_pos are required. NaN clicks are left out of
% the interpolation, so a gap is filled across rather than drawn as a break.
%
% Ctrl+C stops playback. preview_start_frame is an index into all_frames
% (not a video frame number) and is updated each frame so running this
% section again continues from the last one shown.

smooth_window = 3;
trail_length = 10;
preview_start_frame = 60;
preview_pause_s = 0.01;

[sampled_frames, animal_xy, partner_xy] = unique_sorted_tracks(traking_structure);
all_frames = sampled_frames(1):sampled_frames(end);
if preview_start_frame > numel(all_frames)
    preview_start_frame = numel(all_frames);
end
preview_start_frame = max(1, preview_start_frame);
animal_dense = interpolate_track(sampled_frames, animal_xy, all_frames);
partner_dense = interpolate_track(sampled_frames, partner_xy, all_frames);
animal_smooth = smooth_track(animal_dense, smooth_window);
partner_smooth = smooth_track(partner_dense, smooth_window);

vidObj = VideoReader(video_file);
figure;
for frame_index = preview_start_frame:numel(all_frames)
    preview_start_frame = frame_index;
    trail_start = max(1, frame_index - trail_length + 1);

    hold off;
    current_image = read(vidObj, round(all_frames(frame_index)));
    imagesc(current_image);
    hold on;
    plot_closed_polygon(traking_structure, 'outher_border', 'y-');
    plot_closed_polygon(traking_structure, 'inner_border', 'w-');
    plot_closed_polygon(traking_structure, 'inner_box', 'w:');
    plot(animal_smooth(trail_start:frame_index, 1), animal_smooth(trail_start:frame_index, 2), 'r');
    plot(animal_smooth(frame_index, 1), animal_smooth(frame_index, 2), '.r', 'MarkerSize', 10);
    plot(partner_smooth(trail_start:frame_index, 1), partner_smooth(trail_start:frame_index, 2), 'g');
    plot(partner_smooth(frame_index, 1), partner_smooth(frame_index, 2), '.g', 'MarkerSize', 10);
    title(num2str([all_frames(frame_index), all_frames(end)]), 'FontSize', 14);
    pause(preview_pause_s);
end

%% Save the final tracking struct
% Sort by frame, drop repeated frame indices (the later click is kept),
% and store time in seconds using frame_rate. Arena polygons drawn above
% are left on the struct. This file is the one to copy into Traking backups.

[sampled_frames, animal_xy, partner_xy] = unique_sorted_tracks(traking_structure);
frames2stract = sampled_frames;
traking_values_animal = animal_xy;
traking_values_partner = partner_xy;
n_samples = numel(frames2stract);
traking_time = frames2stract / frame_rate;
traking_structure.frames2stract = frames2stract;
traking_structure.animal_pos = animal_xy;
traking_structure.partner_pos = partner_xy;
traking_structure.time = traking_time;

save_tracking_file(traking_structure, fullfile(output_folder, final_save_name));

%% Local functions

function n_frames = video_frame_count(vidObj)
    n_frames = [];
    try
        n_frames = double(vidObj.NumFrames);
    catch
        n_frames = [];
    end
    if isempty(n_frames) || n_frames <= 0
        n_frames = floor(double(vidObj.Duration) * double(vidObj.FrameRate));
    end
end

function fig = open_tracking_figure(position)
    fig = figure('Units', 'normalized', 'OuterPosition', position);
end

function show_video_frame(fig, vidObj, frame_number, brightness_offset, title_text)
    figure(fig);
    hold off;
    imagesc(read(vidObj, round(frame_number)) + brightness_offset);
    hold on;
    title(title_text, 'Interpreter', 'none');
end

function sample_index = samples_for_clicking(xy, start_sample)
    % start_sample <= 1 keeps finished rows and returns only NaN rows.
    % start_sample > 1 re-clicks from that row through the last row.
    n_rows = size(xy, 1);
    if start_sample <= 1
        sample_index = find(any(~isfinite(xy), 2)).';
    elseif start_sample > n_rows
        sample_index = [];
    else
        sample_index = start_sample:n_rows;
    end
end

function vertices = polygon_vertices(title_text)
    title(title_text, 'Interpreter', 'none');
    polygon = drawpolygon;
    vertices = polygon.Position;
    if isempty(vertices)
        error('manual_traking:Border', 'No polygon was drawn.');
    end
end

function title_text = sample_title(frame_number, frame_rate, sample_index, n_samples)
    title_text = sprintf('Frame %d | t = %.2f s | sample %d / %d (%.1f%%)', ...
        frame_number, frame_number / frame_rate, sample_index, n_samples, ...
        100 * sample_index / n_samples);
end

function traking_structure = store_animal_track(traking_structure, frames2stract, animal_pos, track_time)
    traking_structure.frames2stract = frames2stract;
    traking_structure.animal_pos = animal_pos;
    traking_structure.time = track_time;
end

function save_tracking_checkpoint(traking_structure, output_folder, stable_name, timestamp_prefix)
    save_tracking_file(traking_structure, fullfile(output_folder, stable_name));
    stamped_name = [timestamp_prefix, datestr(now, 'yyyy_mm_dd_HH_MM_SS_FFF')];
    save_tracking_file(traking_structure, fullfile(output_folder, stamped_name));
end

function save_tracking_file(traking_structure, file_name)
    [~, ~, extension] = fileparts(file_name);
    if ~strcmpi(extension, '.mat')
        file_name = [file_name, '.mat'];
    end
    save(file_name, 'traking_structure');
    fprintf('Saved %s\n', file_name);
end

function [frames, animal_pos, partner_pos] = unique_sorted_tracks(traking_structure)
    if ~isfield(traking_structure, 'animal_pos') || ~isfield(traking_structure, 'partner_pos')
        error('manual_traking:MissingTracks', ...
            'animal_pos and partner_pos are both required before this step.');
    end
    frames = traking_structure.frames2stract(:).';
    [frames, order] = sort(frames);
    animal_pos = traking_structure.animal_pos(order, :);
    partner_pos = traking_structure.partner_pos(order, :);
    if size(partner_pos, 1) ~= numel(frames) || size(animal_pos, 1) ~= numel(frames)
        error('manual_traking:LengthMismatch', ...
            'frames2stract, animal_pos, and partner_pos must have the same number of rows.');
    end

    % diff == 0 marks the first of two identical frame indices. Dropping
    % those rows keeps the later click.
    drop = find(diff(frames) == 0);
    frames(drop) = [];
    animal_pos(drop, :) = [];
    partner_pos(drop, :) = [];
end

function dense = interpolate_track(sample_frames, xy, all_frames)
    sample_frames = sample_frames(:);
    query_frames = all_frames(:);
    valid = all(isfinite(xy), 2);
    dense = nan(numel(query_frames), 2);
    if nnz(valid) < 2
        return;
    end
    method = 'linear';
    if nnz(valid) >= 4
        method = 'cubic';
    end
    dense(:, 1) = interp1(sample_frames(valid), xy(valid, 1), query_frames, method);
    dense(:, 2) = interp1(sample_frames(valid), xy(valid, 2), query_frames, method);
end

function smoothed = smooth_track(xy, window)
    smoothed = xy;
    finite = all(isfinite(xy), 2);
    if ~any(finite) || window < 2
        return;
    end
    first = find(finite, 1, 'first');
    last = find(finite, 1, 'last');
    if (last - first + 1) < window
        return;
    end
    % Window is the argument after the method. Passing it in the dimension
    % slot (smoothdata(x, 3, 'loess')) does not smooth a column of points.
    % Only the finite span is smoothed; interpolation already left NaNs
    % outside the first and last valid click.
    smoothed(first:last, 1) = smoothdata(xy(first:last, 1), 'loess', window);
    smoothed(first:last, 2) = smoothdata(xy(first:last, 2), 'loess', window);
end

function plot_closed_polygon(traking_structure, field_name, line_style)
    if ~isfield(traking_structure, field_name) || isempty(traking_structure.(field_name))
        return;
    end
    vertices = traking_structure.(field_name);
    plot([vertices(:, 1); vertices(1, 1)], [vertices(:, 2); vertices(1, 2)], line_style);
end
