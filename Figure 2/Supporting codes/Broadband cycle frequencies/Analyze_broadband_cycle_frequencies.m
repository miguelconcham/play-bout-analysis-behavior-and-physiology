%% Analyze_broadband_cycle_frequencies
% Figures and the shareable text summary for low-pass cycle frequencies.
%
% Run Estimate_broadband_cycle_frequencies.m first.

%% Output
% Leave figure_folder empty to write figures inside the results folder.
figure_folder = '';
show_figures = false;

this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes\Broadband cycle frequencies\' ...
        'Analyze_broadband_cycle_frequencies.m'];
end
repo_root = fileparts(fileparts(fileparts(fileparts(this_file))));
run(fullfile(repo_root, 'add_repo_paths.m'));

options = broadband_cycle_options();
saving_folder = options.output_folder;
input_file = fullfile(saving_folder, options.output_filename);
if isempty(figure_folder)
    figure_folder = fullfile(saving_folder, 'figures');
end
if ~isfolder(figure_folder)
    mkdir(figure_folder);
end

fprintf('Loading %s\n', input_file);
load(input_file, 'recordings');
summary_file = fullfile(saving_folder, 'result_summary.txt');
write_broadband_cycle_summary_text(recordings, summary_file);
fprintf('Summary written to %s\n', summary_file);

for recording_n = 1:numel(recordings)
    recording = recordings(recording_n);
    figure_handle = new_figure(show_figures, [1000 700]);
    plot_recording(recording);
    save_vector_and_jpeg(figure_handle, figure_folder, sprintf( ...
        '%s_probe%g_channel%g', file_token(recording.animal_name), ...
        recording.lpag_probe_number, recording.lpag_channel));
    if ~show_figures
        close(figure_handle);
    end
end

figure_handle = new_figure(show_figures, [900 500]);
plot_equal_animal_histogram(recordings);
save_vector_and_jpeg(figure_handle, figure_folder, ...
    'equal_animal_cycle_frequency');
if ~show_figures
    close(figure_handle);
end
fprintf('Figures written to %s\n', figure_folder);

function plot_recording(recording)
subplot(2, 1, 1);
example = recording.example;
if isempty(example.time_s)
    text(0.5, 0.5, 'No accepted cycle', 'HorizontalAlignment', 'center');
    axis off
else
    plot(example.time_s, example.broad, 'Color', [0.15 0.25 0.45]);
    hold on
    plot(example.peak_time_s, example.broad(nearest_sample(example, ...
        example.peak_time_s)), 'k^', 'MarkerFaceColor', 'k');
    trough_y = example.broad(nearest_sample(example, example.trough_time_s));
    plot(example.trough_time_s, trough_y, 'v', ...
        'Color', [0.1 0.45 0.75], 'MarkerFaceColor', [0.1 0.45 0.75]);
    xlabel('Recording time (s)');
    ylabel('40 Hz low-pass');
    title(sprintf('%s probe %g channel %g, example %.2f Hz', ...
        recording.animal_name, recording.lpag_probe_number, ...
        recording.lpag_channel, example.frequency_hz), 'Interpreter', 'none');
    box off
end

subplot(2, 1, 2);
frequency = recording.frequency_hz(recording.accepted);
edges = recording.parameters.histogram_edges_hz;
histogram(frequency, edges, 'FaceColor', [0.2 0.35 0.55], 'EdgeColor', 'none');
xlabel('Cycle frequency (Hz)');
ylabel('Accepted cycles');
title(sprintf('%d accepted of %d candidate cycles', ...
    recording.n_accepted_cycles, recording.n_candidate_cycles));
box off
end

function index = nearest_sample(example, time_s)
index = zeros(size(time_s));
for time_n = 1:numel(time_s)
    [~, index(time_n)] = min(abs(example.time_s - time_s(time_n)));
end
end

function plot_equal_animal_histogram(recordings)
parameters = recordings(1).parameters;
animals = unique({recordings.animal_id}, 'stable');
edges = parameters.histogram_edges_hz;
n_bins = numel(edges) - 1;
density = nan(numel(animals), n_bins);
for animal_n = 1:numel(animals)
    frequency = zeros(0, 1);
    for recording_n = 1:numel(recordings)
        recording = recordings(recording_n);
        if strcmp(recording.animal_id, animals{animal_n})
            frequency = [frequency; recording.frequency_hz(recording.accepted)]; %#ok<AGROW>
        end
    end
    counts = zeros(1, n_bins);
    for bin_n = 1:n_bins
        counts(bin_n) = nnz(frequency >= edges(bin_n) & frequency < edges(bin_n + 1));
    end
    if sum(counts) > 0
        density(animal_n, :) = counts / sum(counts);
    end
end
centers = (edges(1:end - 1) + edges(2:end)) / 2;
mean_density = nan(1, size(density, 2));
for bin_n = 1:size(density, 2)
    column = density(:, bin_n);
    column = column(isfinite(column));
    if ~isempty(column)
        mean_density(bin_n) = mean(column);
    end
end
bar(centers, mean_density, 1, 'FaceColor', [0.2 0.35 0.55], 'EdgeColor', 'none');
xlabel('Cycle frequency (Hz)');
ylabel('Mean fraction of cycles per animal');
title(sprintf('Equal-animal cycle frequencies, %d animals', numel(animals)));
box off
end

function save_vector_and_jpeg(figure_handle, folder_name, base_name)
print(figure_handle, fullfile(folder_name, [base_name, '.svg']), '-dsvg');
print(figure_handle, fullfile(folder_name, [base_name, '.jpg']), ...
    '-djpeg', '-r200');
end

function figure_handle = new_figure(show_figures, size_pixels)
if show_figures
    visibility = 'on';
else
    visibility = 'off';
end
figure_handle = figure('Color', 'w', 'Visible', visibility, ...
    'Position', [80 80 size_pixels]);
set(figure_handle, 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');
end

function name = file_token(text)
name = regexprep(char(text), '[^\w\-]+', '_');
name = regexprep(name, '_+', '_');
name = regexprep(name, '^_|_$', '');
if isempty(name)
    name = 'recording';
end
end
