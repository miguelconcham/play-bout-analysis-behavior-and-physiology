%% Analyze_psth_area_maps
% Streamlined analysis of area-resolved play-bout PSTH maps.
% Keeps only Plot 1 (LPAG vs VLPAG LME), Plot 2 (power percentile maps),
% and Plot 3 (LPAG 4-panel).

this_file = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Analyze_psth_area_maps.m';
repo_root = fileparts(fileparts(fileparts(this_file)));
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');

%% Load precomputed PSTH structures
disp('loading')
load([saving_folder,'\psth_structure_delta_map_corrected_and_updated.mat'],'psth_structure');
load([saving_folder,'\animal_names_delta_map_corrected__and_updated.mat'],'animal_names');
disp('loading ready')

%% Merge PSTH across sessions and compute area-average responses
smooth_wind = 20;
baseline_range = [-2 0];
animal_label = {'B1D1','B1S3','B2S2','B3D2', 'B4S2', 'B4D4'};
mean_response_per_area = [];
n_sample = 1;

for an = 1:numel(animal_label)
    bin_size = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
    psth_ranges = psth_structure(1).hist_range;
    time = psth_ranges(1):bin_size:psth_ranges(2)+bin_size;
    baseline_index = time < baseline_range(2) & time > baseline_range(1);

    all_psth_onset = [];
    session_index = [];
    hard_coded_x_coords = [8 40; 258 290; 508 540; 758 790];

    this_animal = animal_label{an};
    animal2merge = find(cell2mat(cellfun(@(x) contains(x, this_animal), animal_names(:,1), 'UniformOutput', false)))';
    sess_n = 1;

    for j = animal2merge
        if contains(animal_names{j}, animal_label)
            this_animal_playbouts = psth_structure(j).play_bouts_table;
            this_animal_lengths = diff(this_animal_playbouts');
            this_psth_onset = psth_structure(j).play_bout_onset;
            session_index = [session_index; ones(size(this_psth_onset, 2), 1) * sess_n];

            for ch = 1:384
                for trial = 1:size(this_psth_onset, 2)
                    this_psth_onset(ch, trial, :) = (this_psth_onset(ch, trial, :) - mean(this_psth_onset(ch, trial, baseline_index))) / std(this_psth_onset(ch, trial, baseline_index));
                    this_psth_onset(ch, trial, :) = movmean(this_psth_onset(ch, trial, :), smooth_wind);
                    this_psth_onset(ch, trial, time > this_animal_lengths(trial)) = NaN;
                end
            end
            all_psth_onset = cat(2, all_psth_onset, this_psth_onset);
            sess_n = sess_n + 1;
        end
    end

    channel_map = psth_structure(animal2merge(1)).channel_map(2:end, :);
    x_pos = channel_map(:, 1);
    y_pos = channel_map(:, 2);

    probe_list = 1;
    if all(ismember(x_pos, hard_coded_x_coords))
        probe_list = find(any(ismember(hard_coded_x_coords, x_pos), 2));
    end

    index2use = 1:383;
    areas = psth_structure(animal2merge(1)).areas_by_channel(2:end);
    if all(ismember(x_pos, hard_coded_x_coords))
        index2use = index2use(ismember(x_pos, hard_coded_x_coords(probe_list, :)));
        [~, order] = sort(y_pos(index2use));
        index2use = index2use(order);
    end
    areas = areas(index2use);

    mean_probe_response = squeeze(mean(all_psth_onset(index2use, :, :), 2, 'omitmissing'));
    area_list = unique(areas);
    session_list = unique(session_index);

    mean_response_this_probe = nan(numel(area_list), size(mean_probe_response, 2));
    mean_response_this_probe_this_session = nan(numel(area_list), numel(session_list), size(mean_probe_response, 2));

    for area_n = 1:numel(area_list)
        mean_this_area = mean(mean_probe_response(ismember(areas, area_list{area_n}), :));
        mean_response_this_probe(area_n, :) = mean_this_area;

        for sn = 1:numel(session_list)
            mean_probe_response_this_session = squeeze(mean(all_psth_onset(index2use, session_index == session_list(sn), :), 2, 'omitmissing'));
            mean_this_area_this_session = mean(mean_probe_response_this_session(ismember(areas, area_list{area_n}), :));
            mean_response_this_probe_this_session(area_n, sn, :) = mean_this_area_this_session;
        end
    end

    mean_response_per_area(n_sample).mean_response = mean_response_this_probe;
    mean_response_per_area(n_sample).area_list = area_list;
    mean_response_per_area(n_sample).animal = animal_label{an};
    mean_response_per_area(n_sample).probe_n = probe_list;
    mean_response_per_area(n_sample).mean_response_ps = mean_response_this_probe_this_session;
    n_sample = n_sample + 1;
end

%% Aggregate responses and power spectra by brain area
all_areas = [];
for j = 1:numel(mean_response_per_area)
    all_areas = [all_areas; mean_response_per_area(j).area_list];
end
area_list = unique(all_areas);

response_per_area = cell(numel(area_list), 1);
power_per_area_response = cell(numel(area_list), 1);
power_per_area_pctl = cell(numel(area_list), 1);
session_per_area = cell(numel(area_list), 1);

response_time = find(time >= 0 & time <= 6);
fs = 1 / mean(diff(time));
n_rand = 1000;

for j = 1:numel(mean_response_per_area)
    disp(j)
    this_session_responses = mean_response_per_area(j).mean_response_ps;
    this_session_areas = mean_response_per_area(j).area_list;

    for k = 1:numel(this_session_areas)
        if strcmp(this_session_areas{k}, 'isRT')
            this_session_areas{k} = 'isRt';
        end
        idx = ismember(area_list, this_session_areas(k));

        for sn = 1:size(this_session_responses, 2)
            response_per_area{idx} = [response_per_area{idx}; squeeze(this_session_responses(k, sn, :))'];
            [response_p, f] = pspectrum(squeeze(this_session_responses(k, sn, response_time)), fs);

            shufled_response = squeeze(this_session_responses(k, sn, :))';
            rand_pool = find(~isnan(shufled_response));
            shuffled_p = nan(n_rand, numel(response_p));

            for n = 1:n_rand
                idnex4nand = rand_pool(randsample(numel(rand_pool), numel(response_time), false));
                [baseline_p, f] = pspectrum(shufled_response(idnex4nand), fs);
                shuffled_p(n, :) = baseline_p;
            end

            pctls = zeros(1, numel(response_p));
            for t = 1:numel(response_p)
                pctls(t) = mean(shuffled_p(:, t) <= response_p(t)) * 100;
            end

            power_per_area_pctl{idx} = [power_per_area_pctl{idx}; pctls];
            power_per_area_response{idx} = [power_per_area_response{idx}; response_p'];
            session_per_area{idx} = [session_per_area{idx}; {[mean_response_per_area(j).animal, 'P', num2str(mean_response_per_area(j).probe_n(1))]}];
        end
    end
end

%% Plot 1 — compare area PSTHs (LPAG vs VLPAG) with LME
areas2plot = {'LPAG', 'VLPAG'};
y_lim = [-.1 .5];
x_lim = [-1 2];
areas_indexes = find(ismember(area_list, areas2plot))';
areas_this_pair = session_per_area(areas_indexes);

area1index = all(ismember(areas_this_pair{1}, areas_this_pair{2}), 2);
area2index = all(ismember(areas_this_pair{2}, areas_this_pair{1}), 2);

matrix2test = cell(numel(areas_indexes), 1);
areaindexes = {area1index, area2index};

figure
hold on
color_list = generateDistinctColors(numel(areas2plot));
n_col = 1;

for j = areas_indexes
    if size(response_per_area{j}, 1) > 1
        matrix2test{n_col} = response_per_area{j}(areaindexes{n_col}, :);
        plot(time, mean(response_per_area{j}), 'Color', color_list(n_col, :))
    elseif size(response_per_area{j}, 1) == 1
        plot(time, response_per_area{j})
    end
    title(area_list{j})
    xlim(x_lim)
    ylim(y_lim)
    n_col = n_col + 1;
end

repeated_measures = session_per_area{areas_indexes(1)}(area1index);
index2estimate = time >= x_lim(1) & time <= x_lim(2);
array2test = matrix2test{1} - matrix2test{2};
array2test = array2test(:, index2estimate);
[~, T] = size(array2test);
pvals = nan(T, 1);

for t = 1:T
    y = array2test(:, t);
    tbl = table(y, categorical(repeated_measures), 'VariableNames', {'y', 'Subject'});
    lme = fitlme(tbl, 'y ~ 1 + (1|Subject)');
    pvals(t) = lme.Coefficients.pValue(1);
end

h = pvals < 0.05;
legend(areas2plot)
plot([0 0], [-.1 .5], 'k', 'HandleVisibility', 'off')

beg_end = [find(diff([0; h; 0]) == 1) find(diff([0; h; 0]) == -1) - 1];
sub_time = time(index2estimate);
beg_end = sub_time(beg_end);
y_lim = ylim;

for j = 1:size(beg_end, 1)
    fill(beg_end(j, [1 2 2 1]), [-.1 -.1 0 0] + y_lim(2), 'k')
end

%% Plot 2 — power-spectrum percentile maps by area
x_lim = [0 8];
plot_area_list = {'SupCol','DLPAG','LPAG','VLPAG','DR'};

figure
for j = 1:numel(plot_area_list)
    subplot(3, 4, j)
    imagesc(f, 1:size(power_per_area_pctl{j}, 1), power_per_area_pctl{j})
    yticks(1:numel(session_per_area{j}))
    yticklabels(session_per_area{j})
    clim([0 100])
    title(plot_area_list{j})
    xlim(x_lim)
end

figure
for j = 1:numel(plot_area_list)
    subplot(3, 4, j)
    if size(power_per_area_pctl{j}, 1) > 1
        [~, ~, ci] = ttest(power_per_area_pctl{j});
        rand_pool = ~any(isnan(ci));
        plot(f, power_per_area_pctl{j}, ':k')
        hold on
        fill([f(rand_pool)' fliplr(f(rand_pool)')], [ci(1, rand_pool) fliplr(ci(2, rand_pool))], 'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
        plot(f, mean(power_per_area_pctl{j}, 'omitmissing'), 'k', 'LineWidth', 2)
    elseif size(power_per_area_pctl{j}, 1) == 1
        plot(f, power_per_area_pctl{j})
    end
    title(plot_area_list{j})
    xlim(x_lim)
    ylim([0 100])
end

merged_areas = {'LPAG', 'VLPAG'};
area_index = find(ismember(plot_area_list, merged_areas));
all_merged_pow = [];
for j = 1:numel(area_index)
    all_merged_pow = [all_merged_pow; power_per_area_response{area_index(j)}];
end

figure
[~, ~, ci] = ttest(all_merged_pow);
rand_pool = ~any(isnan(ci));
plot(f, all_merged_pow, ':k')
hold on
fill([f(rand_pool)' fliplr(f(rand_pool)')], [ci(1, rand_pool) fliplr(ci(2, rand_pool))], 'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
plot(f, mean(all_merged_pow), 'k')

%% Plot 3 — LPAG detail (PSTH + power spectrum, 4-panel)
figure
x_lim = [-1 2];
j = find(ismember(area_list, 'LPAG'));

subplot(1, 4, 1)
imagesc(time, 1:size(response_per_area{j}(1:end-1, :), 1), response_per_area{j}(1:end-1, :))
yticks(1:numel(session_per_area{j}(1:end-1)))
yticklabels(session_per_area{j}(1:end-1))
clim([-.5 .5])
title(area_list{j})
xlim(x_lim)

subplot(1, 4, 2)
if size(response_per_area{j}, 1) > 1
    [~, ~, ci] = ttest(response_per_area{j}(1:end-1, :));
    rand_pool = ~any(isnan(ci));
    plot(time, response_per_area{j}(1:end-1, :), ':k')
    hold on
    fill([time(rand_pool) fliplr(time(rand_pool))], [ci(1, rand_pool) fliplr(ci(2, rand_pool))], 'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
    plot(time, mean(response_per_area{j}(1:end-1, :)), 'k')
elseif size(response_per_area{j}, 1) == 1
    plot(time, response_per_area{j})
end
title(area_list{j})
xlim(x_lim)
ylim([-.1 .75])
plot([0 0], [-.1 .75], 'r')

x_lim = [0 8];
subplot(1, 4, 3)
matrix2test = power_per_area_pctl{j}(1:end-1, :);
imagesc(f, 1:size(matrix2test, 1), matrix2test)
yticks(1:numel(session_per_area{j}(1:end-1)))
yticklabels(session_per_area{j}(1:end-1))
clim([0 100])
title(area_list{j})
xlim(x_lim)

subplot(1, 4, 4)
if size(power_per_area_response{j}, 1) > 1
    [~, ~, ci] = ttest(matrix2test);
    rand_pool = ~any(isnan(ci));
    plot(f, matrix2test, ':k')
    hold on
    fill([f(rand_pool)' fliplr(f(rand_pool)')], [ci(1, rand_pool) fliplr(ci(2, rand_pool))], 'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
    plot(f, mean(matrix2test), 'k')
elseif size(power_per_area_response{j}, 1) == 1
    plot(f, matrix2test)
end
title(area_list{j})
xlim(x_lim)
