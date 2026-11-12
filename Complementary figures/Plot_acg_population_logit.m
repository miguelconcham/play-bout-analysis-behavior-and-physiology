% Plot_acg_population_logit
% All neurons together, on a logit scale, from a saved ACG shuffle.
%
% Overlap is the standard method. Section 1 loads
% locking_period_acg_shuffle_n500.mat and builds the matrices. Section 2
% only draws. After section 1 has run, edit section 2 and run that
% section alone. The image matrices stay in the workspace.
%
% imagesc_logit.sum.all is neurons by lag for every firing run.
% imagesc_logit.sum.power_high is the high 25% of delta power.
% imagesc_logit.sum.power_low is the low 25% of delta power.
% imagesc_logit.ratio.* is the same layout on the rate-ratio scale.
% percentile.* is the 0-100 percentile in the same row order.
% sort_score.* is the mean percentile from 200 to 400 ms used to sort.
% Row 1 is the highest score and is the top row of imagesc.

%% 1. Load the checkpoint and build the imagesc matrices
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_dir);
addpath(this_dir);
addpath(repo_root);
if exist('percentile_curve_stack', 'file') ~= 2
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

acg_method = 'overlap';
n_shuffles_default = 500;
% Mean percentile inside this lag window sets the row order. 0.2 s is 200 ms.
sort_window_s = [0.200 0.400];

data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\' ...
    'NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data'];
parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
checkpoint = fullfile(parsing_folder, shuffle_checkpoint_filename(n_shuffles_default, acg_method));
if exist(checkpoint, 'file') ~= 2
    error('Plot_acg_population_logit:checkpoint', ...
        'No shuffle checkpoint at %s.', checkpoint);
end

stored = load(checkpoint, 'shuffle_records', 'parameters');
if ~isfield(stored, 'shuffle_records') || ~isfield(stored, 'parameters') ...
        || ~isfield(stored.parameters, 'n_shuffles')
    error('Plot_acg_population_logit:checkpoint', ...
        '%s has no shuffle records.', checkpoint);
end
records = stored.shuffle_records;
if ~iscell(records)
    records = num2cell(records(:));
end
records = records(:);
n_shuffles = stored.parameters.n_shuffles;

column_id = {'all', 'power_high', 'power_low'};
column_title = {'All neurons', 'All neurons, high 25% power', 'All neurons, low 25% power'};
split_name = {'all', 'power', 'power'};
tail_name = {'', 'high', 'low'};
scale_id = {'sum', 'ratio'};
scale_title = {'sum-to-1', 'rate ratio (Poisson = 1)'};
field_name = {'sum_pctl', 'ratio_pctl'};

imagesc_logit = struct();
percentile = struct();
sort_score = struct();
lag = [];
for s = 1:numel(scale_id)
    for c = 1:numel(column_id)
        [this_lag, pctl] = percentile_curve_stack( ...
            records, split_name{c}, tail_name{c}, field_name{s});
        if isempty(pctl)
            imagesc_logit.(scale_id{s}).(column_id{c}) = zeros(0, 0);
            percentile.(scale_id{s}).(column_id{c}) = zeros(0, 0);
            sort_score.(scale_id{s}).(column_id{c}) = zeros(0, 1);
            continue
        end
        if isempty(lag)
            lag = this_lag;
        elseif numel(this_lag) ~= numel(lag) || max(abs(this_lag - lag)) > 1e-9
            error('Plot_acg_population_logit:axis', ...
                'The %s %s curve is not on the same lag axis.', ...
                scale_title{s}, column_title{c});
        end
        [pctl, logit, score] = order_percentile_rows(pctl, lag, sort_window_s, n_shuffles);
        imagesc_logit.(scale_id{s}).(column_id{c}) = logit;
        percentile.(scale_id{s}).(column_id{c}) = pctl;
        sort_score.(scale_id{s}).(column_id{c}) = score;
    end
end

% Plain names for the sum-to-1 imagesc matrices.
imagesc_sum_all = imagesc_logit.sum.all;
imagesc_sum_power_high = imagesc_logit.sum.power_high;
imagesc_sum_power_low = imagesc_logit.sum.power_low;

prep_file = fullfile(parsing_folder, 'acg_population_logit_prep.mat');
save(prep_file, 'lag', 'imagesc_logit', 'percentile', 'sort_score', ...
    'imagesc_sum_all', 'imagesc_sum_power_high', 'imagesc_sum_power_low', ...
    'column_id', 'column_title', 'scale_id', 'scale_title', ...
    'n_shuffles', 'sort_window_s', 'checkpoint', 'prep_file');
disp(['Section 1 saved the matrices to ', prep_file])
disp('imagesc_sum_all is neurons by lag. Run section 2 to draw.')

%% 2. Figure
% Re-run this section after editing it. Section 1 does not need to run
% again while these variables are still in the workspace.
this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
if ~exist('imagesc_logit', 'var')
    data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\' ...
        'NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data'];
    parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
    prep_file = fullfile(parsing_folder, 'acg_population_logit_prep.mat');
    load(prep_file)
end

% Symmetric color range for the logit image. 1 is about the 9th and 91st
% percentiles. 2 is about the 1st and 99th. Change this and re-run section 2.
logit_color_limit = 2;
high_cut = 95;
low_cut = 5;
n_rows = 4;
n_cols = numel(column_id);
n_colors = 256;
half = n_colors / 2;
logit_map = [ ...
    linspace(0.15, 0.96, half)', linspace(0.28, 0.96, half)', linspace(0.70, 0.96, half)'; ...
    linspace(0.96, 0.70, half)', linspace(0.96, 0.18, half)', linspace(0.96, 0.18, half)'];

for s = 1:numel(scale_id)
    figure('Name', ['ACG population logit, ', scale_title{s}], ...
        'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.06 0.08 0.88 0.84])
    for c = 1:n_cols
        % This is the imagesc matrix for the current column.
        this_logit = imagesc_logit.(scale_id{s}).(column_id{c});
        this_pctl = percentile.(scale_id{s}).(column_id{c});
        n = size(this_logit, 1);
        if n < 1
            for row = 1:n_rows
                subplot(n_rows, n_cols, (row - 1) * n_cols + c)
                if row == 1
                    title({column_title{c}, 'n = 0'})
                end
                axis off
            end
            continue
        end

        subplot(n_rows, n_cols, c)
        if n == 1
            imagesc(lag, [0.5 1.5], this_logit)
        else
            imagesc(lag, 1:n, this_logit)
        end
        axis tight
        clim([-logit_color_limit logit_color_limit])
        if n == 1
            ylim([0.5 1.5])
        end
        colormap(gca, logit_map)
        if n > 1
            set(gca, 'YTick', [1 n])
        end
        ylabel('Neuron')
        title({column_title{c}, sprintf('n = %d', n)})
        box off
        if c == n_cols
            colorbar
        end

        subplot(n_rows, n_cols, n_cols + c)
        hold on
        for i = 1:n
            plot(lag, this_logit(i, :), 'Color', [0.25 0.45 0.75 0.28], ...
                'LineWidth', 0.6, 'HandleVisibility', 'off')
        end
        plot(lag, mean(this_logit, 1, 'omitnan'), 'k', 'LineWidth', 2.4)
        plot(lag([1 end]), [0 0], ':', 'Color', [0.45 0.45 0.45], ...
            'HandleVisibility', 'off')
        ylabel('Logit')
        box off

        [high_percent, low_percent] = percentile_tail_fractions(this_pctl, high_cut, low_cut);
        subplot(n_rows, n_cols, 2 * n_cols + c)
        hold on
        plot(lag, high_percent, 'Color', [0.7 0.25 0.25], 'LineWidth', 1.8)
        plot(lag, low_percent, 'Color', [0.25 0.45 0.75], 'LineWidth', 1.8)
        ylim([0 100])
        ylabel('Neurons (%)')
        if c == 1
            legend({sprintf('Percentile >= %g', high_cut), ...
                sprintf('Percentile <= %g', low_cut)}, ...
                'Location', 'northeast', 'Box', 'off')
        end
        box off

        subplot(n_rows, n_cols, 3 * n_cols + c)
        hold on
        plot(lag([1 end]), [0 0], ':', 'Color', [0.45 0.45 0.45], ...
            'HandleVisibility', 'off')
        plot(lag, mean(this_logit, 1, 'omitnan'), 'k', 'LineWidth', 2.4)
        xlabel('Lag (s)')
        ylabel('Mean logit')
        box off
    end
    sgtitle({['All neurons, ', scale_title{s}], ...
        'Rows are sorted by mean percentile from 200 to 400 ms, then shown as logit.'})
end
