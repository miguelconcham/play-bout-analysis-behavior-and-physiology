npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';

saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\Rate_Coinc_Phase';

animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];

animal_file_names = cellfun(@(x) ['B', x], ...
    strsplit([animal_list.name], 'B'), ...
    'UniformOutput', false)';

animal_file_names(1) = [];

animal2exclude = {''};
animal_list(ismember(animal_file_names, animal2exclude)) = [];

animal_names = {};
synch_structure = [];


freq_range_1    = [1 5];
sr              = 2500;
filter_order    = 2000;



Hd_freq = designfilt('bandpassfir', ...
'FilterOrder', filter_order, ...
'CutoffFrequency1', freq_range_1(1), ...
'CutoffFrequency2', freq_range_1(2), ...
'SampleRate', sr, ...
'DesignMethod', 'window', ...
'Window', 'hamming');


n_strctut = 1;

%% PARAMETERS
time_precision = 0.003;

bin_size   = 0.001;
hist_range = [-1.5 1];

areas2analyse = {'DLPAG','DR','LPAG','SupCol','VLPAG'};
behaviors4playbout = {'Pounce','CC','Boxing','Evasion','Pin','Escape','CB','CD'};

%% =========================
% LOOP OVER SESSIONS
%% =========================

tic

for fn = n_strctut:numel(animal_list)

    session_path = [npx_Raw_Data, '\', animal_list(fn).name];

    disp(['Processing: ', animal_list(fn).name])

    % -------------------------
    % RUN CORE FUNCTION
    % -------------------------
    tmp = GENERATE_RATE_COINCIDENCE_PHASE_STRUCT( ...
        session_path, ...
        Hd_freq,...
        bin_size, ...
        hist_range, ...
        time_precision, ...
        areas2analyse);

    % -------------------------
    % STACK STRUCTURES
    % -------------------------
    if isempty(synch_structure)

        synch_structure = tmp;

    else

        synch_structure(n_strctut) = tmp;

    end

    animal_names = [animal_names; ...
        {animal_list(fn).name, 1}];

    n_strctut = n_strctut + 1;

    toc

    % -------------------------
    % INTERMEDIATE SAVE
    % -------------------------
    d = datestr(datetime('now'));
    d = strrep(d, ':', '_');

    save([saving_folder,'\rate_coinc_phase_STRUCT_',d,'.mat'], ...
        'synch_structure','-v7.3');

    save([saving_folder,'\animal_names_',d,'.mat'], ...
        'animal_names');

    disp('saved')

end

%%

%% ============================================================
% PHASE + COINCIDENCE DYNAMICS SUMMARY (PAIR LEVEL)
% Uses your synch_structure fields:
%
% rate_n1, rate_n2
% coin_n1, coin_n2
% phase_n1, phase_n2
% plv_n1, plv_n2
%% ============================================================

pair_phase_summary = struct();

nSessions = numel(synch_structure);

for s = 1:nSessions

    fprintf('Session %d / %d\n', s, nSessions)

    % ------------------------------------------------------------
    % EXTRACT CORE VARIABLES
    % ------------------------------------------------------------
    phase_n1 = synch_structure(s).phase_n1;
    phase_n2 = synch_structure(s).phase_n2;

    coin_n1  = synch_structure(s).coin_n1;
    coin_n2  = synch_structure(s).coin_n2;

    bin_size   = synch_structure(s).bin_size;
    hist_range = synch_structure(s).hist_range;

    nPairs = size(phase_n1,1);
    nBins  = size(phase_n1,2);

    % ------------------------------------------------------------
    % DEFINE PRE / POST WINDOWS
    % assumes bins are centered on hist_range
    % ------------------------------------------------------------
    bin_edges = linspace(hist_range(1), hist_range(2), nBins);

    pre_idx  = bin_edges < 0;
    post_idx = bin_edges >= 0;

    % ------------------------------------------------------------
    % PREALLOCATE OUTPUT
    % ------------------------------------------------------------
    pre_MVL_n1   = nan(nPairs,1);
    post_MVL_n1  = nan(nPairs,1);

    pre_MVL_n2   = nan(nPairs,1);
    post_MVL_n2  = nan(nPairs,1);

    pre_ang_n1   = nan(nPairs,1);
    post_ang_n1  = nan(nPairs,1);

    pre_ang_n2   = nan(nPairs,1);
    post_ang_n2  = nan(nPairs,1);

    dMVL_n1      = nan(nPairs,1);
    dMVL_n2      = nan(nPairs,1);

    dang_n1      = nan(nPairs,1);
    dang_n2      = nan(nPairs,1);

    % ============================================================
    % LOOP OVER PAIRS
    % ============================================================
    for p = 1:nPairs

        %% --------------------------
        % PHASE N1
        %% --------------------------
        ph_pre  = phase_n1(p,pre_idx);
        ph_post = phase_n1(p,post_idx);

        ph_pre  = ph_pre(~isnan(ph_pre));
        ph_post = ph_post(~isnan(ph_post));

        if ~isempty(ph_pre)
            pre_MVL_n1(p)  = abs(mean(exp(1i*ph_pre)));
            pre_ang_n1(p)  = angle(mean(exp(1i*ph_pre)));
        end

        if ~isempty(ph_post)
            post_MVL_n1(p) = abs(mean(exp(1i*ph_post)));
            post_ang_n1(p) = angle(mean(exp(1i*ph_post)));
        end

        %% --------------------------
        % PHASE N2
        %% --------------------------
        ph_pre  = phase_n2(p,pre_idx);
        ph_post = phase_n2(p,post_idx);

        ph_pre  = ph_pre(~isnan(ph_pre));
        ph_post = ph_post(~isnan(ph_post));

        if ~isempty(ph_pre)
            pre_MVL_n2(p)  = abs(mean(exp(1i*ph_pre)));
            pre_ang_n2(p)  = angle(mean(exp(1i*ph_pre)));
        end

        if ~isempty(ph_post)
            post_MVL_n2(p) = abs(mean(exp(1i*ph_post)));
            post_ang_n2(p) = angle(mean(exp(1i*ph_post)));
        end

        %% --------------------------
        % DELTAS (CIRCULAR CORRECT)
        %% --------------------------
        if ~isnan(pre_MVL_n1(p)) && ~isnan(post_MVL_n1(p))
            dMVL_n1(p) = post_MVL_n1(p) - pre_MVL_n1(p);
        end

        if ~isnan(pre_MVL_n2(p)) && ~isnan(post_MVL_n2(p))
            dMVL_n2(p) = post_MVL_n2(p) - pre_MVL_n2(p);
        end

        if ~isnan(pre_ang_n1(p)) && ~isnan(post_ang_n1(p))
            dang_n1(p) = angle(exp(1i*(post_ang_n1(p) - pre_ang_n1(p))));
        end

        if ~isnan(pre_ang_n2(p)) && ~isnan(post_ang_n2(p))
            dang_n2(p) = angle(exp(1i*(post_ang_n2(p) - pre_ang_n2(p))));
        end

    end

    % ============================================================
    % STORE
    % ============================================================
    pair_phase_summary(s).pre_MVL_n1  = pre_MVL_n1;
    pair_phase_summary(s).post_MVL_n1 = post_MVL_n1;
    pair_phase_summary(s).pre_MVL_n2  = pre_MVL_n2;
    pair_phase_summary(s).post_MVL_n2 = post_MVL_n2;

    pair_phase_summary(s).delta_MVL_n1 = dMVL_n1;
    pair_phase_summary(s).delta_MVL_n2 = dMVL_n2;

    pair_phase_summary(s).pre_ang_n1  = pre_ang_n1;
    pair_phase_summary(s).post_ang_n1 = post_ang_n1;
    pair_phase_summary(s).pre_ang_n2  = pre_ang_n2;
    pair_phase_summary(s).post_ang_n2 = post_ang_n2;

    pair_phase_summary(s).delta_ang_n1 = dang_n1;
    pair_phase_summary(s).delta_ang_n2 = dang_n2;

    pair_phase_summary(s).synch_comb   = synch_structure(s).synch_comb;
    pair_phase_summary(s).cluster_info = synch_structure(s).cluster_info;

end
%%

%% ============================================================
% SHIFT-O-GRAM: MVL CHANGE BEFORE vs AFTER PLAY ONSET
%% ============================================================

all_dMVL_n1 = [];
all_dMVL_n2 = [];

for s = 1:numel(pair_phase_summary)

    all_dMVL_n1 = [all_dMVL_n1;
        pair_phase_summary(s).delta_MVL_n1(:)];

    all_dMVL_n2 = [all_dMVL_n2;
        pair_phase_summary(s).delta_MVL_n2(:)];

end

% remove NaNs
all_dMVL_n1 = all_dMVL_n1(~isnan(all_dMVL_n1));
all_dMVL_n2 = all_dMVL_n2(~isnan(all_dMVL_n2));

%% ============================================================
% PLOT
%% ============================================================

figure;

subplot(1,2,1)
histogram(all_dMVL_n1, 50);
hold on;
xline(0,'r','LineWidth',2);

xlabel('\Delta MVL (Neuron 1)');
ylabel('Count');
title('Phase Concentration Change (N1)');
grid on;

subplot(1,2,2)
histogram(all_dMVL_n2, 50);
hold on;
xline(0,'r','LineWidth',2);

xlabel('\Delta MVL (Neuron 2)');
ylabel('Count');
title('Phase Concentration Change (N2)');
grid on;