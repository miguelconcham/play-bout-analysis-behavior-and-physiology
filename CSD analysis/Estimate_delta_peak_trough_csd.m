%% Estimate_delta_peak_trough_csd
% Loop NPX sessions, cut raw LFP around delta peaks/troughs (LPAG center,
% play_session = Type 'Partners session' in the behavior txt), write npy.
%
% Output layout under Data/Analysis results/delta peak trough CSD/:
%   peak/<animal>/                 (NPX1)
%   peak/<animal>/probe<N>/        (NPX2 shank)
%   trough/<animal>/
%   trough/<animal>/probe<N>/
% Each leaf folder: lfp.npy, x_coords.npy, y_coords.npy, T0.npy, fs.npy, meta.mat

data_root     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
npx_Raw_Data  = fullfile(data_root, 'NPX data', 'NPX raw data');
saving_folder = fullfile(data_root, 'Analysis results', 'delta peak trough CSD');
peak_root     = fullfile(saving_folder, 'peak');
trough_root   = fullfile(saving_folder, 'trough');

repo_root = fileparts(mfilename('fullpath'));
repo_root = fileparts(repo_root);
run(fullfile(repo_root, 'add_repo_paths.m'));
addpath(fullfile(repo_root, 'CSD analysis'));

if ~exist(saving_folder, 'dir'); mkdir(saving_folder); end
if ~exist(peak_root, 'dir');     mkdir(peak_root);     end
if ~exist(trough_root, 'dir');   mkdir(trough_root);   end

%% Fixed parameters
sr           = 2500;
filter_order = 2000;
freq_delta   = [1 5];
t1           = -0.50;   % s  (±500 ms; can subsample later in Python)
t2           =  0.50;   % s

Hd_delta = designfilt('bandpassfir', ...
    'FilterOrder', filter_order, ...
    'CutoffFrequency1', freq_delta(1), ...
    'CutoffFrequency2', freq_delta(2), ...
    'SampleRate', sr, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');

%% Animals
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
animal_file_names = cellfun(@(x) ['B', x], strsplit([animal_list.name], 'B'), 'UniformOutput', false)';
animal_file_names(1) = [];
animal2exclude = {''};
animal_list(ismember(animal_file_names, animal2exclude)) = [];

summary = {};
tic
for fn = 1:numel(animal_list)
    animal_name = animal_list(fn).name;
    npx_dir = fullfile(npx_Raw_Data, animal_name);
    disp('============================================================')
    disp(['Estimate delta CSD LFP: ', animal_name])

    try
        results = GENERATE_DELTA_PEAK_TROUGH_CSD_LFP(npx_dir, Hd_delta, t1, t2);
    catch ME
        warning('Estimate_delta_peak_trough_csd:Fail', '%s failed: %s', animal_name, ME.message);
        summary(end + 1, :) = {animal_name, 'ERROR', ME.message, 0, 0}; %#ok<AGROW>
        continue
    end

    if isempty(results)
        summary(end + 1, :) = {animal_name, 'SKIP', 'no LPAG or no play_session', 0, 0}; %#ok<AGROW>
        continue
    end

    for j = 1:numel(results)
        r = results(j);
        if r.npx_type == 1
            peak_dir   = fullfile(peak_root, r.animal_code);
            trough_dir = fullfile(trough_root, r.animal_code);
        else
            peak_dir   = fullfile(peak_root, r.animal_code, r.shank_label);
            trough_dir = fullfile(trough_root, r.animal_code, r.shank_label);
        end

        if ~exist(peak_dir, 'dir');   mkdir(peak_dir);   end
        if ~exist(trough_dir, 'dir'); mkdir(trough_dir); end

        local_save_csd_folder(peak_dir, r.lfp_peak, r.x_coords, r.y_coords, r.T0, r);
        local_save_csd_folder(trough_dir, r.lfp_trough, r.x_coords, r.y_coords, r.T0, r);

        disp(sprintf('  saved peak/trough for %s (%s)  peak=%d  trough=%d  n_chan=%d', ...
            r.animal_code, r.shank_label, r.n_events_peak, r.n_events_trough, numel(r.channel_ids)));
        summary(end + 1, :) = {r.animal_code, r.shank_label, 'OK', ...
            r.n_events_peak, r.n_events_trough}; %#ok<AGROW>
    end
end
toc

save(fullfile(saving_folder, 'estimate_summary.mat'), 'summary', 't1', 't2', 'freq_delta', 'Hd_delta');
disp(['Done. Outputs in ', saving_folder, ' (folders: peak/ and trough/)'])

%% -------------------------------------------------------------------------
function local_save_csd_folder(folder, lfp_stack, x, y, T0, meta_src)
% lfp_stack: (n_events, n_chan, n_time) single/double
if isempty(lfp_stack) || size(lfp_stack, 1) == 0
    warning('Estimate_delta_peak_trough_csd:Empty', 'No events for %s — writing empty placeholder skipped', folder);
    return
end

write_csd_npy(fullfile(folder, 'lfp.npy'), single(lfp_stack));
write_csd_npy(fullfile(folder, 'x_coords.npy'), double(x(:)));
write_csd_npy(fullfile(folder, 'y_coords.npy'), double(y(:)));
write_csd_npy(fullfile(folder, 'T0.npy'), double(T0));
write_csd_npy(fullfile(folder, 'fs.npy'), double(meta_src.fs));

meta = struct();
meta.animal_code       = meta_src.animal_code;
meta.npx_type          = meta_src.npx_type;
meta.probe_num         = meta_src.probe_num;
meta.shank_label       = meta_src.shank_label;
meta.detection_channel = meta_src.detection_channel;
meta.channel_ids       = meta_src.channel_ids;
meta.T0                = meta_src.T0;
meta.fs                = meta_src.fs;
meta.t1                = meta_src.t1;
meta.t2                = meta_src.t2;
meta.n_events          = size(lfp_stack, 1);
meta.n_chan            = size(lfp_stack, 2);
meta.n_time            = size(lfp_stack, 3);
meta.baseline_correct  = false;
meta.rereferenced      = false;
meta.normalize         = 'none';  % reserved for later
save(fullfile(folder, 'meta.mat'), 'meta');
end
