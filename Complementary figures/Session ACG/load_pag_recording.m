function recording = load_pag_recording(session_dir)
%LOAD_PAG_RECORDING  PAG spikes, cluster table, and LFP for one session folder.
%
% Spike times are seconds. LFP is channels by samples, and sample k is at
% time k/2500, the same clock GENERATE_SPIKE_TRAIN_PARSING uses. NPX type
% 2 is read from LFP_PAG.mat and type 1 from LFP_PAG.dat. cluster_info.ch
% stays 0-based. recording.ok is false when this folder is not a session,
% with the reason in recording.reason.

recording = empty_recording(session_dir);
times_file = fullfile(session_dir, 'spike_times_PAG.npy');
clusters_file = fullfile(session_dir, 'spike_clusters_PAG.npy');
info_file = fullfile(session_dir, 'cluster_info_PAG.tsv');
if exist(times_file, 'file') ~= 2 || exist(clusters_file, 'file') ~= 2
    recording.reason = 'PAG spike files are missing';
    return
end
if exist('readNPY', 'file') ~= 2
    error('session_acg:readNPY', ...
        ['readNPY is not on the MATLAB path. GENERATE_SPIKE_TRAIN_PARSING ' ...
        'uses the same helper.']);
end
recording.spikes = double(readNPY(times_file)) / 30000;
recording.clusters = double(readNPY(clusters_file));
recording.spikes = recording.spikes(:);
recording.clusters = recording.clusters(:);
if numel(recording.spikes) ~= numel(recording.clusters)
    recording.reason = 'spike times and spike clusters have different lengths';
    recording.spikes = [];
    recording.clusters = [];
    return
end
if exist(info_file, 'file') == 2
    recording.cluster_info = readtable(info_file, 'FileType', 'text', 'Delimiter', '\t');
    if ismember('group', recording.cluster_info.Properties.VariableNames)
        keep = ismember(recording.cluster_info.group, {'mua', 'good'});
        recording.cluster_info = recording.cluster_info(keep, :);
    end
end
[recording.lfp, recording.reason] = load_pag_lfp(session_dir);
if isempty(recording.lfp)
    recording.spikes = [];
    recording.clusters = [];
    return
end
recording.lfp_time = (1:size(recording.lfp, 2)) / recording.sample_rate;
recording.ok = true;
recording.reason = '';
end


function recording = empty_recording(session_dir)
[~, recording.animal] = fileparts(session_dir);
recording.ok = false;
recording.reason = '';
recording.spikes = [];
recording.clusters = [];
recording.cluster_info = table();
recording.lfp = [];
recording.lfp_time = [];
recording.sample_rate = 2500;
end


function [lfp, reason] = load_pag_lfp(session_dir)
lfp = [];
reason = '';
map_file = fullfile(session_dir, 'chann_map_PAG.mat');
mat_file = fullfile(session_dir, 'LFP_PAG.mat');
dat_file = fullfile(session_dir, 'LFP_PAG.dat');
if exist(map_file, 'file') ~= 2
    reason = 'chann_map_PAG.mat is missing';
    return
end
mapped = load(map_file, 'xcoords');
if ~isfield(mapped, 'xcoords')
    reason = 'chann_map_PAG.mat has no xcoords';
    return
end
hard_coded_x_coords_npx2 = [8 40; 258 290; 508 540; 758 790];
use_mat = any(ismember(mapped.xcoords, hard_coded_x_coords_npx2));
if use_mat
    if exist(mat_file, 'file') ~= 2
        reason = 'LFP_PAG.mat is missing';
        return
    end
    loaded = load(mat_file, 'LFP');
    if ~isfield(loaded, 'LFP')
        reason = 'LFP_PAG.mat has no LFP variable';
        return
    end
    lfp = double(loaded.LFP);
else
    if exist(dat_file, 'file') ~= 2
        reason = 'LFP_PAG.dat is missing';
        return
    end
    file_pointer = fopen(dat_file, 'r');
    if file_pointer < 0
        reason = 'LFP_PAG.dat could not be opened';
        return
    end
    cleanup = onCleanup(@() fclose(file_pointer)); %#ok<NASGU>
    raw = fread(file_pointer, 'int16');
    clear cleanup
    if mod(numel(raw), 384) ~= 0
        reason = 'LFP_PAG.dat length is not a multiple of 384 channels';
        lfp = [];
        return
    end
    lfp = double(reshape(raw, 384, numel(raw) / 384));
end
if isempty(lfp) || ndims(lfp) ~= 2 || size(lfp, 1) < 1 || size(lfp, 2) < 2
    reason = 'PAG LFP is empty';
    lfp = [];
end
end
