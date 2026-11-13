%% Estimate group cross-correlations for all sessions
% Six groups: peak, trough, and unlocked, each split into oscillatory and
% non-oscillatory. Neurons come from oscillatory_neuron_labels.csv. Lock
% type and area come from acg_population_logit_prep.mat, joined on session
% and neuron id. A neuron that is absent from that prep is looked up in
% delta_all_neurons_v2.mat. The studied areas are SC, VLPAG, DLPAG, LPAG,
% and DR.
%
% Inside one recording, spikes of a group are pooled and cross-correlated
% with each of the other groups at 1 ms. Partners-session intervals are
% the time base. Counts from those intervals are added, then divided by
% the reference count and the bin width. Analyze normalizes each session
% by its own target rate before the mean and the confidence interval.
% Run add_repo_paths.m before this script.

this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_dir));
run(fullfile(repo_root, 'add_repo_paths.m'));

data_root = group_ccg_data_root();
npx_raw_data = fullfile(data_root, 'NPX data', 'NPX raw data');
label_file = fullfile(data_root, 'Analysis results', 'acg population logit', ...
    'oscillatory_neuron_labels.csv');
prep_file = fullfile(data_root, 'Analysis results', 'acg population logit', ...
    'acg_population_logit_prep.mat');
delta_file = fullfile(data_root, 'Analysis results', 'phase locking data', ...
    'delta_all_neurons_v2.mat');
saving_folder = fullfile(data_root, 'Analysis results', 'group cross correlation');
if exist(saving_folder, 'dir') ~= 7
    mkdir(saving_folder);
end

generator_options.acg_lag = 1;
generator_options.acg_bin = 0.001;
generator_options.areas_to_keep = {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'};
generator_options.session_type = 'Partners session';
generator_options.min_reference_spikes = 50;
generator_options.alpha = 0.01;
generator_options.non_entrained_p = 0.1;
settings = group_ccg_settings(generator_options);
RECOMPUTE_EXISTING = false;

labels = read_oscillatory_labels(label_file);
identity = empty_identity();
if exist(prep_file, 'file') == 2
    stored = load(prep_file, 'neuron_type', 'neuron_session', 'neuron_area', 'neuron_cluster_id');
    required = {'neuron_type', 'neuron_session', 'neuron_area', 'neuron_cluster_id'};
    if any(~isfield(stored, required))
        error('group_ccg:prep', ...
            '%s does not contain the neuron identity saved by Plot_acg_population_logit.', prep_file);
    end
    identity = identity_from_population(stored.neuron_session, stored.neuron_area, ...
        stored.neuron_cluster_id, stored.neuron_type);
    clear stored
end
if exist(delta_file, 'file') == 2
    loaded = load(delta_file, 'all_neurons');
    if isfield(loaded, 'all_neurons')
        if istable(loaded.all_neurons) && ismember('area', loaded.all_neurons.Properties.VariableNames)
            loaded.all_neurons.area(ismember(loaded.all_neurons.area, {'isRT'})) = {'isRt'};
        end
        extra = identity_from_delta_table(loaded.all_neurons, settings.alpha, settings.non_entrained_p);
        identity = merge_group_identity(identity, extra);
    end
    clear loaded extra
end
if isempty(identity.neuron_id)
    error('group_ccg:identity', ...
        ['No lock type or area was found. Expected %s or %s.'], prep_file, delta_file);
end
members = join_group_membership(labels, identity, settings.areas_to_keep);
member_table = table(members.session, members.neuron_id, members.area, ...
    members.lock_group, members.oscillatory, members.group_id, ...
    'VariableNames', {'session', 'neuron_id', 'area', 'lock_group', 'oscillatory', 'group_id'});
writetable(member_table, fullfile(saving_folder, 'group_ccg_members.csv'));
fprintf(['Labels %d, kept %d, not classified %d, inconsistent %d, unmatched %d, ', ...
    'no lock %d, outside the studied areas %d.\n'], ...
    members.n_labels, members.n_kept, members.n_unclassified, members.n_inconsistent, ...
    members.n_unmatched, members.n_no_lock, members.n_outside_area);

session_names = unique(members.session, 'stable');
failed_sessions = table('Size', [0 2], ...
    'VariableTypes', {'string', 'string'}, ...
    'VariableNames', {'Animal', 'Message'});
skipped_sessions = strings(0, 1);
n_added = 0;

for fn = 1:numel(session_names)
    animal_name = session_names{fn};
    session_dir = fullfile(npx_raw_data, animal_name);
    session_file = fullfile(saving_folder, group_ccg_filename(animal_name));
    fprintf('\nSession %d/%d: %s\n', fn, numel(session_names), animal_name);
    if exist(session_dir, 'dir') ~= 7
        failed_sessions(end + 1, :) = {string(animal_name), "NPX folder was not found."}; %#ok<SAGROW>
        warning('Group CCG: no NPX folder for %s.', animal_name);
        continue
    end
    if ~RECOMPUTE_EXISTING && exist(session_file, 'file') == 2
        fprintf('Keeping existing file: %s\n', session_file);
        skipped_sessions(end + 1, 1) = string(animal_name); %#ok<SAGROW>
        continue
    end
    try
        group_ccg = GENERATE_GROUP_CROSS_CORRELATION(session_dir, members, generator_options);
        save(session_file, 'group_ccg', '-v7.3');
        n_added = n_added + 1;
        clear group_ccg
    catch exception
        warning('Group CCG failed for %s: %s', animal_name, exception.message);
        failed_sessions(end + 1, :) = ...
            {string(animal_name), string(getReport(exception, 'extended', 'hyperlinks', 'off'))}; %#ok<SAGROW>
    end
end

retry_file = fullfile(saving_folder, 'group_ccg_retry.mat');
save(retry_file, 'failed_sessions', 'skipped_sessions', 'settings', 'members', 'n_added');
fprintf('Kept %d existing sessions; added %d; %d failed.\nLog: %s\n', ...
    numel(skipped_sessions), n_added, height(failed_sessions), retry_file);


function identity = empty_identity()
identity.session = {};
identity.area = {};
identity.neuron_id = zeros(0, 1);
identity.lock_group = {};
end
