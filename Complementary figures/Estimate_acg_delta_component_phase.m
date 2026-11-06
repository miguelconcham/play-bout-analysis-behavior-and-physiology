%% Estimate_acg_delta_component_phase
% Loop NPX sessions, call GENERATE_ACG_DELTA_COMPONENT_PHASE, and save one
% file per animal.
%
% delta_component_lag is the autocorrelogram side peak, in seconds. It is
% the same window shaded on the locking-period autocorrelograms (0.2-0.5 s,
% 5-2 Hz). Leave it as a range: the peak lag is not refit per neuron.
% The phase band is the delta band used for the trough/peak labels.

%% Paths, animals, and filter
data_root      = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
npx_Raw_Data   = [data_root, '\NPX data\NPX raw data'];
saving_folder  = [data_root, '\Analysis results\acg delta component phase'];
run('\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');
if ~exist(saving_folder, 'dir')
    mkdir(saving_folder);
end

animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
animal_file_names = cellfun(@(x) ['B', x], strsplit([animal_list.name], 'B'), 'UniformOutput', false)';
animal_file_names(1) = [];
animal2exclude = {''};
animal_list(ismember(animal_file_names, animal2exclude)) = [];
animal_names = {};

delta_component_lag = [0.2 0.5]; % s. Change this to move the side-peak window.
sr           = 2500;
filter_order = 2000;
freq_range   = [1 5];   % delta band of the trough/peak labels
freq_tag     = sprintf('%g-%gHz', freq_range(1), freq_range(2));
lag_tag      = sprintf('lag_%g-%gms', 1000 * delta_component_lag(1), 1000 * delta_component_lag(2));

Hd_freq = designfilt('bandpassfir', ...
    'FilterOrder', filter_order, ...
    'CutoffFrequency1', freq_range(1), ...
    'CutoffFrequency2', freq_range(2), ...
    'SampleRate', sr, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');

%%
tic
for fn = 1:numel(animal_list)
    animal_path = [npx_Raw_Data, '\', animal_list(fn).name];
    animal_name = animal_list(fn).name;
    disp(['Animal ', num2str(fn), ' of ', num2str(numel(animal_list)), ': ', animal_name])

    session = GENERATE_ACG_DELTA_COMPONENT_PHASE(animal_path, Hd_freq, delta_component_lag);
    if isempty(session)
        disp(['  empty result, skipping save: ', animal_name])
        toc
        continue
    end

    out_file = [animal_name, '_', freq_tag, '_', lag_tag, '_acg_delta_component_phase.mat'];
    save([saving_folder, '\', out_file], 'session', '-v7.3');
    animal_names = [animal_names; {animal_name, out_file}]; %#ok<AGROW>
    disp(['  saved ', out_file])
    save([saving_folder, '\acg_delta_component_phase_animal_names.mat'], 'animal_names');
    clear session
    toc
end
