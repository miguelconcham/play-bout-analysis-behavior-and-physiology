data_root                   = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\Codes repository\Data';
hmm_data_folder             = [data_root, '\HMM data\HMM raw data'];
labeled_data_folder         = [data_root, '\HMM data\locomotive behaviors'];
segmented_data_folder       = [data_root, '\Analysis results\locomotive behaviors 2 partners'];
not_labeled_data_folder     = [data_root, '\Analysis results\non labeled behavior 2 partners'];
behavior_folder             = [data_root, '\Behavior backups'];
dir_list = dir(fullfile(behavior_folder, '*.txt'));
spatial_property_names      = {  'Speed','AngleSpeed','AngleAcc','Acc','Wall2CenterPos'...
    'RelativeDistance','RelativeSpeed','RelativeAngleSpeed','RelativeAngleAcc', 'RelativeAcc'};
call_prop_list = {'PrincipalFrequencykHz', 'SlopekHzs', 'Sinuosity', 'DeltaFreqkHz', 'FrequencyStandardDeviationkHz'};


behaviors2check = {'Pin', 'Boxing', 'Evasion', 'Pounce_A','Pounce_B','CD','Escape','CC','CB','Pounce_Ai','Pounce_Bi','Rearing','Sniffing', 'Bite', 'Scratch', 'Grooming'}; %% here you decide what behavior to extract

%% obtain behavioral parameters (Create cell array)
get_non_labeled             = true;
extended_time               = [0 0];
all_behavior                = [];
all_not_labeled_behavior    = [];


for j=1:numel(dir_list)
    animal_code = strrep(dir_list(j).name, '.txt', '');

    [behavior_list,all_non_labeled_list, all_spatial_prop] =get_behavior_properties(animal_code, behaviors2check,call_prop_list, extended_time,get_non_labeled);
 
    all_behavior            = [all_behavior;behavior_list];
    all_not_labeled_behavior = [all_not_labeled_behavior;all_non_labeled_list];
end
size(all_behavior)

all_var_names = [all_spatial_prop,'NumCalls',call_prop_list];


%% save traking info for convolution decoder
folder2save = segmented_data_folder;
starting_n = 0;
for j=1:size(all_behavior,1)
    locomotive_matrix = all_behavior{j,1};
    locomotive_matrix = locomotive_matrix(~any(isnan(locomotive_matrix),2),:);
    if size(locomotive_matrix,2)==numel(all_var_names)
        writeNPY(locomotive_matrix, [folder2save,'\',num2str(j+starting_n), '.npy'])
    end
end
behavior_labels = all_behavior(:,2);
disp('real behavior ready' )
folder2save = not_labeled_data_folder;
starting_n = size(all_behavior,1);
for j=1:size(all_not_labeled_behavior,1)
    locomotive_matrix = all_not_labeled_behavior{j,1};
    locomotive_matrix = locomotive_matrix(~any(isnan(locomotive_matrix),2),:);
    if size(locomotive_matrix,2)==numel(all_var_names)
        writeNPY(locomotive_matrix, [folder2save,'\',num2str(j+starting_n), '.npy'])
    end
end
behavior_labels = [behavior_labels;repmat({''},size(all_not_labeled_behavior,1),1)];
disp('not labaled behavior ready' )