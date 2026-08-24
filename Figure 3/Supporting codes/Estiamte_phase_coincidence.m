


npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Cross_correlogram';
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];


animal_file_names =  cellfun(@(x) ['B', x],strsplit([animal_list.name], 'B'), 'UniformOutput',false)';
animal_file_names(1) = [];
% animal2exclude = {'B4D4 0826 Dual'};
animal2exclude = {''};
animal_list(ismember(animal_file_names,animal2exclude)) = [];
animal_names ={};
n_strctut = 1;

freq_range_1    = [1 5];
freq_range_2    = [6 12];
sr              = 2500;
filter_order    = 2000;



% Parameters for delta
Hd_freq = designfilt('bandpassfir', ...
'FilterOrder', filter_order, ...
'CutoffFrequency1', freq_range_1(1), ...
'CutoffFrequency2', freq_range_1(2), ...
'SampleRate', sr, ...
'DesignMethod', 'window', ...
'Window', 'hamming');
bin_size_freq = 0.01;

% Parameters for theta
% Hd_freq = designfilt('bandpassfir', ...
% 'FilterOrder', filter_order, ...
% 'CutoffFrequency1', freq_range_2(1), ...
% 'CutoffFrequency2', freq_range_2(2), ...
% 'SampleRate', sr, ...
% 'DesignMethod', 'window', ...
% 'Window', 'hamming');
% bin_size_freq = 0.001;
areas2analyse = {'DLPAG'	'DR'	'LPAG'	'SupCol'	'VLPAG'};

time_precision = 0.005;

bin_size = 0.002;
hist_range = [-1.5 1];

psth_edges = hist_range(1):bin_size:hist_range(2);
%%
tic
for fn = n_strctut:numel(animal_list)

    if fn==1
        transt_psth  = GENERATE_PHASE_COINCIDENCE_STRUCTURE([npx_Raw_Data, '\', animal_list(fn).name],Hd_freq,bin_size_freq ,hist_range, time_precision,areas2analyse);
        phase_struct = transt_psth;



        n_strctut = n_strctut+numel(phase_struct);
        animal_names = [animal_names;[repmat(animal_list(fn).name,numel(phase_struct),1) num2cell(1:numel(phase_struct))']]
        
    else
        transt_psth = GENERATE_PHASE_COINCIDENCE_STRUCTURE([npx_Raw_Data, '\', animal_list(fn).name],Hd_freq,bin_size_freq ,hist_range, time_precision,areas2analyse);
        

        for sub_j=1:numel(transt_psth)

            phase_struct(n_strctut) = transt_psth(sub_j);
            n_strctut = n_strctut+1;
        end
        animal_names = [animal_names;[repmat({animal_list(fn).name},numel(transt_psth),1) num2cell(1:numel(transt_psth))' ]]
        save([saving_folder,'\',animal_list(fn).name,'_FreqRange_',num2str([Hd_freq.CutoffFrequency1 Hd_freq.CutoffFrequency2])...
            ,'_phase_coincidence_structure_animal_names.mat'],'animal_names');
        save([saving_folder,'\',animal_list(fn).name,'_FreqRange_',num2str([Hd_freq.CutoffFrequency1 Hd_freq.CutoffFrequency2])...
            ,'_phase_coincidence_structure.mat'],'phase_struct', '-v7.3');

    end
    toc

end


%% save if needed
disp('saving')
save([saving_folder,'\ALL AANIMALS_FreqRange_',num2str([Hd_freq.CutoffFrequency1 Hd_freq.CutoffFrequency2])...
            ,'_phase_coincidence_structure.mat'],'phase_struct', '-v7.3');
save([saving_folder,'\ALL_ANIMALS_FreqRange_',num2str([Hd_freq.CutoffFrequency1 Hd_freq.CutoffFrequency2])...
            ,'_phase_coincidence_structure_animal_names.mat'],'animal_names');

%% load
data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
saving_folder = [data_root, '\Analysis results\phase locking data'];

load([saving_folder,'\theta_phase_couplig_structure_updated_with_non_playbouts.mat'],'phase_struct');
load([saving_folder,'\theta_phase_couplig_animal_names_updated_with_non_playbouts.mat'],'animal_names');
% phase_struct(6 )=[];
% animal_names(6,:)=[];
%%
