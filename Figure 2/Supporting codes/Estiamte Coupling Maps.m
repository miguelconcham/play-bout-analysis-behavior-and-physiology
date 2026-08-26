


npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];


animal_file_names =  cellfun(@(x) ['B', x],strsplit([animal_list.name], 'B'), 'UniformOutput',false)';
animal_file_names(1) = [];
% animal2exclude = {'B4D4 0826 Dual'};
animal2exclude = {''};
animal_list(ismember(animal_file_names,animal2exclude)) = [];
animal_names ={};
n_strctut = 1;

freq_range_1    = [.1 6];
freq_range_2    = [6 12];
sr              = 2500;
filter_order    = 2000;


psth_structure = [];
Hd_freq1 = designfilt('bandpassfir', ...
'FilterOrder', filter_order, ...
'CutoffFrequency1', freq_range_1(1), ...
'CutoffFrequency2', freq_range_1(2), ...
'SampleRate', sr, ...
'DesignMethod', 'window', ...
'Window', 'hamming');

Hd_freq2 = designfilt('bandpassfir', ...
'FilterOrder', filter_order, ...
'CutoffFrequency1', freq_range_2(1), ...
'CutoffFrequency2', freq_range_2(2), ...
'SampleRate', sr, ...
'DesignMethod', 'window', ...
'Window', 'hamming');


%%
for fn = 1:numel(animal_list)
    
    if fn==1
        psth_structure = GENERATE_FREQ_COUPLING_MAPS_STRUCT([npx_Raw_Data, '\', animal_list(fn).name],Hd_freq1,Hd_freq2 );
         
        n_strctut = n_strctut+numel(psth_structure);
        animal_names = [animal_names;[repmat(animal_list(fn).name,numel(psth_structure),1) num2cell(1:numel(psth_structure))']]
    else
        transt_psth = GENERATE_FREQ_COUPLING_MAPS_STRUCT([npx_Raw_Data, '\', animal_list(fn).name],Hd_freq1,Hd_freq2 );
      
        for sub_j=1:numel(transt_psth)
    
            psth_structure(n_strctut) = transt_psth(sub_j);
            n_strctut = n_strctut+1;
        end
        animal_names = [animal_names;[repmat({animal_list(fn).name},numel(transt_psth),1) num2cell(1:numel(transt_psth))' ]]

    end


end
%%

%% save if needed
disp('saving')
save([saving_folder,'\couplig_structure_V2.mat'],'psth_structure', '-v7.3');
save([saving_folder,'\animal_names_coupling_V2.mat'],'animal_names');
