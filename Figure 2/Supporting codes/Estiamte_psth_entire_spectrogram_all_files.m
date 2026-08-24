

npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\psth power by frequency and behavior';
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
figure_3_new_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Codes\Figure codes\Figure 3 Updated';

animal_file_names =  cellfun(@(x) ['B', x],strsplit([animal_list.name], 'B'), 'UniformOutput',false)';
animal_file_names(1) = [];
% animal2exclude = {'B4D4 0826 Dual'};
animal2exclude = {''};
animal_list(ismember(animal_file_names,animal2exclude)) = [];
animal_names ={};
n_strctut = 1;
%  % for delta
% psth_structure = [];
% wind_length     = 1;
% wind_overlap    = .990;
% min_separation = .200;
% f               = .1:.5:6;
% freq_pow_range  = [1 5];

 % for low gamma
% psth_structure = [];
% wind_length     = .100;
% wind_overlap    = .090;
% min_separation = .200;
% f               = 20:50;
% freq_pow_range  = [25 45];



 % for  gamma
psth_structure = [];
wind_length     = .100;
wind_overlap    = .090;
min_separation = .200;
f               = 20:50;
freq_pow_range  = [35 90];


 % for  the entire spectrum
psth_structure = [];
wind_length     = 1;
wind_overlap    = .990;
min_separation = .200;
f               = 0.1:.1:30;
freq_pow_range  = [35 90];


%  % for betha
% psth_structure = [];
% wind_length     = .100;
% wind_overlap    = .090;
% min_separation = .200;
% f               = 10:.5:25;
% freq_pow_range  = [12 20];

%%
for fn = 1:numel(animal_list)
    if fn==1
        psth_structure = GENERATE_THETA_PSTH_SPECTROGRAM([npx_Raw_Data, '\', animal_list(fn).name],wind_length,wind_overlap,min_separation,f,freq_pow_range );
        n_strctut = n_strctut+numel(psth_structure);
        animal_names = [animal_names;[repmat(animal_list(fn).name,numel(psth_structure),1) num2cell(1:numel(psth_structure))']]
    else
        transt_psth = GENERATE_THETA_PSTH_SPECTROGRAM([npx_Raw_Data, '\', animal_list(fn).name],wind_length,wind_overlap,min_separation,f,freq_pow_range );

        for sub_j=1:numel(transt_psth)

            psth_structure(n_strctut) = transt_psth(sub_j);
            n_strctut = n_strctut+1;
        end
        animal_names = [animal_names;[repmat({animal_list(fn).name},numel(transt_psth),1) num2cell(1:numel(transt_psth))' ]]

    end


end
%% save if needed
disp('saving')
% save([saving_folder,'\psth_structure_delta_olny_pouncesdfsdf.mat'],'psth_structure', '-v7.3');
% save([saving_folder,'\animal_names_delta_olny_pouncesdfsdf.mat'],'animal_names');

save([saving_folder,'\psth_structure_delta_full_spectrogram.mat'],'psth_structure', '-v7.3');
save([saving_folder,'\animal_names_delta_full_spectrogram.mat'],'animal_names');
