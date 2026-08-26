
two_animals_data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\TWO ANIMALS';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\Theta psth';
figure_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure Mutual informatin inputs';
animal_list = dir(two_animals_data);
animal_list(1:2) = [];
%%
n_strctut = 1;
cc_delta_structure = [];
animal_names = [];



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
bin_size_freq = 0.1;

%%
for fn = 1:numel(animal_list)

    if fn==1
        cc_delta_structure = GENERATE_STRUCUTRE_animal_synch_cross_correlogram([two_animals_data, '\', animal_list(fn).name],Hd_freq);

        n_strctut = n_strctut+numel(cc_delta_structure);
        animal_names = [animal_names;[repmat(animal_list(fn).name,numel(cc_delta_structure),1) num2cell(1:numel(cc_delta_structure))']]
    else
        transt_psth =GENERATE_STRUCUTRE_animal_synch_cross_correlogram([two_animals_data, '\', animal_list(fn).name],Hd_freq);
        ([two_animals_data, '\', animal_list(fn).name] );

        for sub_j=1:numel(transt_psth)

            cc_delta_structure(n_strctut) = transt_psth(sub_j);
            n_strctut = n_strctut+1;
        end
        animal_names = [animal_names;[repmat(animal_list(fn).name,numel(transt_psth),1) num2cell(1:numel(transt_psth))']]


    end


end

% %%
% 
% behavior_structure = [];
% animal_names_behavior = [];
% n_strctut = 1;
% %%
% for fn = 1:numel(animal_list)
% 
%     if fn==1
%         behavior_structure = GENERATE_STRUCUTRE_animal_synch_only_behavior([two_animals_data, '\', animal_list(fn).name]);
% 
%         n_strctut = n_strctut+numel(behavior_structure);
%         animal_names_behavior = [animal_names_behavior;[repmat(animal_list(fn).name,numel(behavior_structure),1) num2cell(1:numel(behavior_structure))']];
%     else
%         transt_psth = GENERATE_STRUCUTRE_animal_synch_only_behavior([two_animals_data, '\', animal_list(fn).name] );
% 
%         for sub_j=1:numel(transt_psth)
% 
%             behavior_structure(n_strctut) = transt_psth(sub_j);
%             n_strctut = n_strctut+1;
%         end
%         animal_names_behavior = [animal_names_behavior;[repmat(animal_list(fn).name,numel(transt_psth),1) num2cell(1:numel(transt_psth))']]
% 
% 
%     end
% 
% 
% end


%%

% disp('saving')
% save([saving_folder,'\behavior_structure_20bins_delta_no_play.mat'],'behavior_structure');
% save([saving_folder,'\animal_names_behavior_mi_structure_20bins_delta_no_play.mat'],'animal_names_behavior');


%% SAve structure
% 
% disp('saving')
save([saving_folder,'\between_animals_cross_corr_structure_20bins_delta_play.mat'],'cc_delta_structure');
save([saving_folder,'\between_animals_cross_corr_mi_structure_20bins_delta_play.mat'],'animal_names');


