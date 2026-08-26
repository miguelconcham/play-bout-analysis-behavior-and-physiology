
two_animals_data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\TWO ANIMALS';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\Theta psth';
figure_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure Mutual informatin inputs';
animal_list = dir(two_animals_data);
animal_list(1:2) = [];
%% to change powerr ranges and behaviors you need to modify directly the GENERATE function
n_strctut = 1;
mi_structure = [];
animal_names = [];
%%
for fn = 1:numel(animal_list)

    if fn==1
        mi_structure = GENERATE_STRUCUTRE_animal_synch([two_animals_data, '\', animal_list(fn).name]);

        n_strctut = n_strctut+numel(mi_structure);
        animal_names = [animal_names;[repmat(animal_list(fn).name,numel(mi_structure),1) num2cell(1:numel(mi_structure))']]
    else
        transt_psth =GENERATE_STRUCUTRE_animal_synch([two_animals_data, '\', animal_list(fn).name]);
        ([two_animals_data, '\', animal_list(fn).name] );

        for sub_j=1:numel(transt_psth)

            mi_structure(n_strctut) = transt_psth(sub_j);
            n_strctut = n_strctut+1;
        end
        animal_names = [animal_names;[repmat(animal_list(fn).name,numel(transt_psth),1) num2cell(1:numel(transt_psth))']]


    end


end

% %% In older versions i didnt added behavior table in the strucutre, so i
% have to added later in a separate structure
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
save([saving_folder,'\mi_structure_20bins_delta_aggression.mat'],'mi_structure');
save([saving_folder,'\animal_names_mi_structure_20bins_delta_aggression.mat'],'animal_names');

