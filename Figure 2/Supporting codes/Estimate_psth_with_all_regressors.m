


list_of_animals = {'B1D1 1013 Dual','B1S3 1008 Single','B1S3 1009 Single','B2S2 1110 Single2','B2S2 1111 Single2','B3D2 1130 Dual','B4S2 0825 Single'};
list_of_paertner = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2 3]};


saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Theta psth';

% f          = 4:.1:15;      % frequency range for spectrogram
% freq_range = [6 12]; 

% parameters for variance explaiend for theta and dleta (as on figure 2 ans
% supp figure)
% f          = {.1:.1:6, 5:.1:16};      % frequency range for spectrogram
% freq_range = {[1 5],[6 12]}; 
% wind_length = {1, .250};
% spect_bin_size = 0.005;


f          = {10:.5:35, 35:100};      % frequency range for spectrogram
freq_range = {[15 30],[40 90]}; 
wind_length = {.150, .08};
spect_bin_size = 0.005;
n_strctut = 1;

psth_structure = [];
animal_names = [];
%% DONT RUN UNLESS NEEDED
for fn = 1:numel(list_of_animals)

    for pt = list_of_paertner{fn}
        disp([list_of_animals{fn} , ' P', num2str(pt)])

        if fn==1 && pt==1
            psth_structure = GENERATE_THETA_CALL_LOC_REGRESSOR(list_of_animals{fn}, pt, f, freq_range,wind_length,spect_bin_size);
            n_strctut = n_strctut+numel(psth_structure);
            animal_names = [animal_names;[repmat({list_of_animals{fn}},numel(psth_structure),1), repmat({pt},numel(psth_structure),1),num2cell(1:numel(psth_structure))']]
        else

            transt_psth = GENERATE_THETA_CALL_LOC_REGRESSOR(list_of_animals{fn}, pt, f, freq_range,wind_length,spect_bin_size);

            for sub_j=1:numel(transt_psth)

                psth_structure(n_strctut) = transt_psth(sub_j);
                n_strctut = n_strctut+1;
            end
            animal_names = [animal_names;[repmat({list_of_animals{fn}},numel(transt_psth),1), repmat({pt},numel(transt_psth),1), num2cell(1:numel(transt_psth))' ]]

        end

    end


end
%% saving dont save unless you are very sure
% 
% save([saving_folder,'\psth_structure_all_regressors_delta_theta.mat'],'psth_structure', '-v7.3'); 
% save([saving_folder,'\animal_names_all_regressors_delta.mat'],'animal_names');
% save([saving_folder,'\psth_structure_all_regressors_beta_gamma.mat'],'psth_structure', '-v7.3'); 
% save([saving_folder,'\animal_names_all_regressors_beta_gamma.mat'],'animal_names');
