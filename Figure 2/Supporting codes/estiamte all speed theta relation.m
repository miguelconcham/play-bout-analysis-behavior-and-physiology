
list_of_animals = {'B1D1 1013 Dual','B1S3 1008 Single','B1S3 1009 Single','B2S2 1110 Single2','B2S2 1111 Single2','B3D2 1130 Dual','B4S2 0825 Single'};
list_of_paertner = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2 3]};


saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';


n_strctut = 1;

psth_structure = [];
animal_names = [];

% Spectrogram parameters
% wind_length     = 1;    % delta
% wind_overlap    = 0.99;   % delta

wind_length     = .250;    % for theta
wind_overlap    = 0.249;   % for theta

% f          = .1:.1:5;      % frequency range for spectrogram delta
f          = 5:.1:16;       % frequency range for spectrogram thea
% freq_range = [.1 5];       % freq band delta
freq_range = [6 12];       % freq band theta


%%
for fn = 1:numel(list_of_animals)

    for pt = list_of_paertner{fn}
        disp([list_of_animals{fn} , ' P', num2str(pt)])

        if fn==1 && pt==1
            psth_structure = SPEED_THETA_RELATION(list_of_animals{fn}, pt,freq_range,f,wind_length,wind_overlap );
            n_strctut = n_strctut+numel(psth_structure);
            animal_names = [animal_names;[repmat({list_of_animals{fn}},numel(psth_structure),1), repmat({pt},numel(psth_structure),1),num2cell(1:numel(psth_structure))']]
        else

            transt_psth = SPEED_THETA_RELATION(list_of_animals{fn}, pt,freq_range,f,wind_length,wind_overlap );

            for sub_j=1:numel(transt_psth)

                psth_structure(n_strctut) = transt_psth(sub_j);
                n_strctut = n_strctut+1;
            end
            animal_names = [animal_names;[repmat({list_of_animals{fn}},numel(transt_psth),1), repmat({pt},numel(transt_psth),1), num2cell(1:numel(transt_psth))' ]]

        end

    end


end

%% now saving

save([saving_folder,'\psth_structure_speed_theta_v2.mat'],'psth_structure')
save([saving_folder,'\animal_names_speed_theta_v2.mat'],'animal_names')

