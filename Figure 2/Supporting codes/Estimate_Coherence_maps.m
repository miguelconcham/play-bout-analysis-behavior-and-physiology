
npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];

animal_list = animal_list(11:12);
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


pli_struct = [];
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
        pli_struct = GENERATE_COHERNECE_MAPS_STRUCT([npx_Raw_Data, '\', animal_list(fn).name],Hd_freq1,Hd_freq2, 'PAG');
         
        n_strctut = n_strctut+numel(pli_struct);
        animal_names = [animal_names;[repmat(animal_list(fn).name,numel(pli_struct),1) num2cell(1:numel(pli_struct))']]
    else
        transt_psth =  GENERATE_COHERNECE_MAPS_STRUCT([npx_Raw_Data, '\', animal_list(fn).name],Hd_freq1,Hd_freq2, 'PAG');
      
        for sub_j=1:numel(transt_psth)
    
            pli_struct(n_strctut) = transt_psth(sub_j);
            n_strctut = n_strctut+1;
        end
        animal_names = [animal_names;[repmat({animal_list(fn).name},numel(transt_psth),1) num2cell(1:numel(transt_psth))' ]]

    end


end
pli_struct_mpfc = pli_struct;
clear pli_struct
%%
disp('saving')
save([saving_folder,'\coherence_structure_pli_struct_NPX2.mat'],'pli_struct_NPX2', '-v7.3');
save([saving_folder,'\coherence_structure_animal_names_pli_struct_NPX2.mat'],'animal_names');


%% loading 
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';

disp('loading')
load([saving_folder,'\coherence_structure_pli_struct_mpfc.mat'],'pli_struct_mpfc');
load([saving_folder,'\coherence_structure_pli_struct_B1D11007.mat'],'pli_struct_B1D11007');
load([saving_folder,'\coherence_structure.mat'],'pli_struct');
load([saving_folder,'\coherence_structure_animal_names.mat'],'animal_names');
load([saving_folder,'\coherence_structure_pli_struct_NPX2.mat'],'pli_struct_NPX2', '-v7.3');


%%
animal_names(7:numel(pli_struct)+3,1)  = animal_names(4:numel(pli_struct),1);
pli_struct(7:numel(pli_struct)+3) = pli_struct(4:numel(pli_struct));
pli_struct(4) = pli_struct(1);

pli_struct(5:6) = pli_struct(2:3);
pli_struct(2:3) =pli_struct_mpfc;
pli_struct(1) = pli_struct_B1D11007;
pli_struct(14) = pli_struct_NPX2(1);

animal_names{1,1} = 'B1D1 1007 Dual ORIGINALLY mPFC';
animal_names{2,1} = 'B1D1 1012 Dual ORIGINALLY mPFC';
animal_names{3,1} = 'B1D1 1013 Dual ORIGINALLY mPFC';
animal_names{4,1} = 'B1D1 1007 Dual ORIGINALLY PAG';
animal_names{5,1} = 'B1D1 1012 Dual ORIGINALLY PAG';
animal_names{6,1} = 'B1D1 1013 Dual ORIGINALLY PAG';
%% now some plots:  phase lags
figure
colormap(jet)
c_lim = [-2*pi 2*pi]/12;
for j=1:numel(pli_struct)
    subplot(5,3,j)

    [~, order] = sort(pli_struct(j).channel_map(2:end,2));
    matrix  =pli_struct(j).mean_phaseLag_freq2(order ,order);
    matrix((1:size(matrix,1)) + size(matrix,1)*((1:size(matrix,1))-1))
    if j ==4
        a = [ matrix(117:end,:);matrix(1:116,:)];
        a = [ a(:,117:end),a(:,1:116)];
        matrix = a;
    elseif j==3
        matrix = matrix ;
    elseif j==2
        a = [ matrix(116:end,:);matrix(1:115,:)];
        a = [ a(:,116:end),a(:,1:115)];
        matrix = a;
    end



    pcolor(size(matrix,1)+1 - (1:size(matrix,1)),1:size(matrix,1),matrix)
    hold on
    plot([size(matrix,1) 1],[1 size(matrix,1)],'w')
    shading interp
    xticks(1:16:384)
    xticklabels(flipud(pli_struct(j).areas_by_channel(order(1:16:384))))
    yticks(1:16:384)
    yticklabels(pli_struct(j).areas_by_channel(order(1:16:384)))
    axis xy
    clim(c_lim)
    title(animal_names{j,1})
end
%%
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';

disp('loading')
load([saving_folder,'\coherence_structure_pli_struct_mpfc.mat'],'pli_struct_mpfc');
load([saving_folder,'\coherence_structure_pli_struct_B1D11007.mat'],'pli_struct_B1D11007');
load([saving_folder,'\coherence_structure.mat'],'pli_struct');
load([saving_folder,'\coherence_structure_animal_names.mat'],'animal_names');

%%
