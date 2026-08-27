%% Paths, animals, and filter
data_root     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
npx_Raw_Data  = [data_root, '\NPX data\NPX raw data'];
saving_folder = [data_root, '\Analysis results\Theta psth'];
run('\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');
if ~exist(saving_folder, 'dir')
    mkdir(saving_folder);
end
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
freq_range_3    = [25 45]; %actaully detecting 25 hz oscillations
freq_range_4    = [12 20];
sr              = 2500;
filter_order    = 2000;

freq_range_5    = [30 90];


% 
% % Parameters for delta
% Hd_freq = designfilt('bandpassfir', ...
% 'FilterOrder', filter_order, ...
% 'CutoffFrequency1', freq_range_1(1), ...
% 'CutoffFrequency2', freq_range_1(2), ...
% 'SampleRate', sr, ...
% 'DesignMethod', 'window', ...
% 'Window', 'hamming');
% bin_size_freq = 0.01;

% Parameters for theta
% Hd_freq = designfilt('bandpassfir', ...
% 'FilterOrder', filter_order, ...
% 'CutoffFrequency1', freq_range_2(1), ...
% 'CutoffFrequency2', freq_range_2(2), ...
% 'SampleRate', sr, ...
% 'DesignMethod', 'window', ...
% 'Window', 'hamming');
% bin_size_freq = 0.001;


% % Parameters for low gamma (this is actualyl finding betha entraeinment)
% Hd_freq = designfilt('bandpassfir', ...
% 'FilterOrder', filter_order, ...
% 'CutoffFrequency1', freq_range_3(1), ...
% 'CutoffFrequency2', freq_range_3(2), ...
% 'SampleRate', sr, ...
% 'DesignMethod', 'window', ...
% 'Window', 'hamming');
% bin_size_freq = 0.001;


% % Parameters for bet (i dint find anything with this range)
% Hd_freq = designfilt('bandpassfir', ...
% 'FilterOrder', filter_order, ...
% 'CutoffFrequency1', freq_range_3(1), ...
% 'CutoffFrequency2', freq_range_3(2), ...
% 'SampleRate', sr, ...
% 'DesignMethod', 'window', ...
% 'Window', 'hamming');
% bin_size_freq = 0.001; 

% Parameters for low gamma (corrected)
Hd_freq = designfilt('bandpassfir', ...
'FilterOrder', filter_order, ...
'CutoffFrequency1', freq_range_5(1), ...
'CutoffFrequency2', freq_range_5(2), ...
'SampleRate', sr, ...
'DesignMethod', 'window', ...
'Window', 'hamming');
bin_size_freq = 0.001;
%%
tic
for fn = 1:numel(animal_list)

    if fn==1
        transt_psth  = GENERATE_PHASE_COUPLING_STRUCTURE([npx_Raw_Data, '\', animal_list(fn).name],Hd_freq,bin_size_freq );
        phase_struct = transt_psth;

        % save([saving_folder,'\',animal_list(fn).name,'_FreqRange_',num2str([Hd_freq.CutoffFrequency1 Hd_freq.CutoffFrequency2])...
        %     ,'_phase_couplig_structure.mat'],'transt_psth', '-v7.3');


        n_strctut = n_strctut+numel(phase_struct);
        animal_names = [animal_names;[repmat(animal_list(fn).name,numel(phase_struct),1) num2cell(1:numel(phase_struct))']]
    else
        transt_psth = GENERATE_PHASE_COUPLING_STRUCTURE([npx_Raw_Data, '\', animal_list(fn).name],Hd_freq,bin_size_freq );
        % save([saving_folder,'\',animal_list(fn).name,'_FreqRange_',num2str([Hd_freq.CutoffFrequency1 Hd_freq.CutoffFrequency2])...
        %     ,'_phase_couplig_structure.mat'],'transt_psth', '-v7.3');

        for sub_j=1:numel(transt_psth)

            phase_struct(n_strctut) = transt_psth(sub_j);
            n_strctut = n_strctut+1;
        end
        animal_names = [animal_names;[repmat({animal_list(fn).name},numel(transt_psth),1) num2cell(1:numel(transt_psth))' ]]

    end
    toc

end


%% save if needed
disp('saving')
% save([saving_folder,'\delta_phase_couplig_structure_updated_with_non_playbouts_trough_peaks.mat'],'phase_struct', '-v7.3');
% save([saving_folder,'\delta_phase_couplig_animal_names_updated_with_non_playbouts_trough_peaks.mat'],'animal_names');

% save([saving_folder,'\low_gamma_phase_couplig_structure_updated_with_non_playbouts_peak_peaks.mat'],'phase_struct', '-v7.3');
% save([saving_folder,'\low_gamma_phase_couplig_animal_names_updated_with_non_playbouts_peak_peaks.mat'],'animal_names');

save([saving_folder,'\gamma_couplig_structure_updated_with_non_playbouts_peak_peaks.mat'],'phase_struct', '-v7.3');
save([saving_folder,'\gamma_couplig_animal_names_updated_with_non_playbouts_peak_peaks.mat'],'animal_names');

%% load
data_root     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
saving_folder = [data_root, '\Analysis results\Theta psth'];

load([saving_folder,'\theta_phase_couplig_structure_updated_with_non_playbouts.mat'],'phase_struct');
load([saving_folder,'\theta_phase_couplig_animal_names_updated_with_non_playbouts.mat'],'animal_names');
% phase_struct(6 )=[];
% animal_names(6,:)=[];
%%

 phase_prop_names = {'PreferedAngle','MVL','MVLPval','PPC','PPCPval','MeanRate', 'Id'};

all_session_phase_stats = [];
all_session_psth        = [];

all_neurons = [];
for ns = 1:numel(phase_struct)

     all_session_phase_stats = cat(2,all_session_phase_stats, phase_struct(ns).session_phase_stats(1:2,:,:));
     all_session_psth = cat(2,all_session_psth, phase_struct(ns).session_psth(1:2,:,:));
     this_session_cluster_info = phase_struct(ns).cluster_info;

     sub_Table =  array2table(squeeze(phase_struct(ns).session_phase_stats(1,:,:)));
     sub_Table.Properties.VariableNames =phase_prop_names;
     this_session_cluster_info.Partner1 =sub_Table;

     sub_Table =  array2table(squeeze(phase_struct(ns).session_phase_stats(2,:,:)));
     sub_Table.Properties.VariableNames = phase_prop_names;
     this_session_cluster_info.Partner2 =sub_Table;


     sub_Table =  array2table(squeeze(phase_struct(ns).play_phase_stats(1,:,:)));
     sub_Table.Properties.VariableNames =phase_prop_names;
     this_session_cluster_info.Play =sub_Table;

     sub_Table =  array2table(squeeze(phase_struct(ns).pre_play_phase_stats(1,:,:)));
     sub_Table.Properties.VariableNames = phase_prop_names;
     this_session_cluster_info.PrePlay =sub_Table;


      sub_Table =  array2table(squeeze(phase_struct(ns).non_play_phase_stats(1,:,:)));
     sub_Table.Properties.VariableNames =phase_prop_names;
     this_session_cluster_info.NonPlay =sub_Table;

     sub_Table =  array2table(squeeze(phase_struct(ns).pre_non_play_phase_stats(1,:,:)));
     sub_Table.Properties.VariableNames = phase_prop_names;
     this_session_cluster_info.PreNonPlay =sub_Table;


     this_session_cluster_info.session = repmat(animal_names(ns,1),size(this_session_cluster_info,1),1);


     sub_Table =  array2table(squeeze(phase_struct(ns).entire_recording_phase_stats(1,:,:)));
     sub_Table.Properties.VariableNames = phase_prop_names;
     this_session_cluster_info.EntireSession =sub_Table;

     this_session_cluster_info.session = repmat(animal_names(ns,1),size(this_session_cluster_info,1),1);
     all_neurons = [all_neurons; this_session_cluster_info];
end

%%

save([saving_folder,'\theta_all_neurons_v2.mat'],'all_neurons');
%% ploting rate change
all_neurons.area(ismember(all_neurons.area, {'isRT'})) =     {'isRt'  };
% area_list = unique(all_neurons.area)';
condition_1 = 'NonPlay';
condition_2 = 'PreNonPlay';
% area_list =  {'4N'	'DLPAG'	'DMPAG'	'DR'	'InfCol'	'LPAG'	'LSD'	'LSI'	'SupCol'	'VLPAG'	'isRt'	'mlf'};
alpha = 0.01;
area_list = {'SupCol'  'DLPAG'	'LPAG' 'VLPAG'	'DR'};
col = strcmp(phase_prop_names, 'MeanRate');
figure
concatenated_rate_differences = [];
mean_rate_differences = [];
d_prime = [];

for an=1:numel(area_list)
    subplot(1,numel(area_list),an)
    % index = strcmp(all_neurons.area, area_list{an}) & all_neurons.(condition_1).PPCPval<alpha  & all_neurons.(condition_2).PPCPval<alpha & ~isnan(all_neurons.(condition_1).PPC + all_neurons.(condition_2).PPC);
     index = strcmp(all_neurons.area, area_list{an}) & all_neurons.EntireSession.PPCPval<alpha &  ~isnan(all_neurons.(condition_1).PPC + all_neurons.(condition_2).PPC);
    x_jit = (rand(sum(index),2) -.5)*.25;
    matrix2plot = 100*(squeeze(all_session_phase_stats(:,index,col))-repmat(all_neurons.fr(index)',2,1))./repmat(all_neurons.fr(index)',2,1);
    plot((repmat([1 2],sum(index),1)+x_jit)',matrix2plot,':k')
    data = [diff(matrix2plot)' ones(size(matrix2plot,2),1)*an];
    if size(data,1)>1
    mean_rate_differences = [mean_rate_differences;mean(data, 'omitmissing')];    
    else
          mean_rate_differences = [mean_rate_differences;data];
    end


     d_prime = [d_prime;[2*mean(diff(matrix2plot))/sqrt(sum(std(matrix2plot,[],2))) an]];    
    concatenated_rate_differences = [concatenated_rate_differences;[diff(matrix2plot)' ones(size(matrix2plot,2),1)*an all_neurons.EntireSession.PreferedAngle(index)]];
    hold on
        plot((repmat([1 2],sum(index),1)+x_jit)',matrix2plot,'.k', 'MarkerSize',3)

    hold on
    plot([1 2]',mean(matrix2plot,2, 'omitmissing'),'_r', 'MarkerSize',4)
    if all(any(isnan(squeeze(all_session_phase_stats(:,index,col))),1))
        p = 1;
    else
        [p,h,t] = signrank(matrix2plot(1,:), matrix2plot(2,:));
        % [h,p,t] = ttest(matrix2plot(1,:), matrix2plot(2,:))
    end
    ylim([-200 200])
    title([area_list{an},  ' ', num2str(p)])
end
%%

figure

subplot(3,1,2)
subplot(1,2,1)
swarmchart(concatenated_rate_differences(:,2), concatenated_rate_differences(:,1), 'k.')
hold on
plot(mean_rate_differences(:,2),mean_rate_differences(:,1), '_r')

plot([0 numel(area_list)+1], [0 0], 'r')
xticks(1:numel(area_list))
xticklabels(area_list)
xlim tight
ylim([-100 100])

subplot(1,2,2)
bar(d_prime(:,1))
xticks(1:numel(area_list))
xticklabels(area_list)
xlim tight

%% ploting entreinment change per area
area_list = unique(all_neurons.area);
alpha= 0.01;
area_list = {'SupCol'  'DLPAG'	'LPAG' 'VLPAG'	'DR'};

col = strcmp(phase_prop_names, 'PPCPval');
entreinment_per_partner = nan(numel(area_list),6);
figure
for an=1:numel(area_list)
    subplot(1,numel(area_list),an)
    index = strcmp(all_neurons.area, area_list{an});
    x_jit = (rand(sum(index),2) -.5)*.25;
    semilogy((repmat([1 2],sum(index),1)+x_jit)',squeeze(all_session_phase_stats(:,index,col)),':k')
    p1_bool     = all_neurons.(condition_1).PPCPval(index)<alpha;
    p2_bool     = all_neurons.(condition_2).PPCPval(index)<alpha;
    all_bool    = all_neurons.EntireSession.PPCPval(index)<alpha;
    entreinment_per_partner(an,:) = [sum(p1_bool & p2_bool) sum(p1_bool & ~p2_bool) sum(~p1_bool & p2_bool) sum(~p1_bool & ~p2_bool & all_bool) sum(~p1_bool & ~p2_bool & ~all_bool) sum(index)*sum(index)]/sum(index);
    hold on
semilogy([0 2.5],[alpha alpha],'r')
title(area_list{an})
end
%%

%%
all_neurons.area(ismember(all_neurons.area, {'isRT'})) =     {'isRt'  };
% area_list = unique(all_neurons.area)';

% area_list =  {'4N'	'DLPAG'	'DMPAG'	'DR'	'InfCol'	'LPAG'	'LSD'	'LSI'	'SupCol'	'VLPAG'	'isRt'	'mlf'};
% area_list =  {	'DLPAG'	'DMPAG'	'DR'	'LPAG'		'SupCol'	'VLPAG'	'isRt'	'mlf'};
% area_list = {'SupCol'  'DLPAG'	'LPAG' 'VLPAG'	'DR', 'isRt'	};
y_lim = [10^-6 5]
area_list =  {	'DLPAG'	'DMPAG'	'DR'	'LPAG'		'SupCol'	'VLPAG'	};
col = strcmp(phase_prop_names, 'PPC');
figure
concatenated_phase_differences = [];
mean_phase_differences = [];
% condition_1 = 'Partner1';
% condition_2 = 'Partner2';

% condition_1 = 'NonPlay';
% condition_2 = 'PreNonPlay';
condition_1 = 'Play';
condition_2 = 'PrePlay';
for an=1:numel(area_list)
    subplot(1,numel(area_list),an)
    index = strcmp(all_neurons.area, area_list{an}) & all_neurons.EntireSession.PPCPval<alpha &  ~isnan(all_neurons.(condition_1).PPC + all_neurons.(condition_2).PPC);
    x_jit = (rand(sum(index),2) -.5)*.25;
    matrix2plot = [all_neurons.(condition_1).PPC(index) all_neurons.(condition_2).PPC(index)  ];

    data = [diff(matrix2plot')' ones(size(matrix2plot,1),1)*an ] ;
    if size(data,1)>1
        mean_phase_differences = [mean_phase_differences;mean(data, 'omitmissing')];
    else
        mean_phase_differences = [mean_phase_differences;data];
    end




    concatenated_phase_differences = [concatenated_phase_differences;[diff(matrix2plot')' ones(size(matrix2plot,1),1)*an all_neurons.EntireSession.PreferedAngle(index) matrix2plot]];


    semilogy((repmat([1 2],sum(index),1)+x_jit)',matrix2plot',':k')
    hold on
    semilogy((repmat([1 2],sum(index),1)+x_jit)',matrix2plot','.k', 'MarkerSize',3)

    semilogy([1 2]',mean(matrix2plot,1, 'omitmissing'),'_r', 'MarkerSize',4, 'LineWidth',3)
    if all(any(isnan(matrix2plot),1))
        p = 1;
    else
        [p,h,t] = signrank(matrix2plot(1,:)',matrix2plot(2,:)');
        % [h,p] = ttest(matrix2plot(1,:)',matrix2plot(2,:)');
    end
    ylim(y_lim)
    title([area_list{an},  ' ', num2str(p)])
end
%% old PLI plot

%% rate and entraeinmetn change by cell type(peak and trough)



%% same but in paired comparisson (BAR PLOTS)
area_list = {'SupCol'  'DLPAG'	'LPAG' 'VLPAG'	'DR' };
alpha=0.05;
figure
y_lim_ppc = [-0.1 0.4];
y_lim_Rate = [-100 100];
line_width = .25;
names = {'Trough','Peak'};
bar_per_area = [];
for an=1:numel(area_list)
    subplot(1,numel(area_list)+1,an)
    index           = concatenated_phase_differences(:,2)==an;
    angle_index     = concatenated_phase_differences(:,3)<-pi/2 | concatenated_phase_differences(:,3)>pi/2;


    p =ranksum(concatenated_phase_differences(angle_index & index,1),concatenated_phase_differences(~angle_index & index,1));
    p1 = signrank(concatenated_phase_differences(angle_index & index,1));
    % [h,p1]= ttest(concatenated_phase_differences(angle_index & index,end-1),concatenated_phase_differences(angle_index & index,end))

   


    rand_x = (rand(sum(angle_index & index),2)-.5)/2;
    plot((rand_x + ones(size(rand_x))*diag([1 2]))',concatenated_phase_differences(angle_index & index, end-1:end)', 'k:')
    hold on
    plot((rand_x + ones(size(rand_x))*diag([1 2]))',concatenated_phase_differences(angle_index & index, end-1:end)', 'k.')
    this_medians = median(concatenated_phase_differences(angle_index & index, end-1:end))
    if p1<alpha
        text(1.5,.9*y_lim_ppc(2), num2str(round(p1,4)))
        plot([1 2], y_lim_ppc([2 2])*.85, 'r')
    end
    bar(1,this_medians(1) , 'r', 'FaceAlpha',.7, 'EdgeColor','none')
    bar(2,this_medians(2) , 'b', 'FaceAlpha',.7, 'EdgeColor','none')

    yscale log
    axis tight


    rand_x = (rand(sum(~angle_index & index),2)-.5)/2;
    plot((rand_x + ones(size(rand_x))*diag([3 4]))',concatenated_phase_differences(~angle_index & index, end-1:end)', 'k:')
    hold on
    plot((rand_x + ones(size(rand_x))*diag([3 4]))',concatenated_phase_differences(~angle_index & index, end-1:end)', 'k.')


    p2 = signrank(concatenated_phase_differences(~angle_index & index, end-1),concatenated_phase_differences(~angle_index & index, end));
    % [h,p2]= ttest(concatenated_phase_differences(~angle_index & index, end-1)+.2,concatenated_phase_differences(~angle_index & index, end)+.2);
  
     this_medians = median(concatenated_phase_differences(~angle_index & index, end-1:end));
    if p2<alpha
        text(3.5,.9*y_lim_ppc(2), num2str(round(p2,4)))
        plot([3 4], y_lim_ppc([2 2])*.85, 'r')
    end
       bar(3,this_medians(1) , 'r', 'FaceAlpha',.7, 'EdgeColor','none')
        bar(4,this_medians(2) , 'b', 'FaceAlpha',.7, 'EdgeColor','none')
        yscale log
   axis tight
    
    

ylim(y_lim_ppc)
xticks(1:4)
xticklabels({[names{1}, ' ', condition_1],[names{1}, ' ', condition_2],[names{2} , ' ', condition_1],[names{2}, ' ', condition_2]})
   
  
    
 

end


 subplot(1,numel(area_list)+1,numel(area_list)+1)
    index           = true(size(concatenated_phase_differences(:,2)));
    angle_index     = concatenated_phase_differences(:,3)<-pi/2 | concatenated_phase_differences(:,3)>pi/2;


    p =ranksum(concatenated_phase_differences(angle_index & index,1),concatenated_phase_differences(~angle_index & index,1));
    p1 = signrank(concatenated_phase_differences(angle_index & index,1));
    % [h,p1]= ttest(concatenated_phase_differences(angle_index & index,end-1),concatenated_phase_differences(angle_index & index,end))

   
    y2plot = atanh(concatenated_phase_differences(angle_index & index, end-1:end)');


    rand_x = (rand(sum(angle_index & index),2)-.5)/2;
    plot((rand_x + ones(size(rand_x))*diag([1 2]))',y2plot, 'k:')
    hold on
    plot((rand_x + ones(size(rand_x))*diag([1 2]))',y2plot, 'k.')
    this_medians = median(y2plot)
    if p1<alpha
        text(1.5,.9*y_lim_ppc(2), num2str(round(p1,8)))
        plot([1 2], y_lim_ppc([2 2])*.85, 'r')
    end
    bar(1,this_medians(1) , 'r', 'FaceAlpha',.7, 'EdgeColor','none')
    bar(2,this_medians(2) , 'b', 'FaceAlpha',.7, 'EdgeColor','none')

    yscale log
    axis tight


    rand_x = (rand(sum(~angle_index & index),2)-.5)/2;
    y2plot =atanh(concatenated_phase_differences(~angle_index & index, end-1:end)');
    plot((rand_x + ones(size(rand_x))*diag([3 4]))',y2plot, 'k:')
    hold on
    plot((rand_x + ones(size(rand_x))*diag([3 4]))',y2plot, 'k.')


    p2 = signrank(concatenated_phase_differences(~angle_index & index, end-1),concatenated_phase_differences(~angle_index & index, end));
    % [h,p2]= ttest(concatenated_phase_differences(~angle_index & index, end-1)+.2,concatenated_phase_differences(~angle_index & index, end)+.2);
  
     this_medians = median(y2plot);
    if p2<alpha
        text(3.5,.9*y_lim_ppc(2), num2str(round(p2,8)))
        plot([3 4], y_lim_ppc([2 2])*.85, 'r')
    end
       bar(3,this_medians(1) , 'r', 'FaceAlpha',.7, 'EdgeColor','none')
        bar(4,this_medians(2) , 'b', 'FaceAlpha',.7, 'EdgeColor','none')
        % yscale log
   axis tight
    
    

ylim(y_lim_ppc)
xticks(1:4)
xticklabels({[names{1}, ' ', condition_1],[names{1}, ' ', condition_2],[names{2} , ' ', condition_1],[names{2}, ' ', condition_2]})
   
  %%
  %% same but in paired comparisson (SQUARE PLOTS)

alpha=0.05;
figure('units','normalized','outerposition',[0 0 1 1]);
y_lim_ppc = [-0.005 0.07];
y_lim_Rate = [-100 100];
line_width = .25;
names = {'Trough','Peak'};
bar_per_area = [];
for an=1:numel(area_list)
    subplot(2,numel(area_list)+1,an)
    index           = concatenated_phase_differences(:,2)==an;
    angle_index     = concatenated_phase_differences(:,3)<-pi/2 | concatenated_phase_differences(:,3)>pi/2;


    p =ranksum(concatenated_phase_differences(angle_index & index,1),concatenated_phase_differences(~angle_index & index,1));
    p1 = signrank(concatenated_phase_differences(angle_index & index,1));
    % [h,p1]= ttest(concatenated_phase_differences(angle_index & index,end-1),concatenated_phase_differences(angle_index & index,end))

   


    y2plot = atanh(concatenated_phase_differences(angle_index & index, end-1:end));
    plot(y2plot(:,1),y2plot(:,2), 'k.')
    hold on
    this_medians = median(y2plot);
     plot(y_lim_ppc, y_lim_ppc, 'r')

    axis([y_lim_ppc y_lim_ppc])
    plot(this_medians(1), this_medians(2), 'xr')

        title([area_list{an},' ' num2str(round(p1,4))])
   

axis square
if an==1
        ylabel({'Trough', condition_2})
end
xlabel(condition_1)




      subplot(2,numel(area_list)+1,numel(area_list)+1 +an)
  
     y2plot = atanh(concatenated_phase_differences(~angle_index & index, end-1:end));
    plot(y2plot(:,1),y2plot(:,2), 'k.')
    hold on
    this_medians = median(y2plot);
    plot(this_medians(1), this_medians(2), 'xr')
   plot(y_lim_ppc, y_lim_ppc, 'r')

    axis([y_lim_ppc y_lim_ppc])

    p2 = signrank(concatenated_phase_differences(~angle_index & index, end-1),concatenated_phase_differences(~angle_index & index, end));
    % [h,p2]= ttest(concatenated_phase_differences(~angle_index & index, end-1)+.2,concatenated_phase_differences(~angle_index & index, end)+.2);
  
     this_medians = median(concatenated_phase_differences(~angle_index & index, end-1:end));
    if p2<alpha
       title(num2str(round(p2,4)))
       
    end
      
axis square
if an==1
    ylabel({'Peak', condition_2})
end
xlabel(condition_1)


    
    

% ylim(y_lim_ppc)
% xticks(1:4)
% xticklabels({[names{1}, ' ', condition_1],[names{1}, ' ', condition_2],[names{2} , ' ', condition_1],[names{2}, ' ', condition_2]})
% 
  
    
 

end


 subplot(2,numel(area_list)+1,numel(area_list)+1)
    index           = true(size(concatenated_phase_differences(:,2)));
    angle_index     = concatenated_phase_differences(:,3)<-pi/2 | concatenated_phase_differences(:,3)>pi/2;


 
    p1 = signrank(concatenated_phase_differences(angle_index & index, end),concatenated_phase_differences(angle_index & index, end-1));
    % [h,p1]= ttest(concatenated_phase_differences(angle_index & index,end-1),concatenated_phase_differences(angle_index & index,end))

   

    y2plot = atanh(concatenated_phase_differences(angle_index & index, end-1:end));
    plot(y2plot(:,1),y2plot(:,2), 'k.')
    hold on
    this_medians = median(y2plot);
   plot(y_lim_ppc, y_lim_ppc, 'r')

    axis([y_lim_ppc y_lim_ppc])
    plot(this_medians(1), this_medians(2), 'xr')
    if p1<alpha
        title(num2str(round(p1,6)))
    end
axis square
 subplot(2,numel(area_list)+1,2*(numel(area_list)+1))
    index           = true(size(concatenated_phase_differences(:,2)));
    angle_index     = concatenated_phase_differences(:,3)<-pi/2 | concatenated_phase_differences(:,3)>pi/2;



    p1 = signrank(concatenated_phase_differences(~angle_index & index, end),concatenated_phase_differences(~angle_index & index, end-1));
    % [h,p1]= ttest(concatenated_phase_differences(angle_index & index,end-1),concatenated_phase_differences(angle_index & index,end))

   

    y2plot = atanh(concatenated_phase_differences(~angle_index & index, end-1:end));
    plot(y2plot(:,1),y2plot(:,2), 'k.')
    hold on
    this_medians = median(y2plot);
   plot(y_lim_ppc, y_lim_ppc, 'r')

    axis([y_lim_ppc y_lim_ppc])
    plot(this_medians(1), this_medians(2), 'xr')
    if p1<alpha
        title(num2str(round(p1,6)))
    end

axis square
%%


y = [concatenated_phase_differences(angle_index , end-1);concatenated_phase_differences(angle_index , end)];
factor1 = [concatenated_phase_differences(angle_index , 2);concatenated_phase_differences(angle_index , 2)];
factor2 = [concatenated_phase_differences(angle_index , end-1)*0;concatenated_phase_differences(angle_index , end)*0+1];

no_nan = y>0;
[p, tbl, stats] = anovan(log(y(no_nan)), {factor1(no_nan), factor2(no_nan)}, 'model', 'interaction', ...
    'varnames', {'Factor1', 'Factor2'});
%%

area_list = {'SupCol'};

figure
col = 4;
    subplot(1,3,1)
    index = strcmp(all_neurons.area, area_list{1});
    x_jit = (rand(sum(index),2) -.5)*.25;

     matrix2plot = squeeze(all_session_phase_stats(:,index,col));
    plot((repmat([1 2],sum(index),1)+x_jit)',matrix2plot,':k')
    hold on
    plot((repmat([1 2],sum(index),1)+x_jit)',matrix2plot,'.k', 'MarkerSize',3)

    plot([1 2]',mean(matrix2plot, 2,'omitmissing'),'_r', 'MarkerSize',4)
    plot([1 2]',mean(matrix2plot, 2,'omitmissing'),'r', 'LineWidth',3)

    if all(any(isnan(matrix2plot),1))
        p = 1;
    else
        [p,h,t] = signrank(matrix2plot(1,:)',matrix2plot(2,:)');
    end
    ylim([-.25 .25])
    title([area_list{1},  ' ', num2str(p)])

    col = 6;
    subplot(1,3,2)
    matrix2plot = squeeze(all_session_phase_stats(:,index,col))-repmat(all_neurons.fr(index)',2,1);
    plot((repmat([1 2],sum(index),1)+x_jit)',matrix2plot,':k')
    hold on
    plot((repmat([1 2],sum(index),1)+x_jit)',matrix2plot,'.k', 'MarkerSize',3)
    plot([1 2]',mean(matrix2plot, 2,'omitmissing'),'r', 'LineWidth',3)

   
    plot([1 2]',mean(matrix2plot,2, 'omitmissing'),'_r', 'MarkerSize',4)
    if all(any(isnan(squeeze(all_session_phase_stats(:,index,col))),1))
        p = 1;
    else
        [p,h,t] = signrank(matrix2plot(1,:), matrix2plot(2,:))
    end
    ylim([-5 5])
    title([area_list{1},  ' ', num2str(p)])




%%
figure
z_score_limit = 1;
index2include = abs(zscore(concatenated_rate_differences(:,1)))<z_score_limit & abs((concatenated_phase_differences(:,1) - mean(concatenated_phase_differences(:,1), 'omitmissing'))/std(concatenated_phase_differences(:,1), 'omitmissing'))<z_score_limit;
plot(concatenated_rate_differences(index2include,1), concatenated_phase_differences(index2include,1), 'k.')
no_nan = ~any(isnan([concatenated_rate_differences(:,1) concatenated_phase_differences(:,1)]),2);
[c,p]=corr(concatenated_rate_differences(no_nan & index2include,1), concatenated_phase_differences(no_nan & index2include,1))
xlabel('% Rate Change (Partner2 - Partner1)')
ylabel('PPC change (Partner2 - Partner1)')


%% entreinment per region

%% delta and theta comparisson

load([saving_folder,'\theta_all_neuronsv2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) =     {'isRt'  };
all_neurons_theta = all_neurons;
load([saving_folder,'\delta_all_neurons.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) =     {'isRt'  };
all_neurons_delta = all_neurons;


theta_partner_1 = all_neurons_theta.(condition_1).Pval;



alpha = 0.01;


% area_list = unique(all_neurons.area)';

% area_list =  {'4N'	'DLPAG'	'DMPAG'	'DR'	'InfCol'	'LPAG'	'LSD'	'LSI'	'SupCol'	'VLPAG'	'isRt'	'mlf'};

area_list = {'DLPAG'	'DMPAG'	'DR'	'InfCol'	'LPAG'		'SupCol'	'VLPAG'	'isRt'	'mlf'};


percentages_partner_1 = zeros(numel(area_list),4);
percentages_partner_2 = zeros(numel(area_list),4);


for an=1:numel(area_list) 
    index = strcmp(all_neurons.area, area_list{an});
    percentages_partner_1(an,1) = sum(all_neurons_theta.(condition_1).Pval(index)<alpha & all_neurons_delta.(condition_1).Pval(index)>alpha)/sum(~isnan(all_neurons_theta.(condition_1).Pval(index)));
    percentages_partner_1(an,2) = sum(all_neurons_delta.(condition_1).Pval(index)<alpha & all_neurons_theta.(condition_1).Pval(index)>alpha)/sum(~isnan(all_neurons_delta.(condition_1).Pval(index)));
    percentages_partner_1(an,3) =  sum(all_neurons_delta.(condition_1).Pval(index)<alpha & all_neurons_theta.(condition_1).Pval(index)<alpha)/sum(~isnan(all_neurons_delta.(condition_1).Pval(index)));

    percentages_partner_2(an,1) = sum(all_neurons_theta.(condition_2).Pval(index)<alpha & all_neurons_delta.(condition_2).Pval(index)>alpha)/sum(~isnan(all_neurons_theta.(condition_2).Pval(index)));
    percentages_partner_2(an,2) = sum(all_neurons_delta.(condition_2).Pval(index)<alpha & all_neurons_theta.(condition_2).Pval(index)>alpha)/sum(~isnan(all_neurons_delta.(condition_2).Pval(index)));
    percentages_partner_2(an,3) =  sum(all_neurons_delta.(condition_2).Pval(index)<alpha & all_neurons_theta.(condition_2).Pval(index)<alpha)/sum(~isnan(all_neurons_delta.(condition_2).Pval(index)));
end





percentages_partner_1(:,4) = 1-sum(percentages_partner_1,2);
percentages_partner_2(:,4) = 1-sum(percentages_partner_2,2);


figure
subplot(1,2,1)
bar(percentages_partner_1, 'stacked')
xticks(1:numel(area_list))
xticklabels(area_list)


subplot(1,2,2)
bar(percentages_partner_2, 'stacked')
xticks(1:numel(area_list))
xticklabels(area_list)







