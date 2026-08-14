%%  define folders
repo_root       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\Codes repository';
data_root       = [repo_root, '\Data'];
figure_dir      = [repo_root, '\Figure 1\outputs'];
synch_folder    = [data_root, '\Synch data'];
hmm_raw_data    = [data_root, '\HMM data\HMM raw data'];
call_folder     = [data_root, '\CallDetectionBackup'];
behavior_folder = [data_root, '\Behavior backups'];
analyssis_folder = [data_root, '\Analysis results\HMM 2 and 3 states 2 partners'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig1B Example call rate speed for  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

animal_name     = 'B1D1 1013 Dual';
repeated_animal = strsplit(animal_name, ' ');
repeated_animal = repeated_animal{end};
beh_bin = 0.01; %for hmm stimate
start_event = 480;
end_event   = 490;

load([synch_folder, '\',animal_name, '\synch_model_video2audio.mat'], 'synch_model_video2audio')

Behavior            =   readtable([behavior_folder,'\',animal_name, '.txt']); 
Behavior(:,2)       = [];
Behavior.Properties.VariableNames = {'Animal', 'Start', 'End', 'Length', 'Type'};
Behavior.Start      = predict(synch_model_video2audio, Behavior.Start);
Behavior.End        = predict(synch_model_video2audio, Behavior.End);
Behavior.Type2      = Behavior.Type;
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_A', 'Pounce_B'})) = {'Pounce'}; %% Merging behaviors to Type2
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
Behavior.Type2(ismember( Behavior.Type2,'')) = {'Other'};
Behavior(ismember(Behavior.Animal, 'Reversal'),:) = [];


min_time2analysis = Behavior.Start(ismember(Behavior.Type,'Partners session'))  ;
max_time2analysis = Behavior.End(ismember(Behavior.Type,'Partners session'))  ;

ListOfPartners = [min_time2analysis max_time2analysis];
[~, time_ordered] = sort(ListOfPartners(:,1));
ListOfPartners = ListOfPartners(time_ordered,:);

pt = find(ListOfPartners(:,1)<=start_event & ListOfPartners(:,2)>=start_event);

load([hmm_raw_data,'\',animal_name, ' P', num2str(pt), '_PropAndTime'], 'adjusted_binned_time', 'all_properties', 'ALL_VARIABLE_NAMES')

variable_name2states        =  [hmm_raw_data, '\',animal_name, ' P', num2str(pt),'_states_K2.npy'];
% variable_behavior_prop      =  [hmm_raw_data, '\',animal_name, ' P', num2str(pt),'.npy'];
time_limit = (adjusted_binned_time>= ListOfPartners(pt,1) & adjusted_binned_time<=  ListOfPartners(pt,2) );
% time_limit = (adjusted_binned_time>= min(min(play_bouts_table)) & adjusted_binned_time<=   max(max(play_bouts_table)) );

hmm_states      =  readNPY(variable_name2states);
% behavior_prop   = readNPY(variable_behavior_prop);


sr          = 250000;

CallStats   = readtable([call_folder, '\',animal_name,'_Stats.xlsx']);
CallStats.Properties.VariableNames = cellfun(@(x) strrep(x, '_', ''),CallStats.Properties.VariableNames, 'UniformOutput',false );

[audio_data, sr] =  audioread([call_folder, '\', animal_name, '.wav'], round([start_event*sr  end_event*sr ]));



time_index = adjusted_binned_time>=start_event & adjusted_binned_time<=end_event;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% estiamte call spectrogram (take some time) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 250000;
f = 20000:100000;
window_insec        = 0.02;
overlap_insec       = 0.018;

window              = round(window_insec*fs);
noverlap            = round(overlap_insec*fs);

[s,f,t] = spectrogram(audio_data,window,noverlap,f,fs);
play_song([],[],[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% smooth image (also take some time) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pow = log10(abs(s));
% Parameters
sigma = std(pow(:)) ;       % Standard deviation
kernel_size = ceil(2*6*sigma);  % Kernel size (usually ~6*sigma for full support)

% Create grid
[x, y] = meshgrid(-kernel_size:kernel_size, -kernel_size:kernel_size);
gaussian_kernel = exp(-(x.^2 + y.^2) / (2*sigma^2));
gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));  % Normalize

smoothed_data = conv2(pow, gaussian_kernel, 'same');
I_sharp = imsharpen(smoothed_data);

h = fspecial('unsharp');         % Unsharp mask (a high-pass filter)
I_highcontrast = imfilter(smoothed_data, h);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% add call subplot to figure (PANEL B) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

behaviors_in_range = find((Behavior.Start>=start_event & Behavior.Start<=end_event) | (Behavior.End>=start_event & Behavior.End<=end_event))';
figure
colormap(1-gray)
subplot(4,1,1)

mean_call_length = median(CallStats.CallLengths);

imagesc(t+start_event, f, I_sharp)
axis xy
clim([-3.5 -.5])
ylim([35 80]*1000)
xlim([start_event end_event])
set(gca, 'TickDir', 'out')
xticks([start_event start_event+1])
xticklabels([])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% call rate, and other variables (PANEL B) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mov_win_sec = .5;
mov_win = mov_win_sec/beh_bin;
subplot(4,1,2)
hold off
var_2_plot      = {'NumCalls'};
var2_plot_index = find(ismember(ALL_VARIABLE_NAMES,var_2_plot));
plot(adjusted_binned_time(time_index), movmean(all_properties(time_index,var2_plot_index(1)),mov_win)/(mean_call_length/mov_win_sec), 'k')
ylabel('Call Rate (Hz)')
xlim([start_event end_event])
set(gca, 'TickDir', 'out')
xticks([start_event start_event+1])
xticklabels([])

var_2_plot      = {'AnimalSpeed' , 'PartnerSpeed'};
var2_plot_index = find(ismember(ALL_VARIABLE_NAMES,var_2_plot));

subplot(4,1,3)
set(gca, 'TickDir', 'out')
hold off

plot(adjusted_binned_time(time_index), all_properties(time_index,var2_plot_index(1)), 'k')
hold on
y_lim1 = ylim;
ylabel('Speed (a.u.)')


subplot(4,1,4)
hold off
plot(adjusted_binned_time(time_index), all_properties(time_index,var2_plot_index(2)), 'k')
hold on
y_lim2 = ylim;
xlim([start_event end_event])
ylabel('Speed (a.u.)')


hold on

for bn=behaviors_in_range
        beh_start   = Behavior.Start(bn);
        beh_end     = Behavior.End(bn);
        animal_type = Behavior.Animal(bn);
        behavior_type = Behavior.Type{bn};

        if ismember(animal_type,repeated_animal)
            subplot(4,1,3)
            fill([beh_start beh_end beh_end beh_start], y_lim1([1 1 2 2]), 'k', 'FaceAlpha',.5, 'EdgeColor','k')
            text((beh_end+beh_start)/2, mean(y_lim1), behavior_type, 'Color', 'r')
            xlim([start_event end_event])
            set(gca, 'TickDir', 'out')
            xticks([start_event start_event+1])
            xticklabels([])
        else
            subplot(4,1,4)
            fill([beh_start beh_end beh_end beh_start], y_lim2([1 1 2 2]), 'k', 'FaceAlpha',.5, 'EdgeColor','k')
            text((beh_end+beh_start)/2, mean(y_lim1), behavior_type,  'Color', 'r')
            xlim([start_event end_event])
            set(gca, 'TickDir', 'out')
            xticks([start_event start_event+1])    

        end
end

% SAVE Fig1B %%
% print(gcf,'-vector','-dsvg',[figure_dir,'\raw data scheme.svg'])
%% loading basic hmm data
 last_hmm_dir = analyssis_folder;


prediction_struct_files = dir([last_hmm_dir, '\* prediction_struct*']);
is_there_play_bout      = [];
is_there_hmm            = [];
is_this_hmm             = [];
filled_play_bouts       = [];
tripple_states          = [];
re_assignment           = [];
is_there_play_beh       = [];
session_index           = [];



total_number_of_hmm = 0;
percentage_of_bouts_with_play =nan(numel(prediction_struct_files),1);

for fn= 1:numel(prediction_struct_files)

    figure
     load([last_hmm_dir,'\',prediction_struct_files(fn).name]) 
     psth_edges = prediction_struct.psth_edges;
    [hmm_length_ordered, pb_order] = sort(diff(prediction_struct.filled_play_bouts'));
    imagesc(psth_edges, 1:numel(hmm_length_ordered), 1-cat(3,prediction_struct.is_there_hmm(pb_order,:), prediction_struct.is_there_play_beh(pb_order,:), prediction_struct.is_there_play_beh(pb_order,:) ))
    title(sum(sum(prediction_struct.is_this_hmm & prediction_struct.is_there_play_beh,2)>0))
    percentage_of_bouts_with_play(fn) = sum(sum(prediction_struct.is_this_hmm & prediction_struct.is_there_play_beh,2)>0)/size(prediction_struct.is_this_hmm,1);
    axis xy
    hold on
    plot(hmm_length_ordered,1:numel(hmm_length_ordered),'k')
    plot(hmm_length_ordered*0,1:numel(hmm_length_ordered),'k')
    
     is_there_play_bout      = [is_there_play_bout;prediction_struct.is_there_play_bout];
     is_there_play_beh       = [is_there_play_beh;prediction_struct.is_there_play_beh];
     is_there_hmm            = [is_there_hmm;prediction_struct.is_there_hmm];
     is_this_hmm             = [is_this_hmm;prediction_struct.is_this_hmm];
     filled_play_bouts       = [filled_play_bouts;prediction_struct.filled_play_bouts];
     tripple_states          = [tripple_states;prediction_struct.what_3states_is];
     re_assignment           = [re_assignment;prediction_struct.re_assignment]
     session_index           = [session_index;ones(size(prediction_struct.filled_play_bouts,1),1)*fn];
end
% current_hmm = tripple_states==3;
psth_edges = prediction_struct.psth_edges;


%% Fig S2E

filled_hmm_states = filled_play_bouts;

A = is_there_hmm;
B = is_there_play_beh;

[m, n] = size(A);

C = zeros(size(A));

for r = 1:m
    a = A(r,:);
    b = B(r,:);

    % Find start + end indices of sequences of 1s in A
    starts = find(diff([0 a]) == 1);
    ends   = find(diff([a 0]) == -1) - 1;

    for k = 1:length(starts)
        idx = starts(k):ends(k);   % the contiguous segment in A

        % If ANY 1 appears in B at these positions → set segment to 1
        if any(b(idx))
            C(r, idx) = 1;
        end
    end
end


% Assume binary matrices A, B, C of the same size
[m, n] = size(A);

% Define RGB colors
color1 = [0.0 0.8 0.2];       % onlyAC (red)
color2 = [.9 .9 .9];     % onlyAnoC (dark red)
color3 =  [1.0 0.3 1.0];       % B==1 (pink)
color4 = [1 1 1];       % noneAB (white) (b is play, A is AHMM, C is engaged)

% Initialize RGB image
RGB = ones(m, n, 3);  % default white

% Logical masks
onlyAC    = (A==1)   & (C==1);
onlyAnoC  = (A==1)   & ~(C==1);
B_is_1    = (B==1);          % replaces previous bothAB
noneAB    = ~(A );
% 1) onlyAC -> color1
for k = 1:3
    RGB(:,:,k) = RGB(:,:,k) .* ~onlyAC + color1(k) * onlyAC;
end

% 2) onlyAnoC -> color2
for k = 1:3
    RGB(:,:,k) = RGB(:,:,k) .* ~onlyAnoC + color2(k) * onlyAnoC;
end

% 3) B == 1 -> color3
for k = 1:3
    RGB(:,:,k) = RGB(:,:,k) .* ~B_is_1 + color3(k) * B_is_1;
end

% 4) noneAB -> color4 (white, but included for completeness)
for k = 1:3
    RGB(:,:,k) = RGB(:,:,k) .* ~noneAB + color4(k) * noneAB;
end

% Display
figure; imagesc(RGB);
axis xy; axis tight; axis equal;

% Plot
figure('units','normalized','outerposition',[0 0 .5 1]);
[hmm_length_ordered, pb_order] = sort(diff(filled_hmm_states'));
image(psth_edges,1:size(RGB,1), RGB(pb_order,:,:));
axis xy
xlim([-5 10])


%% Fig S2 (not shwon any more) Estiamte cumulative engagned and unegnaged lenght distributions
hmm_length_edges = 0:.1:15;
hmm_length_edges_centers = .5*(hmm_length_edges(2:end) + hmm_length_edges(1:end-1));


distribution_lengths_play         =nan(numel(prediction_struct_files),numel(hmm_length_edges)-1);
distribution_lengths_without_play =nan(numel(prediction_struct_files),numel(hmm_length_edges)-1);

cum_distribution_lengths_play         =nan(numel(prediction_struct_files),numel(hmm_length_edges)-1);
cum_distribution_lengths_without_play =nan(numel(prediction_struct_files),numel(hmm_length_edges)-1);

hmm_lenghts = diff(filled_play_bouts')';

for sn=1:numel(prediction_struct_files)
        

hmm_lengths_with_play       = hmm_lenghts(sum(is_this_hmm & is_there_play_beh,2)>0 & session_index==sn);
hmm_lengths_without_play    = hmm_lenghts(sum(is_this_hmm & is_there_play_beh,2)<=0 & session_index==sn);
n_events = sum(session_index==sn);
distribution_lengths_play(sn,:) = movmean((histcounts(hmm_lengths_with_play,hmm_length_edges))/numel(hmm_lengths_with_play),10);
distribution_lengths_without_play(sn,:) = movmean((histcounts(hmm_lengths_without_play,hmm_length_edges))/numel(hmm_lengths_without_play),10);

cum_distribution_lengths_play(sn,:) = movmean(cumsum(histcounts(hmm_lengths_with_play,hmm_length_edges))/numel(hmm_lengths_with_play),10);
cum_distribution_lengths_without_play(sn,:) = movmean(cumsum(histcounts(hmm_lengths_without_play,hmm_length_edges))/numel(hmm_lengths_without_play),10);


end

% plot engaged and unengaged state length 
line_width = 2.5;
figure
subplot(1,2,1)
plot(hmm_length_edges_centers,distribution_lengths_without_play, ':k')
hold on
plot(hmm_length_edges_centers,mean(distribution_lengths_without_play), 'k', 'LineWidth', line_width)
plot(hmm_length_edges_centers,distribution_lengths_play, ':r')
plot(hmm_length_edges_centers,mean(distribution_lengths_play), 'r', 'LineWidth', line_width)

subplot(1,2,2)
plot(hmm_length_edges_centers,cum_distribution_lengths_without_play, ':k')
hold on
plot(hmm_length_edges_centers,mean(cum_distribution_lengths_without_play), 'k', 'LineWidth', line_width)
plot(hmm_length_edges_centers,cum_distribution_lengths_play, ':r')
plot(hmm_length_edges_centers,mean(cum_distribution_lengths_play), 'r', 'LineWidth', line_width)


%% Fig S2F load varaible onset to hmm and plot time warpped data
% It will load and time wrapp all variables, but will only plot the ones in var_list
% if you also select plot_bool = true, you will also get plot for
% each variable (false by default, since we don tuse these figures)

plot_bool = false; %change this to true if you want to plot not timewraped variables


variable_onset_struct_files = dir([last_hmm_dir,'\* variable_onset_struct*']);
var_list = {'AnimalSpeed','NumCalls','RelativeDistance'};

variable_names = [];
total_number_of_hmm         = 0;
total_number_of_hmm_3states = 0;

for fn= 1:numel(variable_onset_struct_files)
     load([last_hmm_dir,'\',variable_onset_struct_files(fn).name]) 
        
     if fn==1
        
         variable_names = variable_onset_struct.variable_types';
     end
     total_number_of_hmm = total_number_of_hmm+size(variable_onset_struct.beh_properties_onset,2);
     total_number_of_hmm_3states = total_number_of_hmm_3states+size(variable_onset_struct.beh_properties_onset_3states  ,3);
end

all_variable_onsets= nan(numel(variable_names),total_number_of_hmm,size(variable_onset_struct.beh_properties_onset,3));
all_variable_offsets= nan(numel(variable_names),total_number_of_hmm,size(variable_onset_struct.beh_properties_onset,3));


all_hmm_lengths = [];
total_number_of_hmm         = 0;

for fn= 1:numel(variable_onset_struct_files)
      load([last_hmm_dir,'\',variable_onset_struct_files(fn).name]) 

        n_hmm = size(variable_onset_struct.beh_properties_onset,2);
      
        all_hmm_lengths = [all_hmm_lengths; diff(variable_onset_struct.filled_play_bouts')'];     
        variables_this_session = variable_onset_struct.beh_properties_onset ;
        for vn=1:size(variables_this_session,1)
            this_var = variables_this_session(vn,:,:);
            this_var = this_var(:);
            variables_this_session(vn,:,:) = (variables_this_session(vn,:,:) - mean(this_var, 'omitmissing'))/std(this_var, 'omitmissing');
        end
        all_variable_onsets(:,total_number_of_hmm+1:total_number_of_hmm+n_hmm,:)  =variables_this_session  ;
         
        variables_this_session = variable_onset_struct.beh_properties_offset ;
        for vn=1:size(variables_this_session,1)
            this_var = variables_this_session(vn,:,:);
            this_var = this_var(:);
           variables_this_session(vn,:,:) = (variables_this_session(vn,:,:)- mean(this_var, 'omitmissing'))/std(this_var, 'omitmissing');
        end
        all_variable_offsets(:,total_number_of_hmm+1:total_number_of_hmm+n_hmm,:)  = variables_this_session;
      
        total_number_of_hmm = total_number_of_hmm+n_hmm
end
psth_edges = variable_onset_struct.psth_edges;
% select state with and witouh play and estiamte their lengths (needed to match length of engnaged and unegnaged)

is_there_play  =  any(is_this_hmm & is_there_play_beh,2);
there_is_no_play =~any(is_this_hmm & is_there_play_beh,2);

play_lengths    = all_hmm_lengths(is_there_play);
noplay_lengths  = all_hmm_lengths(there_is_no_play);

Cost = (play_lengths-noplay_lengths').^2;

% find paried play and non pay states wih same length distributiomn (load if exist, estiamte otherwise)

if exist([last_hmm_dir,'\matching_lengths.mat'], 'file')~=2
    disp('Estimating')
    costUnmatched =4;
    figure
    p = 0;
    while p<0.05
        M = matchpairs(Cost,costUnmatched); %first colum is withplay second column witohuhtplay
        plot(play_lengths( M(:,1)), noplay_lengths( M(:,2)), '.')
        pause(.1)
        [h,p]= kstest2(play_lengths( M(:,1)), noplay_lengths( M(:,2)));
        costUnmatched = costUnmatched/1.05;
    end

    save([last_hmm_dir,'\matching_lengths.mat'], 'play_lengths','noplay_lengths','M','is_there_play','there_is_no_play')
else
    disp('Loading')
    load([last_hmm_dir,'\matching_lengths.mat'], 'play_lengths','noplay_lengths','M','is_there_play','there_is_no_play')
end
% defined indexes of states with matched lengths

is_there_play = find(is_there_play);
is_there_play = is_there_play(M(:,1));


there_is_no_play = find(there_is_no_play);
there_is_no_play = there_is_no_play(M(:,2));


original_play = any(is_this_hmm & is_there_play_beh,2);
aux_false = false(size(original_play));
aux_false(is_there_play) = 1;
is_there_play = aux_false==1;

aux_false = false(size(original_play));
aux_false(there_is_no_play) = 1;
there_is_no_play = aux_false==1;

original_play(is_there_play)
original_play(there_is_no_play)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot (select 'plot_bool=1' if needed) and time wrap %%%%
%%%% varaibles aligned to hmm onset                      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_lim = [-2 4];
is_there_play  = is_there_play & all_hmm_lengths>0;
there_is_no_play = there_is_no_play & all_hmm_lengths>0;

[hmm_length_ordered, pb_order] = sort(all_hmm_lengths);

bin_size = mean(diff(psth_edges));
time_before_after = [-10 10];

wrapped_bins_porp = 1/3; %proportion of total amount of bins
Total       =  round(range(time_before_after)/((1-wrapped_bins_porp)*bin_size));
wrapped_n   = Total -round((range(time_before_after)/bin_size)) ;



time_wraped_varaibles_play = nan(size(all_variable_onsets,1), sum(is_there_play), Total);
time_wraped_varaibles_noplay = nan(size(all_variable_onsets,1), sum(there_is_no_play), Total);
time_wrapped_time = [-(time_before_after(1):bin_size:-bin_size)/time_before_after(1), linspace(0,1,wrapped_n ), 1+(bin_size:bin_size:time_before_after(2))/time_before_after(2)];

for variable_n=1:numel(variable_names)

   
    play_lengths = hmm_length_ordered(is_there_play(pb_order));
    no_play_lengths = hmm_length_ordered(there_is_no_play(pb_order));
    matrix2plot_withplay = squeeze(all_variable_onsets(variable_n,pb_order(is_there_play(pb_order)),:));
    matrix2plot_withoutplay = squeeze(all_variable_onsets(variable_n,pb_order(there_is_no_play(pb_order)),:));
     for j=1:size(matrix2plot_withplay,1)
         data_before = psth_edges>=time_before_after(1) & psth_edges<0;
         data_during = psth_edges>=0 & psth_edges<=play_lengths(j);
         data_after = psth_edges>play_lengths(j) & psth_edges<=play_lengths(j)+time_before_after(2);


         data_before = matrix2plot_withplay(j,data_before);
         time_during = psth_edges(data_during);
         data_during = matrix2plot_withplay(j,data_during);
         data_during = interp1(time_during, data_during,linspace(time_during(1), time_during(end), wrapped_n));
         data_after = matrix2plot_withplay(j,data_after);

         if numel(data_after)<round(time_before_after(2)/bin_size)
            data_after = [data_after, nan(1, round(time_before_after(2)/bin_size)-numel(data_after))];
         end

         time_wraped_varaibles_play(variable_n,j,:) = [data_before,data_during,data_after];
     end


      for j=1:size(matrix2plot_withoutplay,1)
         data_before = psth_edges>=time_before_after(1) & psth_edges<0;
         data_during = psth_edges>=0 & psth_edges<=play_lengths(j);
         data_after = psth_edges>play_lengths(j) & psth_edges<=play_lengths(j)+time_before_after(2);


         data_before = matrix2plot_withoutplay(j,data_before);
         time_during = psth_edges(data_during);
         data_during = matrix2plot_withoutplay(j,data_during);
         data_during = interp1(time_during, data_during,linspace(time_during(1), time_during(end), wrapped_n));
         data_after = matrix2plot_withoutplay(j,data_after);

         if numel(data_after)<round(time_before_after(2)/bin_size)
            data_after = [data_after, nan(1, round(time_before_after(2)/bin_size)-numel(data_after))];
         end

         time_wraped_varaibles_noplay(variable_n,j,:) = [data_before,data_during,data_after];
     end



  
    if plot_bool
        figure('units','normalized','outerposition',[0 0 .5 1]);
        colormap(1-gray)
        subplot(5,2,1)


        imagesc(psth_edges,1:numel(play_lengths), matrix2plot_withplay)
        hold on
        plot([0 0],[1 numel(play_lengths)],'r')
        hold on
        plot(play_lengths,1:numel(play_lengths),'r')
        axis xy
        yticks([])
        ylabel('With Play')
        xlim(x_lim)
        title(variable_names{variable_n})


        subplot(5,2,3)

        imagesc(psth_edges,1:numel(no_play_lengths), matrix2plot_withoutplay)
        hold on
        plot([0 0],[1 numel(no_play_lengths)],'r')
        hold on
        plot(no_play_lengths,1:numel(no_play_lengths),'r')
        axis xy
        xlim(x_lim)
        yticks([])
        ylabel('Without Play')


        subplot(5,2,5)
        mean2plot = mean(matrix2plot_withoutplay);
        [~, ~, ci] = ttest(matrix2plot_withoutplay);
        fill([psth_edges fliplr(psth_edges )], [ci(1,:) fliplr(ci(2,:)) ], 'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
        hold on
        plot(psth_edges,mean2plot, 'k' )
        mean2plot = mean(matrix2plot_withplay);
        [~, ~, ci] = ttest(matrix2plot_withplay);
        fill([psth_edges fliplr(psth_edges )], [ci(1,:) fliplr(ci(2,:)) ], 'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
        hold on
        plot(psth_edges,mean2plot, 'r' )
        xlim(x_lim)
        ylim tight
        y_lim = ylim;
        hold on
        plot([0 0],y_lim,'b', 'HandleVisibility','off')




        subplot(5,2,[7 9])
        mean2plot      = mean(matrix2plot_withoutplay);
        [~, ~, ci] = ttest(matrix2plot_withoutplay);
        for j=find(psth_edges>=0)
            lengths_to_include = no_play_lengths>=psth_edges(j);
            mean2plot(j) = mean(matrix2plot_withoutplay(lengths_to_include,j));
            ci(:,j) = mean2plot(j) + 1.96*std(matrix2plot_withoutplay(lengths_to_include,j))*[-1 1]/sqrt(sum(lengths_to_include));
        end
        no_nan = ~any(isnan(ci));
        fill([psth_edges(no_nan) fliplr(psth_edges(no_nan) )], [ci(1,no_nan) fliplr(ci(2,no_nan)) ], 'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
        hold on
        plot(psth_edges,mean2plot, 'k' )

        mean2plot      = mean(matrix2plot_withplay);
        [~, ~, ci] = ttest(matrix2plot_withplay);
        for j=find(psth_edges>=0)
            lengths_to_include = play_lengths>=psth_edges(j);
            mean2plot(j) = mean(matrix2plot_withplay(lengths_to_include,j));
            ci(:,j) = mean2plot(j) + 1.96*std(matrix2plot_withplay(lengths_to_include,j))*[-1 1]/sqrt(sum(lengths_to_include));
        end
        no_nan = ~any(isnan(ci));
        fill([psth_edges(no_nan) fliplr(psth_edges(no_nan) )], [ci(1,no_nan) fliplr(ci(2,no_nan)) ], 'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
        hold on
        plot(psth_edges,mean2plot, 'r' )
        y_lim2 = y_lim;
        y_lim2(2) = 1.5*y_lim(2);
        hold on
        plot([0 0],y_lim2,'b', 'HandleVisibility','off')
        axis tight
        legend({ 'Without Play','Play'})
        xlim(x_lim)



        subplot(5,2,2)
        matrix2plot_withplay = squeeze(all_variable_offsets(variable_n,pb_order(is_there_play(pb_order)),:));
        imagesc(psth_edges,1:numel(play_lengths), matrix2plot_withplay)
        hold on
        plot([0 0],[1 numel(play_lengths)],'r')
        hold on
        plot(-play_lengths,1:numel(play_lengths),'r')
        axis xy
        xlim(x_lim)
        yticks([])
        ylabel('With Play')
        title(variable_names{variable_n})

        subplot(5,2,4)
        matrix2plot_withoutplay = squeeze(all_variable_offsets(variable_n,pb_order(there_is_no_play(pb_order)),:));
        imagesc(psth_edges,1:numel(no_play_lengths), matrix2plot_withoutplay)
        hold on
        plot([0 0],[1 numel(no_play_lengths)],'r')
        hold on
        plot(-no_play_lengths,1:numel(no_play_lengths),'r')
        axis xy
        xlim(x_lim)
        yticks([])
        ylabel('Without Play')


        subplot(5,2,6)
        mean2plot = mean(matrix2plot_withoutplay);
        [~, ~, ci] = ttest(matrix2plot_withoutplay);
        fill([psth_edges fliplr(psth_edges )], [ci(1,:) fliplr(ci(2,:)) ], 'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
        hold on
        plot(psth_edges,mean2plot, 'k' )
        mean2plot = mean(matrix2plot_withplay);
        [~, ~, ci] = ttest(matrix2plot_withplay);
        fill([psth_edges fliplr(psth_edges )], [ci(1,:) fliplr(ci(2,:)) ], 'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
        hold on
        plot(psth_edges,mean2plot, 'r' )

        y_lim = ylim;
        hold on
        plot([0 0],y_lim,'b', 'HandleVisibility','off')
        axis tight
        xlim(x_lim)

        subplot(5,2,[8 10])
        mean2plot      = mean(matrix2plot_withoutplay);
        [~, ~, ci] = ttest(matrix2plot_withoutplay);
        for j=find(psth_edges<=0)
            lengths_to_include = no_play_lengths>=-psth_edges(j);
            mean2plot(j) = mean(matrix2plot_withoutplay(lengths_to_include,j));
            ci(:,j) = mean2plot(j) + 1.96*std(matrix2plot_withoutplay(lengths_to_include,j))*[-1 1]/sqrt(sum(lengths_to_include));
        end
        no_nan = ~any(isnan(ci));
        fill([psth_edges(no_nan) fliplr(psth_edges(no_nan) )], [ci(1,no_nan) fliplr(ci(2,no_nan)) ], 'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
        hold on
        plot(psth_edges,mean2plot, 'k' )

        mean2plot      = mean(matrix2plot_withplay);
        [~, ~, ci] = ttest(matrix2plot_withplay);
        for j=find(psth_edges<=0)
            lengths_to_include = play_lengths>=-psth_edges(j);
            mean2plot(j) = mean(matrix2plot_withplay(lengths_to_include,j));
            ci(:,j) = mean2plot(j) + 1.96*std(matrix2plot_withplay(lengths_to_include,j))*[-1 1]/sqrt(sum(lengths_to_include));
        end
        no_nan = ~any(isnan(ci));
        fill([psth_edges(no_nan) fliplr(psth_edges(no_nan) )], [ci(1,no_nan) fliplr(ci(2,no_nan)) ], 'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
        hold on
        plot(psth_edges,mean2plot, 'r' )
        y_lim2 = y_lim;
        y_lim2(2) = 1.5*y_lim(2);
        plot([0 0],y_lim2,'b', 'HandleVisibility','off')
        legend({ 'Without Play','Play'})
        axis tight
        xlim(x_lim)


        pause(.1)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot time wrapped data (PANEL H-J) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

variables2plot = find(ismember(variable_names,var_list))';
% time_wraped_varaibles_play 
% time_wraped_varaibles_noplay 
for vn=variables2plot
    figure('units','normalized','outerposition',[0 0 .2 1]);
    subplot(2,1,1)
    colormap(1-gray)
    imagesc(time_wrapped_time,1:(size(time_wraped_varaibles_play,2)*2), [squeeze(time_wraped_varaibles_noplay(vn,:,:));squeeze(time_wraped_varaibles_play(vn,:,:))] )
    axis xy
    hold on
    plot([time_wrapped_time(1) time_wrapped_time(end)], size(time_wraped_varaibles_play,2)*[1 1], 'r')
    title(variable_names{vn})
    subplot(2,1,2)
    hold on
    matrix2plot      = squeeze(time_wraped_varaibles_noplay(vn,:,:));
    [~, ~, ci] = ttest(matrix2plot);
    no_nan = ~any(isnan(ci));
    fill([time_wrapped_time(no_nan) fliplr(time_wrapped_time(no_nan) )], [ci(1,no_nan) fliplr(ci(2,no_nan)) ], 'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
    plot(time_wrapped_time,mean(matrix2plot), 'k' )

    matrix2plot      = squeeze(time_wraped_varaibles_play(vn,:,:));
    [~, ~, ci] = ttest(matrix2plot);
    no_nan = ~any(isnan(ci));
    fill([time_wrapped_time(no_nan) fliplr(time_wrapped_time(no_nan) )], [ci(1,no_nan) fliplr(ci(2,no_nan)) ], 'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
    plot(time_wrapped_time,mean(matrix2plot), 'r' )


    % print(gcf,'-vector','-dsvg',[figure_dir,'\', variable_names{vn}, ' between states tw.svg'])
    pause(.1)
end

%% Fig S2H

% Make sure time_val is a row vector
time_val = psth_edges(:)'; % Convert to 1 x n_timepoints if needed
% Duration bounds
% Parameters
T1 = -5;  % Pre-event duration (e.g., seconds)
T2 =  5;  % Post-event duration
hmm_lengths = diff(filled_hmm_states');
val2use = hmm_lengths<max(psth_edges)-T2 & hmm_lengths>0;
hmm_lengths = hmm_lengths(val2use);
play_matrix = is_there_play_beh(val2use,:);

n_bins = 50;

% Number of trials/samples
n_trials = size(play_matrix, 1);
n_timepoints = size(play_matrix, 2);

% Make sure time_val is a row vector matching columns of is_there_play_beh
time_val = time_val(:)';

% Preallocate output matrix: rows = trials, cols = 60 interpolated time points
mz = zeros(n_trials, n_bins * 3);
last_play_latency = nan(numel(n_trials),1);
next_play_latency = nan(numel(n_trials),1);

for trial_i = 1:n_trials
    event_length = hmm_lengths(trial_i);
    play_starts   = psth_edges(find(diff(play_matrix(trial_i, :))==1));
    play_ends    =  psth_edges(find(diff(play_matrix(trial_i, :))==-1));
   

    last_Event      = max(play_ends(play_ends<0));
    following_Event = min(play_starts(play_starts>event_length))-event_length;
    if isempty(last_Event)        
        last_Event  = NaN;
    end
    if isempty(following_Event)        
        following_Event  = NaN;
    end

    last_play_latency(trial_i) = last_Event;
    next_play_latency(trial_i) = following_Event;

    % Define the three intervals for this trial
    pre_interval = [T1, 0];
    event_interval = [0, event_length];
    post_interval = [event_length, event_length + T2];

    % Extract indices and data for pre-event
    idx_pre = find(time_val >= pre_interval(1) & time_val <= pre_interval(2));
    t_pre_orig = time_val(idx_pre);
    signal_pre = play_matrix(trial_i, idx_pre);
    t_pre_interp = linspace(pre_interval(1), pre_interval(2), n_bins);
    interp_pre = interp1(t_pre_orig, signal_pre, t_pre_interp, 'linear', 'extrap');

    % Extract indices and data for event
    idx_event = find(time_val >= event_interval(1) & time_val <= event_interval(2));
    t_event_orig = time_val(idx_event);
    signal_event = play_matrix(trial_i, idx_event);
    t_event_interp = linspace(event_interval(1), event_interval(2), n_bins);
    interp_event = interp1(t_event_orig, signal_event, t_event_interp, 'linear', 'extrap');

    % Extract indices and data for post-event
    idx_post = find(time_val >= post_interval(1) & time_val <= post_interval(2));
    t_post_orig = time_val(idx_post);
    signal_post = play_matrix(trial_i, idx_post);
    t_post_interp = linspace(post_interval(1), post_interval(2), n_bins);
    interp_post = interp1(t_post_orig, signal_post, t_post_interp, 'linear', 'extrap');

    % Concatenate interpolated parts (each length 20)
    mz(trial_i, :) = [interp_pre, interp_event, interp_post];
end

[~, lenght_order] = sort(hmm_lengths);




mz_play = sum(mz(:, (n_bins+1):(2*n_bins)),2)>0;
figure
n = 256; % Number of colors in the colormap

% Create a colormap that starts at white [1 1 1] and ends at red [1 0 0]
redColormap = [linspace(1,.5,n)' linspace(1,0,n)' linspace(1,0,n)'];

% Apply the colormap
colormap(redColormap);

mz(mz>1) = 1;
hmm_with_play = ceil(mz(mz_play,:));
hmm_length_with_play = hmm_lengths(mz_play);
[~,lenght_order_with_play] = sort(hmm_length_with_play);
session_index_play          = session_index(mz_play);

hmm_without_play        = ceil(mz(~mz_play,:));
hmm_length_without_play = hmm_lengths(~mz_play);
[~,lenght_order_without_play] = sort(hmm_length_without_play);
session_index_without_play      = session_index(~mz_play); 

subplot(5,1,1:2)
play_color_array = repmat(ceil([hmm_without_play(lenght_order_without_play,:);hmm_with_play(lenght_order_with_play,:)]),1,1,3);
for k=1:3
    play_color_array(:,:,k) =  play_color_array(:,:,k) *color3(k);
end

for j=1:size(play_color_array,1)
    for   k=1:size(play_color_array,2)
        if all(play_color_array(j,k,:)==0)
            play_color_array(j,k,:)=1;
        end
    end
end
imagesc(1:size(hmm_with_play,2), 1:size(mz,1), play_color_array)
clim([0 1])
hold on
xticks([0:3]*n_bins)
plot([0 n_bins*3], size(hmm_without_play,1)*[1 1], 'r', 'LineWidth',2)
yticks([size(hmm_without_play,1)/2 size(hmm_without_play,1)+size(hmm_with_play,1)/2])
ptcg_wo_play    = round(100*size(hmm_without_play,1)/size(mz,1));
ptcg_with_play  = round(100*size(hmm_with_play,1)/size(mz,1));
yticklabels({num2str(ptcg_wo_play), num2str(ptcg_with_play)})
xticklabels([])
plot([n_bins n_bins], [1 size(mz,1)], 'c')
plot([2*n_bins 2*n_bins], [1 size(mz,1)], 'c')
set(gca, 'TickDir', 'out')

axis xy

ax = subplot(5,1,3:5)

session_list = unique(session_index_play)';
staked_mean_play = nan(numel(session_list),size(mz,2));

for sn = session_list
    staked_mean_play(sn,:) = mean(hmm_with_play(session_index_play==sn,:));
end

plot(1:size(hmm_with_play,2),100*staked_mean_play, ':k')
hold on
plot(1:size(hmm_with_play,2),100*mean(staked_mean_play), 'k', 'LineWidth',2)
hold on


session_list = unique(session_index_without_play)';
staked_mean_without_play = nan(numel(session_list),size(mz,2));
for sn = session_list
    staked_mean_without_play(sn,:) = mean(hmm_without_play(session_index_without_play==sn,:));
end
% second_red = [153 51 51]/255;
second_red = 'c';
plot(1:size(hmm_without_play,2), 100*staked_mean_without_play, ':', 'Color', second_red)
hold on
plot(1:size(hmm_with_play,2),100*mean(staked_mean_without_play), 'Color',second_red, 'LineWidth',2)
set(ax, 'TickDir', 'out')

%% Fig S2D load behavior aligned to hmm
% first load, then time warp, then plot

behavior_onset_offset_struct_files = dir([last_hmm_dir, '\*behavior_onset_offset_struct*']);

behavior_type_list = [];
total_number_of_hmm = 0;
total_number_of_hmm_3states = 0;

for fn= 1:numel(behavior_onset_offset_struct_files)
    load([last_hmm_dir, '\',behavior_onset_offset_struct_files(fn).name])
    behavior_type_list = [behavior_type_list; behavior_onset_offset_struct.behavior_tpes];
    total_number_of_hmm = total_number_of_hmm+size(behavior_onset_offset_struct.filled_play_bouts,1);
    total_number_of_hmm_3states = total_number_of_hmm_3states++size(behavior_onset_offset_struct.filled_hmm_3states,1);

end
behavior_type_list = unique(behavior_type_list);
behavior_type_list(ismember(behavior_type_list, {'Partners session', 'SA'})) = [];


merged_behaviors_onset  = zeros(numel(behavior_type_list),total_number_of_hmm,size(behavior_onset_offset_struct.behavior_offset,3));


merged_behaviors_onset_3states  = zeros(3,numel(behavior_type_list),total_number_of_hmm_3states,size(behavior_onset_offset_struct.behavior_offset,3));
all_hmm_3states = [];
all_hmm_lengths = [];
total_number_of_hmm = 0;
total_number_of_hmm_3states = 0;
session_index_again = [];
for fn= 1:numel(behavior_onset_offset_struct_files)
      disp([ 'Loading ', behavior_onset_offset_struct_files(fn).name])
    load([last_hmm_dir, '\',behavior_onset_offset_struct_files(fn).name])
    disp('Processing')
  
    n_hmm = size(behavior_onset_offset_struct.filled_play_bouts,1);
    n_hmm_3states = size(behavior_onset_offset_struct.filled_hmm_3states,1);
    all_hmm_lengths = [all_hmm_lengths; diff(behavior_onset_offset_struct.filled_play_bouts')'];
    session_index_again = [session_index_again;ones(size(behavior_onset_offset_struct.filled_play_bouts,1),1)*fn];
    this_3_States_matrix = behavior_onset_offset_struct.filled_hmm_3states;

    for j=1:3
        this_3_States_matrix (behavior_onset_offset_struct.filled_hmm_3states(:,1)==j-1,1) = re_assignment(fn,j)-1;
    end
    all_hmm_3states = [all_hmm_3states;this_3_States_matrix];

    current_behaviors = find(ismember(behavior_type_list,behavior_onset_offset_struct.behavior_tpes));

    for beavior_present = current_behaviors'
        beh_index = ismember(behavior_onset_offset_struct.behavior_tpes,behavior_type_list(beavior_present));
        merged_behaviors_onset(beavior_present,total_number_of_hmm+1:total_number_of_hmm+n_hmm,:) = behavior_onset_offset_struct.behavior_onset(beh_index,:,:);

        merged_behaviors_onset_3states(:,beavior_present,total_number_of_hmm_3states+1:total_number_of_hmm_3states+n_hmm_3states,:) = behavior_onset_offset_struct.behavior_onset_3states(re_assignment(fn,:), beh_index,:,:);


    end
    total_number_of_hmm=total_number_of_hmm+n_hmm;
    total_number_of_hmm_3states = total_number_of_hmm_3states+n_hmm_3states;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% now make timewrap %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure time_val is a row vector
time_val = psth_edges(:)'; % Convert to 1 x n_timepoints if needed
% Duration bounds
% Parameters
T1 = -5;  % Pre-event duration (e.g., seconds)
T2 =  5;  % Post-event duration
hmm_lengths = diff(filled_hmm_states');
val2use = hmm_lengths<max(psth_edges)-T2 & hmm_lengths>0;
hmm_lengths = hmm_lengths(val2use);
play_behavior_matrix = merged_behaviors_onset(:,val2use,:);

n_bins = 50;

% Number of trials/samples
n_behaviors     = size(play_behavior_matrix, 1) ;
n_trials        = size(play_behavior_matrix, 2);
n_timepoints    = size(play_behavior_matrix, 3);

% Make sure time_val is a row vector matching columns of is_there_play_beh
time_val = time_val(:)';

% Preallocate output matrix: rows = trials, cols = 60 interpolated time points
mz_behavior = zeros(n_behaviors,n_trials, n_bins * 3);
for bn = 1:n_behaviors
    disp(bn)
    for trial_i = 1:n_trials
        event_length = hmm_lengths(trial_i);




        % Define the three intervals for this trial
        pre_interval = [T1, 0];
        event_interval = [0, event_length];
        post_interval = [event_length, event_length + T2];

        % Extract indices and data for pre-event
        idx_pre = find(time_val >= pre_interval(1) & time_val <= pre_interval(2));
        t_pre_orig = time_val(idx_pre);
        signal_pre = squeeze(play_behavior_matrix(bn,trial_i, idx_pre));
        t_pre_interp = linspace(pre_interval(1), pre_interval(2), n_bins);
        interp_pre = interp1(t_pre_orig, signal_pre, t_pre_interp, 'linear', 'extrap');

        % Extract indices and data for event
        idx_event = find(time_val >= event_interval(1) & time_val <= event_interval(2));
        t_event_orig = time_val(idx_event);
        signal_event = squeeze(play_behavior_matrix(bn,trial_i, idx_event));
        t_event_interp = linspace(event_interval(1), event_interval(2), n_bins);
        interp_event = interp1(t_event_orig, signal_event, t_event_interp, 'linear', 'extrap');

        % Extract indices and data for post-event
        idx_post = find(time_val >= post_interval(1) & time_val <= post_interval(2));
        t_post_orig = time_val(idx_post);
        signal_post = squeeze(play_behavior_matrix(bn,trial_i, idx_post));
        t_post_interp = linspace(post_interval(1), post_interval(2), n_bins);
        interp_post = interp1(t_post_orig, signal_post, t_post_interp, 'linear', 'extrap');

        % Concatenate interpolated parts (each length 20)
        mz_behavior(bn,trial_i, :) = [interp_pre, interp_event, interp_post];
    end
end


%%%%%%%%%%%%%%%%%
%%%% now plot %%%%
%%%%%%%%%%%%%%%%%%

figure
ax = subplot(1,1,1)
plot(1:size(norm_mean_prob,2),  norm_mean_prob(latency_order,:) + repmat((1:numel(latency_order))',1,size(norm_mean_prob,2)), 'r')
hold on
plot(repmat([1 3*n_bins],numel(latency_order)+1,1)', repmat((1:numel(latency_order)+1)',1,2)', ':k')
plot([n_bins n_bins],[1 numel(latency_order)+1], 'k')
plot(2*[n_bins n_bins],[1 numel(latency_order)+1], 'k')
yticks((1:numel(latency_order)) +.5)
yticklabels(sub_behavior_labes)
axis tight
set(ax, 'TickDir', 'out')
%% %% analysie time serieos of hmm


%% Fig 2F (create shufled distribution and plot)
figure('units','normalized','outerposition',[0 0 1 1]);
n_rand = 1000;
bin_size = 0.01;
time2check = 120;
n_lags                      = round(time2check/bin_size);
all_play_ac                 = nan(numel(unique(session_index)), n_lags+1);
all_upper_bound95           = nan(numel(unique(session_index)), n_lags+1);
all_upper_bound99           = nan(numel(unique(session_index)), n_lags+1);
pcrtl_play                  = nan(numel(unique(session_index)), n_lags+1);
random_peaks_distribution   = [];

for sn = unique(session_index)'

    this_session_hmm            = filled_play_bouts(session_index==sn,:);
    this_session_playful_hmm    = is_this_a_playful_hmm(session_index==sn);
    all_time                    = (floor(this_session_hmm(1,1)/bin_size)*bin_size):bin_size:(ceil(this_session_hmm(end,2)/bin_size)*bin_size);

    playful_time_series         = any(all_time>=this_session_hmm(this_session_playful_hmm,1) & all_time<=this_session_hmm(this_session_playful_hmm,2));
  
    [play_acf, lags] = autocorr(double(playful_time_series), 'NumLags',n_lags);

    all_play_ac(sn,:) = play_acf;

    shufflled_play_acf = nan(n_rand,n_lags+1);
    for nr = 1:n_rand
        re_order = randsample(size(this_session_hmm,1),size(this_session_hmm,1));

        this_session_hmm = filled_play_bouts(session_index==sn,:);
        this_session_intervals = this_session_hmm(2:end,1)-this_session_hmm(1:end-1,2);
        this_session_intervals = [this_session_intervals;mean(this_session_intervals)];
        this_session_playful_hmm = is_this_a_playful_hmm(session_index==sn);

        shuffled_hmm = this_session_hmm(re_order,:);
        shuffled_hmm_lengths = diff(shuffled_hmm')';
        shuffled_intervals = this_session_intervals(re_order);
        shuffled_playful = this_session_playful_hmm(re_order,:);

        shuffled_hmm(1,1) = this_session_hmm(1,1);
        for j=1:numel(shuffled_hmm_lengths)-1

            shuffled_hmm(j,2)   = shuffled_hmm(j,1) + shuffled_hmm_lengths(j);
            shuffled_hmm(j+1,1) = shuffled_hmm(j,1) + shuffled_hmm_lengths(j) + shuffled_intervals(j);
        end
        shuffled_hmm(end,2 ) =  shuffled_hmm(end,1) + shuffled_hmm_lengths(end);
       all_time                    = (floor(shuffled_hmm(1,1)/bin_size)*bin_size):bin_size:(ceil(shuffled_hmm(end,2)/bin_size)*bin_size);

        shuffled_playful_time_series         = any(all_time>=shuffled_hmm(shuffled_playful,1) & all_time<=shuffled_hmm(shuffled_playful,2));

        [shufflled_play_acf(nr,:), ~] = autocorr(double(shuffled_playful_time_series), 'NumLags',round(time2check/bin_size));
    end
  
    matrix2plot = shufflled_play_acf;
    this_pctl = 100 * mean(repmat(play_acf,n_rand,1)<= shufflled_play_acf);

    for nr = 1:25:n_rand

        shufled_pctl =  100 * mean(repmat(shufflled_play_acf(nr,:),n_rand,1)<= shufflled_play_acf);
        y = 100  -shufled_pctl;
        y = smoothdata(y,'gaussian',10/bin_size);
        [pks,locs,w,p] = findpeaks(y,'MinPeakHeight',80,'MinPeakProminence', .4*range(y), 'MinPeakDistance',10/bin_size);
        random_peaks_distribution =  [random_peaks_distribution; [(lags(locs)*bin_size)'  (lags(locs)*bin_size)'*0+sn (lags(locs)*bin_size)'*0+nr]];

    end
    



    pcrtl_play(sn,:) = this_pctl;
    upper_bound95 = prctile(matrix2plot, 95, 1);   % 5th percentile across rows (dim 1)
    upper_bound99 = prctile(matrix2plot, 99, 1);


    subplot(5,3, sn)
    plot(lags*bin_size,play_acf , 'r')
    hold on
    plot(lags*bin_size, upper_bound95, 'b')
    plot(lags*bin_size, upper_bound99, 'g')


    all_upper_bound95(sn,:) =upper_bound95;
    all_upper_bound99(sn,:) = upper_bound99;
    xlim([0 20])
    pause(.1)
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot autocorrelogram excess play %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

excces_play_time    = 5;
x_lim               = [0 10];
y_lim_pctl          = [75 100];
y_lim_p_val         = [.9 1];

figure
colormap(jet)
subplot(5 ,2,[1 3 5])
imagesc(lags*bin_size, 1:15,all_play_ac-all_upper_bound95)
axis xy
clim([-.2 .2])
xlim(x_lim)

rank_signficiant = nan(1, size(all_play_ac,2));
for t=1:numel(rank_signficiant)
rank_signficiant(t) = signrank(all_play_ac(:,t),all_upper_bound95(:,t), 'tail','right');
end

subplot(5,2,[7 9])
plot(x_lim, [0 0], 'k')
xlim(x_lim)
hold on
plot(lags*bin_size,all_play_ac-all_upper_bound95, ':k')
plot(lags*bin_size,mean(all_play_ac-all_upper_bound95), 'k', 'LineWidth',3)
xlim(x_lim)
ylim([-.2 .2])

yyaxis right
semilogy(lags*bin_size,1-rank_signficiant, 'Color',  [0.5, 0, 0.5])
hold on
semilogy(lags*bin_size,lags*0 + 1-0.05, ':', 'Color',  [0.5, 0, 0.5])
ylim(y_lim_p_val)

subplot(5 ,2,[1 3 5]+1)
imagesc(lags*bin_size, 1:15,100-pcrtl_play)
clim(y_lim_pctl)
axis xy
xlim(x_lim)


subplot(5,2,[7 9]+1)
plot(x_lim, [95 95], 'k')
hold on
plot(lags*bin_size,100-pcrtl_play, ':k')
plot(lags*bin_size,100-mean(pcrtl_play), 'k', 'LineWidth',3)
matrix2plot = 100-pcrtl_play;
ci = [-1; 1]*std(matrix2plot)/sqrt(size(matrix2plot,1)) + [1;1]*mean(matrix2plot);
fill([lags*bin_size fliplr(lags*bin_size)], [ci(1,:) fliplr(ci(2,:))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
xlim(x_lim)
ylim(y_lim_pctl)


%% Fig S2I find second peaks  
figure
all_time_events = [];
for sn=1:15
y = 100  -pcrtl_play(sn,:);
y = smoothdata(y,'gaussian',10/bin_size);
[pks,locs,w,p] = findpeaks(y,'MinPeakHeight',80,'MinPeakProminence', .4*range(y), 'MinPeakDistance',20/bin_size);





subplot(5,3,sn)
plot(lags*bin_size,y)
hold on
plot(lags(locs)*bin_size, y(locs), '.r')

all_time_events =[all_time_events; [(lags(locs)*bin_size)'  (lags(locs)*bin_size)'*0+sn]];
end

% plot distributiuomn of second play peak and random distribution (PLAY BOUT FIGURE C)


[real_f, real_x] = histcounts(all_time_events(:,1), 0:1:120,'Normalization','cdf');  

n_sessions = unique(random_peaks_distribution(:,2));
distribution_random = nan(numel(n_sessions), numel(real_f));
Xboots = cell(numel(n_sessions),1);
all_ks = [];

for nn =1:numel(n_sessions)
    Xboots{nn} =random_peaks_distribution(random_peaks_distribution(:,2)==nn,1); 
[random_f, random_x] = histcounts(random_peaks_distribution(random_peaks_distribution(:,2)==nn,1), 0:1:120,'Normalization','cdf'); 
distribution_random(nn,:) =random_f;
[ks,p] = kstest2(random_peaks_distribution(random_peaks_distribution(:,2)==nn,1),all_time_events(:,1));
all_ks = [all_ks;p];
end


figure
bin_centers = .5*(real_x(1:end-1) + real_x(2:end));
plot(bin_centers, real_f, 'r')
% f = CDF values, x = sorted data
hold on
plot(bin_centers, distribution_random, 'k')
xlabel('Value'); ylabel('CDF');

%% Fig S2 B

cd(hmm_raw_data)
model_comp_data= dir('model_comp_*');


staked_ll = [];
staked_aic = [];
staked_BIC = [];

for j=1:numel(model_comp_data)

    a = readNPY(model_comp_data(j).name);


staked_ll = [staked_ll,a(:,1)];
staked_aic = [staked_aic,a(:,2)];
staked_BIC = [staked_BIC,a(:,3)];
end

figure

subplot(1,3,1)
plot(2:5,staked_aic, 'k:')
xticks(2:5)
hold on
plot(2:5,mean(staked_aic'), 'k')
xlim([1.5 5.5])
xlabel('Number of States')
ylabel('AIC')
set(gca, 'FontSize',12)
set(gca,'TickDir','out');

subplot(1,3,2)
plot(2:5,staked_BIC, 'k:')
xlabel('Number of States')
ylabel('BIC')
hold on
plot(2:5,mean(staked_BIC'), 'k')
set(gca, 'FontSize',12)
xlim([1.5 5.5])
xticks(2:5)
set(gca,'TickDir','out');


subplot(1,3,3)
plot(2:5,staked_ll, 'k:')
hold on
plot(2:5,mean(staked_ll'), 'k')
xlabel('Number of States')
ylabel('LL')
set(gca, 'FontSize',12)
xlim([1.5 5.5])
xticks(2:5)
set(gca,'TickDir','out');



%% Fig 2E
last_hmm_dir = analyssis_folder;


play_behavior_struct_files = dir([last_hmm_dir, '\*play_behavior_struct*']);


all_behavior_types = [];

for fn= 1:numel(play_behavior_struct_files)
     load([last_hmm_dir, '\',play_behavior_struct_files(fn).name]) 
     all_behavior_types = [all_behavior_types; play_behavior_struct.behavior_types];

end
all_behavior_types = unique(all_behavior_types);

all_numbers         = zeros(numel(all_behavior_types),1);
all_behavior_count  = zeros(numel(all_behavior_types),3);

for fn= 1:numel(play_behavior_struct_files)
    load([last_hmm_dir, '\',play_behavior_struct_files(fn).name]) 
    behavior_count = diag(play_behavior_struct.numbers)*play_behavior_struct.proportions;

    for bn = find(ismember(all_behavior_types,play_behavior_struct.behavior_types))'
        beh_index = ismember(play_behavior_struct.behavior_types, all_behavior_types(bn));
        all_numbers(bn) = all_numbers(bn)+play_behavior_struct.numbers(beh_index);

        all_behavior_count(bn,:) =  all_behavior_count(bn,:)+behavior_count(beh_index,:);
    end
end

behaviors2merge = {'Pounce_A','Pounce_B'};
re_name = {'Pounce'};
value = sum(all_behavior_count(ismember(all_behavior_types,behaviors2merge),:));
all_behavior_count(ismember(all_behavior_types,behaviors2merge),:) = repmat(value,numel(behaviors2merge), 1 );
value = sum(all_numbers(ismember(all_behavior_types,behaviors2merge),:));
all_numbers(ismember(all_behavior_types,behaviors2merge)) = value;


all_behavior_count(ismember(all_behavior_types,behaviors2merge(1)),:) = [];

all_numbers(ismember(all_behavior_types,behaviors2merge(1))) = [];
all_behavior_types(ismember(all_behavior_types,behaviors2merge(1))) = [];
all_behavior_types(ismember(all_behavior_types,behaviors2merge)) = re_name;

size(all_numbers)
size(all_behavior_count)
size(all_behavior_types)

behaviors2merge = {'Pounce_Ai','Pounce_Bi'};
re_name = {'PounceI'};
value = sum(all_behavior_count(ismember(all_behavior_types,behaviors2merge),:));
all_behavior_count(ismember(all_behavior_types,behaviors2merge),:) = repmat(value,numel(behaviors2merge), 1 );
value = sum(all_numbers(ismember(all_behavior_types,behaviors2merge),:));
all_numbers(ismember(all_behavior_types,behaviors2merge)) = value;

all_behavior_count(ismember(all_behavior_types,behaviors2merge(1)),:) = [];
all_numbers(ismember(all_behavior_types,behaviors2merge(1))) = [];
all_behavior_types(ismember(all_behavior_types,behaviors2merge(1))) = [];

all_behavior_types(ismember(all_behavior_types,behaviors2merge)) = re_name;


behaviors2delete = {'','Sniffing_C'};

all_behavior_count(ismember(all_behavior_types,behaviors2delete),:) = [];
all_numbers(ismember(all_behavior_types,behaviors2delete)) = [];
all_behavior_types(ismember(all_behavior_types,behaviors2delete)) = [];


all_behavior_proportions = diag(1./all_numbers)*all_behavior_count;

[~,bar_order] = sort(all_behavior_proportions(:,2), 'descend');
figure

toincludemove = ~ismember(all_behavior_types(bar_order), {'Partners session','Sniffing_C',''});
subplot(5,1,1)
bar(all_numbers( bar_order(toincludemove)))
xticks(1:numel(bar_order(toincludemove)))
xticklabels([])
axis tight


subplot(5,1,2:4)
bar(all_behavior_proportions(bar_order(toincludemove),[2,1,3]), 'stacked')
xticks(1:numel(bar_order(toincludemove)))


all_behavior_types_realNames = {'Unlabeled','Bite','Boxing','PlayfulApproach','Chasing','nonPlayfulApproach','Escape','Evasion','Grooming',...
    'Pin','PounceNeck','PounceNeckImmob','PounceBack','PounceBackImmob','Rearing','Scratch','Sniffing'};
% all_behavior_types_realNames= all_behavior_types
xticklabels(all_behavior_types(bar_order(toincludemove)))
axis tight
%%