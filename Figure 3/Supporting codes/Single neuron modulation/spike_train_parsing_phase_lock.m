% Session-level generator (same metrics, GENERATE_PHASE_COUPLING folder loading, no plots):
%   GENERATE_SPIKE_TRAIN_PARSING.m
audio_alignment = true;
%%

% load('LFP_PAG')
% LFP = LFP_PAG;
% clear LFP_PAG

load('LFP')


sound(sin(pi*(440*(2^(5/12)))*(1:200000)/200000)', 200000)
%%
sampling_freq = 2500;
% chanel_map_name = 'ChanelMap_regular.mat';
% load(chanel_map_name)
% 
% behavior_files = dir('*.txt');
% Behavior =  readtable('_NMMTT_230929 23-10-12 13-28-07.txt');
% Behavior =  readtable(behavior_files(1).name);
lfp_time_shared_timeline = (1:size(LFP,2))/sampling_freq;


% disp(animal_types)
%% define play session ranges
partner_names = {'Dummy2','Dummy4'};

% animal_n = ismember(animal_types, 'Single'); % 3 is single for 1110
animal_sessions = nan(numel(partner_names),2);

for j=1:2
    animal_sessions(j,:) = [min(Behavior.Start(ismember(Behavior.Animal,partner_names{j} ))) max(Behavior.End(ismember(Behavior.Animal,partner_names{j} )))];
end
[~,session_order] = sort(animal_sessions(:,1));
partner_names = partner_names(session_order);
animal_sessions = animal_sessions(session_order,:);
%%
spike_times     = readNPY('spike_times.npy');
spike_clusters  = readNPY('spike_clusters.npy');
cluster_group   = tdfread('cluster_group.tsv');
cluster_data    = tdfread('cluster_info.tsv');
clusters_list   = cluster_data.cluster_id( ismember(cluster_data.group, {'good', 'mua'}));
cluster_type    = cluster_data.group(      ismember(cluster_data.group, {'good', 'mua'}),:);

chanel_list    = cluster_data.ch(         ismember(cluster_data.group, {'good', 'mua'}));
% ycoords_all_clusters     = ycoords(chanel_list);

last_reelvant_event = double(max(spike_times))/30000;
%%
n_rand = 500;
sampling_freq = 2500;
theta_range = [6 12];
theta_filter = designfilt('bandpassfir', ...
    'FilterOrder',2500,'CutoffFrequency1',theta_range(1), ...
    'CutoffFrequency2',theta_range(2),'SampleRate',sampling_freq);
delta_range = [1 5];
delta_filter = designfilt('bandpassfir', ...
    'FilterOrder',2500,'CutoffFrequency1',delta_range(1), ...
    'CutoffFrequency2',delta_range(2),'SampleRate',sampling_freq);

bin_size = 0.001;
smoothing_time = 0.04;
smooth_window = round(smoothing_time/bin_size);
max_coverage = .9;
rate_wind       = .001;
norm_range = [-3.5 3.5];

conv_function   = 500*normpdf(linspace(norm_range(1), norm_range(2), 101))/sum(normpdf(linspace(norm_range(1), norm_range(2), 101)));


figure
plot(-.05:0.001:.05, conv_function)

% smooth_window = 1;
hist_range          = bin_size*ceil(1/(theta_range(1)*bin_size))*[-1 1];
hist_edges          = hist_range(1):bin_size:hist_range(2);
hist_edges_centers  = hist_edges(2:end)-bin_size/2;

angle_hist_range            = [-pi pi];
angle_bin_size              = pi/16;
angle_hist_edges            = angle_hist_range(1):angle_bin_size:angle_hist_range(2);
angle_hist_edges_centers    = angle_hist_edges(2:end)-angle_bin_size/2;

temporal_shift_range        = [-.250 .250];
temporal_shift_step         = 0.01;
temporal_shift_values       = temporal_shift_range(1):temporal_shift_step:temporal_shift_range(2);



fit_tresh = 0;
autocorr_range = [-.5 .5];
autocorr_bin_size = .01;
autocorr_bin_edges = autocorr_range(1):autocorr_bin_size:autocorr_range(2);
autocorr_centers = autocorr_bin_edges(2:end) - autocorr_bin_size/2;
treshold = 2;
min_run_length = .3;
min_run_spikes = 4;
n_psth_samples      = 50;

run_events_table_per_cell = cell(numel(clusters_list),1);
y_labels = {'0', 'pi/2', 'pi', '2pi/2'};
%% save parameters as config.

config = [];
config.n_rand =n_rand;
config.theta_range =theta_range;
config.delta_range = delta_range;
config.bin_size =bin_size;
config.smoothing_time =smoothing_time;
config.smooth_window =smooth_window;
config.max_coverage =max_coverage;
config.rate_wind =rate_wind;
config.norm_range =norm_range;
config.conv_function =conv_function;
config.hist_range = hist_range;
config.hist_edges = hist_edges;
config.hist_edges_centers =hist_edges_centers;
config.angle_hist_range =angle_hist_range;
config.angle_bin_size =angle_bin_size;
config.angle_hist_edges =angle_hist_edges;
config.angle_hist_edges_centers =angle_hist_edges_centers;
config.fit_tresh =fit_tresh;
config.autocorr_range =autocorr_range;
config.autocorr_bin_size =autocorr_bin_size;
config.autocorr_bin_edges =autocorr_bin_edges;
config.autocorr_centers =autocorr_centers;
config.treshold =treshold;
config.min_run_length =min_run_length;
config.min_run_spikes =min_run_spikes;
config.temporal_shift_range =temporal_shift_range;
config.temporal_shift_step =temporal_shift_step;
config.temporal_shift_values =temporal_shift_values;
config.n_psth_samples = n_psth_samples;

%%



for dn=129:numel(clusters_list)

    id =clusters_list(dn);
    ch = cluster_data.ch(cluster_data.cluster_id==id);
    ct = cluster_data.group(cluster_data.cluster_id==id,:);
    chanel_LFP = LFP(ch,:);
    spike_time_sec = double(spike_times(spike_clusters==id))/30000;
    filtered_LFP = filtfilt(theta_filter,chanel_LFP(lfp_time_shared_timeline<last_reelvant_event));
    session_lfp_time = lfp_time_shared_timeline(lfp_time_shared_timeline<last_reelvant_event);
     hilbert_transformed = hilbert(filtered_LFP);
    current_angles      = angle(hilbert_transformed);
    current_envelope    =  abs(hilbert_transformed);


    peak_prominence    = .1*std(filtered_LFP);
    peak_distance_sec  = 1/theta_range(2);
    [~, min_locs]      = findpeaks(-filtered_LFP, 'MinPeakDistance',peak_distance_sec*sampling_freq , 'MinPeakProminence', peak_prominence);


    % sound(sin(pi*(440*(2^(5/12)))*(1:200000)/200000)', 200000)
    % sound(sin(pi*(880*(2^(0/12)))*(1:200000)/200000)', 200000)
    % pause(0.1)
    % sound(sin(pi*(880*(2^(3/12)))*(1:200000)/200000)', 200000)
    % pause(0.1)
    % sound(sin(pi*(880*(2^(5/12)))*(1:200000)/200000)', 200000)
    % pause(0.1)
    % sound(sin(pi*(880*(2^(12/12)))*(1:200000)/200000)', 200000)
    % 

    [f,x,~,~] = ecdf(current_angles);

    f(find(diff(x)==0)) = [];
    x(find(diff(x)==0)) = [];

    uniform_phases = 2*pi*interp1(x,f,current_angles) - pi;


    spike_times_sec =  double(spike_times(ismember(spike_clusters , id)))/30000;
   

    this_neuron_phases = interp1(session_lfp_time, uniform_phases, spike_times_sec);
    this_neuron_phases = this_neuron_phases(~isnan(this_neuron_phases));
    this_neuron_spikes = spike_times_sec(~isnan(this_neuron_phases));


    all_rate = zeros(ceil(last_reelvant_event/rate_wind),1);
    all_rate_time = (1:numel(all_rate))*rate_wind;


    [~, spike_index] = groupcounts(ceil(this_neuron_spikes/rate_wind));
    spike_index(spike_index>numel(all_rate))=[];
    all_rate(spike_index) =1;

    all_rate = conv(all_rate', conv_function);
    all_rate = all_rate(51:end-50);

    run_events_table = find_runs(all_rate',treshold,rate_wind,this_neuron_spikes, min_run_spikes, min_run_length, max_coverage);



   
    polar_histogram                 = nan(size(run_events_table,1), numel(angle_hist_edges_centers));    
    mean_theta_alignment              = nan(size(run_events_table,1),numel(hist_edges_centers));
    real_theta_alignment            = nan(size(run_events_table,1),numel(hist_edges_centers));
    mean_raw_theta_alignment        = nan(size(run_events_table,1),numel(hist_edges_centers));
    real_raw_theta_alignment        = nan(size(run_events_table,1),numel(hist_edges_centers));
    fitted_theta_alignment          = nan(size(run_events_table,1),numel(hist_edges_centers));
    sin_fit_param                   = nan(size(run_events_table,1),6);
    sin_fit_param_no_trend          = nan(size(run_events_table,1),6);
    sin_fit_param_no_mean           = nan(size(run_events_table,1),6);
    lfp_param                       = nan(size(run_events_table,1),3);
    this_neuron_autoccorrelograms   = nan(size(run_events_table,1),numel(autocorr_bin_edges)-1);
    phase_entrainment               = nan(size(run_events_table,1),4); 
    temporal_shifts                 = nan(size(run_events_table,1), numel(temporal_shift_values));
            
  
    for j=1:size(run_events_table,1)
        run_start           = run_events_table.RunStartTime(j);
        run_end             = run_events_table.RunEndTime(j);
        this_run_peaks      = min_locs(session_lfp_time(min_locs)>=run_start & ...
            session_lfp_time(min_locs)<=run_end);
        this_run_peaks      = session_lfp_time(this_run_peaks);
        this_run_lfp        = filtered_LFP(session_lfp_time>=run_start & session_lfp_time<=run_end);
        this_run_evnelope   = current_envelope(session_lfp_time>=run_start & session_lfp_time<=run_end);
        if isempty(this_run_lfp)
            lfp_range = NaN;
            mean_envelope = NaN;
            mean_freq = NaN;

        else
            lfp_range = range(this_run_lfp)/2;
            mean_envelope= median(this_run_evnelope);
            if numel(this_run_peaks)>1
                mean_freq = 1/mean(diff(this_run_peaks));
            else
                mean_freq = NaN;
            end
        end
        lfp_param(j,:)  = [lfp_range mean_envelope mean_freq];   
        this_run_psth   = nan(numel(this_run_peaks),numel(hist_edges_centers));
    
        this_run_psth_samples   = min(n_psth_samples, 2^numel(this_run_peaks));
        mean_responses          = nan(this_run_psth_samples,numel(hist_edges_centers));
        raw_responses           = nan(this_run_psth_samples,numel(hist_edges_centers));
        fited_responses         = nan(this_run_psth_samples,numel(hist_edges_centers));

        for nps =1:this_run_psth_samples
            spikes_this_psth = this_neuron_spikes;

            order2collect = randsample(numel(this_run_peaks),numel(this_run_peaks));
            for peak_n = order2collect'
                peak_time = this_run_peaks(peak_n);
                spikes2plot = spikes_this_psth(spikes_this_psth>=peak_time+hist_range(1) & ...
                    spikes_this_psth<=peak_time+hist_range(2))-peak_time;
                spikes_this_psth(spikes_this_psth>=peak_time+hist_range(1) & ...
                    spikes_this_psth<=peak_time+hist_range(2))=[];
                this_run_psth(peak_n,:) = histcounts(spikes2plot,hist_edges);
            end
            this_run_spikes = this_neuron_spikes(this_neuron_spikes>=run_start & this_neuron_spikes<=run_end);
            all_lags = autocorrelogram(this_run_spikes,autocorr_range);
            this_neuron_autoccorrelograms(j,:) = histcounts(all_lags,autocorr_bin_edges);
            if size(this_run_psth,1)>1
                mean_response   = movmean(mean(this_run_psth),smooth_window)/bin_size;
                raw_response    = mean(this_run_psth)/bin_size;
                mean_responses(nps,:)   = mean_response;
                raw_responses(nps,:)    = raw_response;
            elseif size(this_run_psth,1)==1
                mean_response = movmean((this_run_psth),smooth_window)/bin_size;
                raw_response = this_run_psth/bin_size;
                mean_responses(nps,:)   = mean_response;
                raw_responses(nps,:)    = raw_response;
            end
           
        end
        if this_run_psth_samples>1
            mean_response   = mean(mean_responses,"omitnan");
            raw_response    = mean(raw_responses,"omitnan");
        else
            mean_response   = mean_responses;
            raw_response    = raw_responses;
        end

        
        real_index = min(find(~(any(isnan(mean_responses),2))));
        if ~isempty(real_index)
            real_theta_values = mean_responses(real_index,:);
            real_theta_alignment(j,:)     = mean_responses(real_index,:);
            real_raw_theta_alignment(j,:) = raw_responses(real_index,:);

            y = real_theta_values(ceil(.5*smooth_window):end-ceil(.5*smooth_window));
            rate_range = range(y);
            min_rate = min(y);
            if ~any(isnan(y)) & ~isempty(y)
                y = (y-min_rate)/rate_range;
                x = hist_edges_centers(ceil(.5*smooth_window):end-ceil(.5*smooth_window));
                [~,peaks_mean_positive] = findpeaks(y, 'MinPeakDistance',(1/theta_range(2))/bin_size, 'MinPeakProminence', std(y));
                [~,peaks_mean_negative] = findpeaks(-y, 'MinPeakDistance',(1/theta_range(2))/bin_size, 'MinPeakProminence', std(y));
                consecutive_peaks = sort([peaks_mean_positive,peaks_mean_negative]);
                per = 2*mean(diff(consecutive_peaks))*bin_size;
                if numel(consecutive_peaks) <2 || per==0
                    per = 1/mean(theta_range);
                end
                mdl = fittype('a*sin(b*x + c) + d*x + e','indep','x');
                [fittedmdl,gof] = fit(x',y',mdl,'start',[rand(),1/(per/(2*pi)),rand(), rand(), rand()]);


                sin_fit_param_no_mean(j,1) = fittedmdl.a*rate_range;
                sin_fit_param_no_mean(j,2) = fittedmdl.b;
                sin_fit_param_no_mean(j,3) = fittedmdl.c;
                sin_fit_param_no_mean(j,4) = fittedmdl.d*rate_range;
                sin_fit_param_no_mean(j,5) = fittedmdl.e*rate_range + min_rate;
                sin_fit_param_no_mean(j,6) = gof.rsquare;
                                
            end
        end

        y = mean_response(ceil(.5*smooth_window):end-ceil(.5*smooth_window));
        rate_range = range(y);
        min_rate = min(y);
        if ~any(isnan(y)) & ~isempty(y)
            y = (y-min_rate)/rate_range;
            x = hist_edges_centers(ceil(.5*smooth_window):end-ceil(.5*smooth_window));
            [~,peaks_mean_positive] = findpeaks(y, 'MinPeakDistance',(1/theta_range(2))/bin_size, 'MinPeakProminence', std(y));
            [~,peaks_mean_negative] = findpeaks(-y, 'MinPeakDistance',(1/theta_range(2))/bin_size, 'MinPeakProminence', std(y));
            consecutive_peaks = sort([peaks_mean_positive,peaks_mean_negative]);
            per = 2*mean(diff(consecutive_peaks))*bin_size;
            if numel(consecutive_peaks) <2 || per==0
                per = 1/mean(theta_range);
            end
            mdl = fittype('a*sin(b*x + c) + d*x + e','indep','x');
            [fittedmdl,gof] = fit(x',y',mdl,'start',[rand(),1/(per/(2*pi)),rand(), rand(), rand()]);
       

            sin_fit_param(j,1) = fittedmdl.a*rate_range;
            sin_fit_param(j,2) = fittedmdl.b;
            sin_fit_param(j,3) = fittedmdl.c;
            sin_fit_param(j,4) = fittedmdl.d*rate_range;
            sin_fit_param(j,5) = fittedmdl.e*rate_range + min_rate;
            sin_fit_param(j,6) = gof.rsquare;
            fitted_theta_alignment(j,ceil(.5*smooth_window):end-ceil(.5*smooth_window)) = ...
            fittedmdl(x)*rate_range + min_rate;

            [fittedmdl,gof] = fit_sin(x',y',[rand(),1/(per/(2*pi)),rand(), rand()]);
            
            sin_fit_param_no_trend(j,1) = fittedmdl.a*rate_range;
            sin_fit_param_no_trend(j,2) = fittedmdl.b;
            sin_fit_param_no_trend(j,3) = fittedmdl.c;
            sin_fit_param_no_trend(j,4) = fittedmdl.d*rate_range + min_rate;
            sin_fit_param_no_trend(j,5) = gof.rsquare;

        end
        mean_theta_alignment(j,ceil(.5*smooth_window):(end-ceil(.5*smooth_window))) = ...
            mean_response(ceil(.5*smooth_window):end-ceil(.5*smooth_window));
        mean_raw_theta_alignment(j,:) = raw_response;
        
        this_phases                 = this_neuron_phases(this_neuron_spikes>=run_start  & ...
            this_neuron_spikes<=run_end);


        time_shift_index= 1;
        for time_shift=temporal_shift_values
           if numel(this_run_spikes)>0 &&  ~any(floor((this_run_spikes+time_shift)*2500)>size(uniform_phases,2)) && ...
                   ~any(floor((this_run_spikes+time_shift)*2500)<=0) 
            temporal_shifts(j,time_shift_index) = circ_r(uniform_phases(floor((this_run_spikes+time_shift)*2500))' );
           end
            time_shift_index = time_shift_index+1;
        end

        r                           = circ_r(this_phases);
        phase                       = circ_mean(this_phases);
        if numel(this_run_spikes)>1
            shifted_spiketrain          = shift_spiketrain(this_run_spikes,run_start, run_end, rand(n_rand,1)*(run_end-run_start));
        else
            shifted_spiketrain = [];
        end
        if ~isempty(shifted_spiketrain) && ~any(any(floor(shifted_spiketrain*2500)>size(uniform_phases,2)))
            this_neuron_random_phases   = uniform_phases(floor(shifted_spiketrain*2500));

            rand_mvl                    = circ_r(this_neuron_random_phases');
            rand_mvl                    = sort(rand_mvl);
            [~,loc]                     = min(abs(r-rand_mvl));
            r_stat                      = mean(loc)/n_rand;
            [~,loc]                     = min(abs(max(temporal_shift_values)-rand_mvl));
            r_stat_max                  = mean(loc)/n_rand; 
            phase_entrainment(j,:)      = [r r_stat r_stat_max phase];
        end

        polar_histogram(j,:)         = histcounts(this_phases,angle_hist_edges);
    end
    % sound(sin(pi*(440*(2^(12/12)))*(1:200000)/200000)', 200000)   
    % pause(0.1)
    % sound(sin(pi*(880*(2^(5/12)))*(1:200000)/200000)', 200000)
    % pause(0.1)
    % sound(sin(pi*(880*(2^(3/12)))*(1:200000)/200000)', 200000)
    % pause(0.1)
    % sound(sin(pi*(880*(2^(0/12)))*(1:200000)/200000)', 200000)
  
    % xlabel('Time')
    % title(['ID =', num2str(id), ' CH = ' , num2str(ch), ' ',  neuorn_area{1}])
    % saveas(gcf, ['ID =', num2str(id), ' CH = ' , num2str(ch), ' ',  neuorn_area{1}, ' Firing episodes.fig'])
    % close(gcf)
    run_events_table.SinFit                 = sin_fit_param;    
    run_events_table.SinFit_NOTREND         = sin_fit_param_no_trend;    
    run_events_table.SinFit_NOMEAN          = sin_fit_param_no_mean;    
    run_events_table.ThetaPsth              = mean_theta_alignment;
    run_events_table.ThetaPsthReal          = real_theta_alignment;
    run_events_table.RowThetaPsth           = mean_raw_theta_alignment;
    run_events_table.RowThetaPsthReal       = real_raw_theta_alignment;
    run_events_table.autoccorrelograms      = this_neuron_autoccorrelograms;
    run_events_table.fitted_theta_alignment = fitted_theta_alignment;
    run_events_table.lfp_param              = lfp_param;
    run_events_table.phase_entrainment      = phase_entrainment;
    run_events_table.temporal_shifts        = temporal_shifts;
    run_events_table.Id                     = ones(size(run_events_table,1),1)*id;
    run_events_table.Ch                     = ones(size(run_events_table,1),1)*ch;
    run_events_table.ct                     = repmat(ct, size(run_events_table,1),1);


    figure('units','normalized','outerposition',[0 0 1 1])
    colormap(1-gray)
    if size(run_events_table,1)>1
        subplot(1,9,1:2)
        matrix2plot = real_raw_theta_alignment;
        phase_shift = mod(sin_fit_param(:,3), 2*pi);
        amplitud    = sin_fit_param(: ,1);
        phase_shift = phase_shift + (amplitud<0)*pi;
        phase_shift = mod(phase_shift, 2*pi);
        if size(matrix2plot,1)>1
            [sorted_phase_shift,matrix2plot_order] = sort(phase_shift);

            imagesc(hist_edges_centers, 1:(2*numel(phase_shift)), [matrix2plot(matrix2plot_order,:);matrix2plot(matrix2plot_order,:)])
            [~,pi2_loc] =min(abs(sorted_phase_shift- pi/2));
            [~,pi_loc] =min(abs(sorted_phase_shift-pi));
            [~,pi32_loc] =min(abs(sorted_phase_shift-3*pi/2));
            [~,pi0_loc] =min(abs(sorted_phase_shift));
            y_ticks = [pi0_loc pi2_loc pi_loc pi32_loc numel(phase_shift)+[pi0_loc pi2_loc pi_loc pi32_loc]];
            y_ticks = unique(y_ticks);
            y_tick_labels = [y_labels,y_labels];
            if size(matrix2plot,1)>20
                yticks(y_ticks)
                yticklabels(y_tick_labels)
            end
            axis xy
        end
        xlabel('Time (s)')
        clim([0 1])
        title(['ID =', num2str(id), ' CH = ' , num2str(ch), ' ',  ct])

        subplot(1,9,3:4)
        matrix2plot = mean_theta_alignment;

        if size(matrix2plot,1)>1
            imagesc(hist_edges_centers, 1:(2*numel(phase_shift)), [matrix2plot(matrix2plot_order,:);matrix2plot(matrix2plot_order,:)])
            if size(matrix2plot,1)>20
                yticks(y_ticks)
                yticklabels(y_tick_labels)
            end
            axis xy
        end
        xlabel('Time (s)')
        yticks([])

        subplot(1,9,5:6)
        matrix2plot = fitted_theta_alignment;
        for j=1:size(matrix2plot,1)

            matrix2plot(j,ceil(.5*smooth_window):(end-ceil(.5*smooth_window))) = ...
                zscore( matrix2plot(j,ceil(.5*smooth_window):(end-ceil(.5*smooth_window))));
        end

        if size(matrix2plot,1)>1
            imagesc(hist_edges_centers, 1:(2*numel(phase_shift)), [matrix2plot(matrix2plot_order,:);matrix2plot(matrix2plot_order,:)])
            if size(matrix2plot,1)>20
                yticks(y_ticks)
                yticklabels(y_tick_labels)
            end
            axis xy
        end
        xlabel('Time (s)')
        yticks([])




        subplot(4,9,7)
        [~,TRMF1]               = rmoutliers(abs(sin_fit_param(:,1)));
        [~,TRMF2]               = rmoutliers(abs(sin_fit_param(:,5)));
        osc_amplitudes2plot     =  abs(sin_fit_param(~TRMF1 & ~TRMF2,1));
        base_amplitudes2plot    = abs(sin_fit_param(~TRMF1 & ~TRMF2,5));

        histogram(osc_amplitudes2plot)
        title('Oscillation Rate')
        xlabel('Rate (Hz)')
        subplot(4,9,8)
        histogram(base_amplitudes2plot)
        xlabel('Rate (Hz)')
        title('Base Rate')
        subplot(4,9,9)
        amplitude_baseline_rate = osc_amplitudes2plot./base_amplitudes2plot;
        histogram(amplitude_baseline_rate*100)
        xlabel('%')
        title('Osc %of Baseline Rate')

        subplot(4,9,16)

        non_nan = ~isnan(phase_entrainment(:,1)) & ~isnan(abs(sin_fit_param(:,1)./sin_fit_param(:,5)));
        [~, TRMF1] = rmoutliers(100*abs(sin_fit_param(:,1)./sin_fit_param(:,5)));

        if sum(non_nan & ~TRMF1)>0
            hold on
            [~, TRMF1] = rmoutliers(100*abs(sin_fit_param(:,1)./sin_fit_param(:,5)));
            plot(100*abs(sin_fit_param(~TRMF1,1)./sin_fit_param(~TRMF1,5)), phase_entrainment(~TRMF1,1),'k.' )
            [c, p]=corr(phase_entrainment(non_nan & ~TRMF1,1),abs(sin_fit_param(non_nan& ~TRMF1,1)./sin_fit_param(non_nan& ~TRMF1,5)));
            sig_index = phase_entrainment(:,2)>0.95;

            plot(100*abs(sin_fit_param(sig_index & non_nan & ~TRMF1,1)./sin_fit_param(sig_index & non_nan & ~TRMF1,5)),...
                phase_entrainment(sig_index & non_nan & ~TRMF1,1),'r.', 'MarkerSize', 16)
            xlabel('Osc %')
            ylabel('MVL')
            title(num2str([c p]))
        end

        subplot(4,9,17)
        % 
        % [~, loc] = max(fitted_theta_alignment, [],2);
        % max_phase_value = mod(2*pi*hist_edges_centers(loc)/.125, 2*pi);
        % el2ex = abs(loc-ceil(.5*smooth_window))<3 | abs(size(fitted_theta_alignment,2)-ceil(.5*smooth_window) - loc)<3;
        amplitud = run_events_table.SinFit(:,1);
        phase_shift = run_events_table.SinFit(:,3);
        angle1 = phase_shift + (amplitud<0)*pi;
        angle1 = mod(angle1+pi/2 ,2*pi);
        angle2 = mod(-run_events_table.phase_entrainment(:,4), 2*pi);
        plot(angle1/pi, angle2/pi, '.k')

        hold on
        plot([0 2], [0 2], 'r')
        % plot(max_phase_value(sig_index & ~el2ex), phase_entrainment(sig_index & ~el2ex,3), '.r', 'MarkerSize', 16)
        % [c, p]= circ_corrcc(phase_entrainment(~el2ex,3), max_phase_value(~el2ex));
        [c, p]= circ_corrcc(angle1, angle2);

        title(num2str([c p]))
        xlabel('Phase max oscilation')
        ylabel('Mean Phase')



        subplot(4,9,18)
        yticks([])
        yyaxis right
        temporal_shifts_norm = diag(1./phase_entrainment(:,1))*temporal_shifts;
        [max_value,loc] = max(temporal_shifts_norm,[],2);
        matrix2plot = temporal_shifts_norm(max_value>1,:);
        matrix2plot = matrix2plot(~any(isnan(matrix2plot),2),:);
        mean_response = mean(matrix2plot);
        ci =1.96*std(matrix2plot)/sqrt(size(matrix2plot,1));
        fill([temporal_shift_values fliplr(temporal_shift_values)],[mean_response-ci fliplr(mean_response+ci)],'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
        hold on
        plot(temporal_shift_values, mean_response, 'r')
        xlabel('Temporal Phase Shift')
        ylabel('MVL')





        subplot(4,9,25:27)
        % variable2sort = amplitude_baseline_rate(:,1); %amplitu_baseline_rate
        variable2sort = phase_entrainment(:,1) ;
        hold on

        pctls2use = [0 100; 50 100; 80 100; 90 100];
        ranges2plot = pctls2use;
        for nn = 1:size(pctls2use,1)
            ranges2plot(nn,:)= prctile(variable2sort,pctls2use(nn,:));
        end
        legends= cell(size(ranges2plot,1)-1,1);

        for r_index = 1:size(ranges2plot,1)
            index = variable2sort>=ranges2plot(r_index,1) &  variable2sort<=ranges2plot(r_index,2);
            if sum(index)>1
                rate_osc = mean(this_neuron_autoccorrelograms(index,:));
            elseif sum(index)==1
                rate_osc = this_neuron_autoccorrelograms(index,:);
            else
                rate_osc = nan(size(autocorr_centers));
            end
            rate_osc(autocorr_centers>.02) = smooth(rate_osc(autocorr_centers>.02), .04/autocorr_bin_size);
            rate_osc(autocorr_centers<-.02) = smooth(rate_osc(autocorr_centers<-.02), .04/autocorr_bin_size);
            rate_osc = zscore(rate_osc);
            plot(autocorr_centers,rate_osc, 'Color', [1 0 0] + [-1 0 1]*r_index/size(ranges2plot,1), 'LineWidth',1)
            legends{r_index} = ['Ptcl ', num2str(pctls2use(r_index,1)) , ' to ', num2str(pctls2use(r_index,2))];
        end
        legend(legends)
        xticks([autocorr_range(1) -1/8 1/8 autocorr_range(2)])
        xlabel('Time (s)')


        subplot(4,9,34:36)

        median_fit  = median(sin_fit_param(:,6));
        [~, median_loc] = min(abs(median_fit-sin_fit_param(:,6)));
        median_loc = median_loc(1);

        % plot(hist_edges_centers, real_raw_theta_alignment(median_loc,:)*bin_size/smoothing_time, '.k')
        plot(hist_edges_centers, mean_raw_theta_alignment(median_loc,:)*bin_size/smoothing_time, 'k')
        hold on
        plot(hist_edges_centers, mean_theta_alignment(median_loc,:), 'b')
        plot(hist_edges_centers, fitted_theta_alignment(median_loc,:), 'r')
        legend({'Raw', 'Mov Average', 'Fit'})
        ylabel('Rate (Hz)')
        xlabel('Time (s)')
        legend('boxoff')   

        % 
        % subplot(4,9,36)
        % polarhistogram(phase_entrainment(:,3), -pi:(pi/12):pi, 'FaceColor', 'r', 'EdgeColor','none', 'Normalization','percentage')
        % hold on
        % polarhistogram(phase_shift+pi/2, -pi:(pi/12):pi, 'FaceColor', 'b', 'EdgeColor','none', 'Normalization','percentage')
        % legend({'MV angle', 'phase_shift'}, 'Location','best')
        legend('boxoff')




    else


        title(['ID =', num2str(id), ' CH = ' , num2str(ch), ' ',  ct])
    end
    saveas(gcf,['ID =', num2str(id), ' CH = ' , num2str(ch), ' ',  ct, '.jpg'] )
    saveas(gcf,['ID =', num2str(id), ' CH = ' , num2str(ch), ' ',  ct, '.fig'] )

    close gcf
    run_events_table_per_cell{dn} =run_events_table;


end
save('run_events_table_per_cell', 'run_events_table_per_cell','config')

%% puting tables together
ALL_TABLES = [];
trans_table = run_events_table_per_cell{1};
trans_table_variables = trans_table.Properties.VariableNames;
trans_table_variables(ismember(trans_table_variables,'RunMaxate')) = {'RunMaxRate'};
trans_table.Properties.VariableNames =trans_table_variables;
ALL_TABLES = [ALL_TABLES;trans_table];
for j=2:size(run_events_table_per_cell,1)
    trans_table = run_events_table_per_cell{j};
    if size(trans_table,1)>0
        trans_table_variables = trans_table.Properties.VariableNames;
        trans_table_variables(ismember(trans_table_variables,'RunMaxate')) = {'RunMaxRate'};
        trans_table.Properties.VariableNames =trans_table_variables;

        ALL_TABLES = [ALL_TABLES;trans_table];
    end
end

phase_shift = mod(ALL_TABLES.SinFit(:,3), 2*pi);
amplitud    = ALL_TABLES.SinFit(:,1);
phase_shift = phase_shift + (amplitud<0)*pi;
phase_shift = mod(phase_shift, 2*pi);

ALL_TABLES.SinFit(:,3) =phase_shift;
ALL_TABLES.SinFit(:,1) = abs(ALL_TABLES.SinFit(:,1));
%%



%%
% a(1) sin(b(2)x +c(3)) + d(4)x + e(5)

figure
% plot(abs(ALL_TABLES.SinFit(:,1)./ALL_TABLES.SinFit(:,5)), phase_shift, 'k.')
hold on
swarmchart(round(abs(ALL_TABLES.SinFit(:,1)./ALL_TABLES.SinFit(:,5)),1), phase_shift+2*pi, 'b.')


data2pca = [phase_shift,abs(ALL_TABLES.SinFit(:,1)./ALL_TABLES.SinFit(:,5)), abs(ALL_TABLES.SinFit(:,1)), ALL_TABLES.SinFit(:,[2 4 5]) ALL_TABLES.lfp_param];

variable_names =  {'phase_shift','amp_ratio','osc_rate','freq_neuron','slope_linear_trend', 'base_rate', 'lfp_range', 'mean_evnelope', 'mean_freq'}   ;
data2pca_table = array2table(data2pca);
data2pca_table.Properties.VariableNames = variable_names;

for j=1:size(data2pca,2)
    data2pca(:,j) = (data2pca(:,j)-nanmean(data2pca(:,j)))/nanstd(data2pca(:,j));
end


[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(data2pca);


cases2exclude = ALL_TABLES.SinFit(:,5)<0 | isnan( data2pca_table.mean_freq) | data2pca_table.freq_neuron<6 ;
osc_trhes = .7;
amp_rate = ALL_TABLES.phase_entrainment(:,1);
% amp_rate = data2pca_table.osc_rate./data2pca_table.base_rate;

z_socre_range = -2:.05:10;
z_socre_range_centers = (z_socre_range(1:end-1) + z_socre_range(2:end))/2;
figure
subplot(2,2,1)
[c_no,p_no] = corr(data2pca_table.base_rate(amp_rate<=osc_trhes & ~cases2exclude), data2pca_table.lfp_range(amp_rate<=osc_trhes &~cases2exclude), ...
    'type','Spearman');
loglog(data2pca_table.base_rate(amp_rate<=osc_trhes & ~cases2exclude), data2pca_table.lfp_range(amp_rate<=osc_trhes &~cases2exclude), 'k.')
hold on
[c_o,p_o] = corr(data2pca_table.base_rate(amp_rate>osc_trhes & ~cases2exclude), data2pca_table.lfp_range(amp_rate>osc_trhes &~cases2exclude), ...
    'type','Spearman');
loglog(data2pca_table.base_rate(amp_rate>osc_trhes & ~cases2exclude), data2pca_table.lfp_range(amp_rate>osc_trhes &~cases2exclude), 'r.')
xlabel('Neuron Base Rate (Hz)')
ylabel('LFP range')
title(num2str([c_no,p_no,c_o,p_o]))

subplot(2,2,2)
[c_no,p_no] = corr(data2pca_table.osc_rate(amp_rate<=osc_trhes & ~cases2exclude), data2pca_table.lfp_range(amp_rate<=osc_trhes & ~cases2exclude), ...
    'type','Spearman');
loglog(data2pca_table.osc_rate(amp_rate<=osc_trhes & ~cases2exclude), data2pca_table.lfp_range(amp_rate<=osc_trhes & ~cases2exclude), 'k.')
hold on
[c_o,p_o] = corr(data2pca_table.osc_rate(amp_rate>osc_trhes & ~cases2exclude), data2pca_table.lfp_range(amp_rate>osc_trhes & ~cases2exclude), ...
    'type','Spearman');
loglog(data2pca_table.osc_rate(amp_rate>osc_trhes & ~cases2exclude), data2pca_table.lfp_range(amp_rate>osc_trhes & ~cases2exclude), 'r.')
title(num2str([c_no,p_no,c_o,p_o]))
xlabel('Neuron Osc Rate (Hz)')
ylabel('LFP range')


subplot(2,2,3)
[c_no,p_no] = corr(data2pca_table.freq_neuron(amp_rate<=osc_trhes & ~cases2exclude), data2pca_table.mean_freq(amp_rate<=osc_trhes & ~cases2exclude), ...
    'type','Spearman');
loglog(data2pca_table.freq_neuron(amp_rate<=osc_trhes & ~cases2exclude), data2pca_table.mean_freq(amp_rate<=osc_trhes & ~cases2exclude), 'k.')
hold on
[c_o,p_o] = corr(data2pca_table.freq_neuron(amp_rate>osc_trhes & ~cases2exclude), data2pca_table.mean_freq(amp_rate>osc_trhes & ~cases2exclude), ...
    'type','Spearman');
loglog(data2pca_table.freq_neuron(amp_rate>osc_trhes & ~cases2exclude), data2pca_table.mean_freq(amp_rate>osc_trhes & ~cases2exclude), 'r.')
title(num2str([c_no,p_no,c_o,p_o]))
xlabel('Neuron Osc Rate (Hz)')
ylabel('LFP Freq (Hz)')


subplot(2,2,4)
[c_no,p_no] = corr(data2pca_table.base_rate(amp_rate<=osc_trhes & ~cases2exclude), data2pca_table.osc_rate(amp_rate<=osc_trhes & ~cases2exclude), ...
    'type','Spearman');
loglog(data2pca_table.base_rate(amp_rate<=osc_trhes & ~cases2exclude), data2pca_table.osc_rate(amp_rate<=osc_trhes & ~cases2exclude), 'k.')
hold on
[c_o,p_o] = corr(data2pca_table.base_rate(amp_rate>osc_trhes & ~cases2exclude), data2pca_table.osc_rate(amp_rate>osc_trhes & ~cases2exclude), ...
    'type','Spearman');
loglog(data2pca_table.base_rate(amp_rate>osc_trhes & ~cases2exclude), data2pca_table.osc_rate(amp_rate>osc_trhes & ~cases2exclude), 'r.')
title(num2str([c_no,p_no,c_o,p_o]))
xlabel('Neuron Base Rate (Hz)')
ylabel('Neuron Osc Rate (Hz)')



%%
angle_bin_size = pi/8;
angle_edges         =  0:angle_bin_size:(2*pi);
angle_edges_centers = (angle_edges(1:end-1) + angle_edges(2:end))/2;
cell_ids            = unique(ALL_TABLES.Id);
angle_histograms    = nan(numel(cell_ids), numel(angle_edges)-1,2);
oscillatory_events = false(size(ALL_TABLES,1),1);
cell_firing_rate = nan(numel(cell_ids),2);

id_index = 1;
for id = cell_ids'

    this_cell_indexes       = find(ALL_TABLES.Id==id & ~isnan(ALL_TABLES.SinFit(:,3)));
    this_cell_table         = ALL_TABLES(this_cell_indexes,:);
    % this_cell_amplitud_rate = abs(this_cell_table.SinFit(:,1)./this_cell_table.SinFit(:,5));
    this_cell_amplitud_rate =  this_cell_table.phase_entrainment(:,2);
    this_cell_angles = ALL_TABLES.phase_entrainment(ALL_TABLES.Id==id & ~isnan(ALL_TABLES.phase_entrainment(:,3)),3);
    % this_cell_angles =  this_cell_table.phase_entrainment(:,3);
    this_cell_rate = ALL_TABLES.RunNumSpikes(this_cell_indexes)./ALL_TABLES.RunLength(this_cell_indexes);
        this_cell_angles = mod(this_cell_angles, 2*pi);
    % figure
    % plot(this_cell_table.SinFit(osc_index,1), this_cell_table.lfp_param(osc_index,1), 'r.')
    % hold on
    %  plot(this_cell_table.SinFit(~osc_index,1), this_cell_table.lfp_param(~osc_index,1), 'k.')
    %
    % plot(this_cell_table.SinFit(:,5), this_cell_table.lfp_param(:,1), '.')
    % pct80       = prctile(this_cell_amplitud_rate,80);
    % osc_index   = this_cell_amplitud_rate>.95;
    osc_index =  this_cell_table.phase_entrainment(:,1)>.7;

    oscillatory_events(this_cell_indexes(osc_index)) = true;
    if sum(osc_index)>1
        % angle_histograms(id_index,:,1)  = histcounts(this_cell_table.SinFit(osc_index,3), angle_edges)/sum(osc_index);
        % angle_histograms(id_index,:,2)  = histcounts(this_cell_table.SinFit(~osc_index,3), angle_edges)/sum(~osc_index);
        angle_histograms(id_index,:,1)  = histcounts(this_cell_angles(osc_index), angle_edges)/sum(osc_index);
        angle_histograms(id_index,:,2)  = histcounts(this_cell_angles(~osc_index), angle_edges)/sum(~osc_index);


        % cell_firing_rate(id_index,1)    = nanmean(ALL_TABLES.SinFit(ALL_TABLES.Id==id,5));
    end
    cell_firing_rate(id_index,2)     = cluster_data.fr(cluster_data.cluster_id==id );
    cell_firing_rate(id_index,1)     = mean(this_cell_rate);

    id_index = id_index+1;
end

os_mean = circ_mean(ALL_TABLES.SinFit(oscillatory_events,3));
nonos_mean = circ_mean(ALL_TABLES.SinFit(~oscillatory_events & ~isnan(ALL_TABLES.SinFit(:,3)),3));

figure
subplot(1,2,1)
matrix2plot = squeeze(angle_histograms(:,:,1));
matrix2plot = matrix2plot(~any(isnan(matrix2plot),2),:);
mean_response = mean(matrix2plot);
ci =1.96*std(matrix2plot)/sqrt(size(matrix2plot,1));
fill(180*[angle_edges_centers fliplr(angle_edges_centers)]/pi,[mean_response-ci fliplr(mean_response+ci)],'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
hold on
plot(180*angle_edges_centers/pi, mean_response, 'r')
fill(360-angle_bin_size*180/pi+ 180*[angle_edges_centers fliplr(angle_edges_centers)]/pi,[mean_response-ci fliplr(mean_response+ci)],'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
hold on
plot(360-angle_bin_size*180/pi+ 180*angle_edges_centers/pi, mean_response, 'r')

matrix2plot = squeeze(angle_histograms(:,:,2));
matrix2plot = matrix2plot(~any(isnan(matrix2plot),2),:);
mean_response = mean(matrix2plot);
ci =1.96*std(matrix2plot)/sqrt(size(matrix2plot,1));
fill(180*[angle_edges_centers fliplr(angle_edges_centers)]/pi,[mean_response-ci fliplr(mean_response+ci)],'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
hold on
plot(180*angle_edges_centers/pi, mean_response, 'k')
fill(360-angle_bin_size*180/pi+ 180*[angle_edges_centers fliplr(angle_edges_centers)]/pi,[mean_response-ci fliplr(mean_response+ci)],'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
hold on
plot(360-angle_bin_size*180/pi+ 180*angle_edges_centers/pi, mean_response, 'k')
axis tight
y_lim = ylim;

plot([os_mean os_mean]*180/pi, y_lim, 'r')
plot([os_mean os_mean]*180/pi  + 360, y_lim, 'r')

plot([nonos_mean nonos_mean]*180/pi, y_lim, 'k')
plot([nonos_mean nonos_mean]*180/pi  + 360, y_lim, 'k')


plot((0:.01:4*pi)*180/pi, (sin((0:.01:4*pi) -pi/2)*range(y_lim)/2) + mean(y_lim), ':k', 'LineWidth',2)

subplot(1,2,2)
polarhistogram(phase_shift(oscillatory_events), -pi:(pi/18):pi, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', .5, 'Normalization','percentage')
hold on
polarhistogram(phase_shift(~oscillatory_events), -pi:(pi/18):pi, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', .5, 'Normalization','percentage')
axis_pointer= gca;
r_lim  = axis_pointer.RLim;

polarplot([os_mean,os_mean],r_lim, 'r', 'LineWidth',2 )
polarplot([nonos_mean,nonos_mean],r_lim, 'k', 'LineWidth',2 )
central_os     = circ_rtest(phase_shift(oscillatory_events));
central_nonos  = circ_rtest(phase_shift(~oscillatory_events & ~isnan(phase_shift)));
[pval, table] = circ_wwtest(phase_shift(oscillatory_events), phase_shift(~oscillatory_events & ~isnan(phase_shift)));
title(num2str([central_os central_nonos pval table{2,5}]))

%% estiamte phase entreinment per neuron
figure

phase_bins = 0:(pi/2):2*pi;
ALL_TABLES.normalize_rates = abs(ALL_TABLES.SinFit(:, [1 5]));
cell_id = unique(ALL_TABLES.Id)';
non_nan = ~isnan(ALL_TABLES.normalize_rates);
line_types = {'','--', ':'};
for this_id = cell_id
    index = ALL_TABLES.Id == this_id;
    ALL_TABLES.normalize_rates(index & non_nan(:,1),1) = zscore(ALL_TABLES.normalize_rates(index & non_nan(:,1),1)) ;
    ALL_TABLES.normalize_rates(index & non_nan(:,2),2) = zscore(ALL_TABLES.normalize_rates(index & non_nan(:,2),2)) ;
end

for j=1:numel(phase_bins)-1

    these_phases_index =  ALL_TABLES.SinFit(:,3)>=phase_bins(j) & ALL_TABLES.SinFit(:,3)<phase_bins(j+1);
    subplot(1,3,1)
    hold on
    sub_r = ALL_TABLES.phase_entrainment(these_phases_index,1);
    [f,x] = ecdf(sub_r);
    plot(x,f, line_types{mod(j,3)+1}, 'Color',[1 0 0 ] + [-1 0 0]*j/(numel(phase_bins)-1))


    subplot(1,3,2)
    hold on
    sub_r =  ALL_TABLES.normalize_rates(these_phases_index,1);
    [f,x] = ecdf(sub_r);
    plot(x,f, line_types{mod(j,3)+1}, 'Color',[1 0 0 ] + [-1 0 0]*j/(numel(phase_bins)-1))


    subplot(1,3,3)
    hold on
    sub_r =  ALL_TABLES.normalize_rates(these_phases_index,2);
    [f,x] = ecdf(sub_r);
    plot(x,f, line_types{mod(j,3)+1}, 'Color',[1 0 0 ] + [-1 0 0]*j/(numel(phase_bins)-1))

    % plot(0.005:0.01:1, r_counts/sum(r_counts))

end

%%
mean_temporal_shifts = nan(numel(id_list), 2,size(ALL_TABLES.temporal_shifts,2));
mean_theta_lok = nan(numel(id_list), 2,2*size(ALL_TABLES.ThetaPsth,2));
CLUSTER_TYPES = string(ALL_TABLES.ct);
osc_param = nan(numel(id_list),22);
id_list = unique(ALL_TABLES.Id)';
%%
for nn = 1:numel(id_list)
    id      = id_list(nn);
    index   = ALL_TABLES.Id == id;
    ch      = unique(ALL_TABLES.Ch(index));
    c_type  = unique(CLUSTER_TYPES(index,:));

    if sum(index)>1
        
        phases = ALL_TABLES.SinFit(index,3);
        psth = ALL_TABLES.ThetaPsth(index,:);
        SinFit = ALL_TABLES.SinFit(index,:);
        [phases_values, sorted_phases] = sort(phases);
        figure('units','normalized','outerposition',[0 0 1 1])
        colormap(1-gray)
        subplot(1,5,1)
        imagesc(hist_edges_centers,1:size(psth,1), psth )
        axis xy

        subplot(1,5,2)
        phaseintime1 = -.125*phases_values/(2*pi);
        phaseintime2 = -.125*(phases_values-2*pi)/(2*pi);


        opositephaseintime1 = -.125*mod(pi+phases_values, 2*pi)/(2*pi);
        opositephaseintime2 = -.125*(-pi+phases_values)/(2*pi);

        opositephaseintime = max(opositephaseintime1,opositephaseintime2);
        imagesc(hist_edges_centers,1:size(psth,1), psth(sorted_phases,:) )
        axis xy
        hold on
        plot(phaseintime1, 1:size(psth,1), '.r')
        plot(phaseintime2, 1:size(psth,1), '.r')
        plot(opositephaseintime, 1:size(psth,1), '.b')

        xlim([-.125 .125])
        peak_values = nan(numel(sorted_phases),2);

        for phase_event = 1:numel(sorted_phases)
            this_psth = psth(sorted_phases(phase_event),:);
            peak_index = hist_edges_centers>=phaseintime1(phase_event) & ...
                hist_edges_centers<=opositephaseintime(phase_event);
            trough_index =  hist_edges_centers>=opositephaseintime(phase_event) & ...
                hist_edges_centers<=phaseintime2(phase_event);
            peak_rate   = mean(this_psth(peak_index));
            trough_rate = mean(this_psth(trough_index));
            peak_values(phase_event,:) = [peak_rate trough_rate];
        end





        subplot(1,5,3)
        % osc_amplitud        = abs(SinFit(sorted_phases,1));
        % sorted_osc_amplitud = abs(osc_amplitud(sorted_phases));
        % all_events          = 1:numel(sorted_osc_amplitud);
        % all_events_woo      = all_events;
        % [~, OtR]            = rmoutliers(sorted_osc_amplitud);
        % all_events_woo(OtR) = [];
        % sorted_osc_amplitud(OtR)   = [];
        % sorted_osc_amplitud = interp1(all_events_woo,sorted_osc_amplitud,all_events);
        % plot(abs(sorted_osc_amplitud),1:size(psth,1))
        plot(peak_values(:,1),1:size(psth,1))




        shifted_psth = nan(size(psth,1), 2*size(psth,2));
        all_closes_points = [];
        for j=1:numel(phases)
            [~,closes_point] = min(abs(hist_edges_centers -opositephaseintime(j) ));
            all_closes_points = [all_closes_points;closes_point]   ;
            shifted_psth(sorted_phases(j),(numel(hist_edges_centers)/2-(closes_point-numel(hist_edges_centers)/2) + 1):...
                (numel(hist_edges_centers)/2-(closes_point-numel(hist_edges_centers)/2)  + size(psth,2))) = psth(sorted_phases(j),:);
        end
        mvl             = ALL_TABLES.phase_entrainment(index,1);
        [~,sorted_by_r] = sort(mvl);
        subplot(5,5,3:5:14)
        double_time = [hist_edges_centers,hist_edges_centers + range(hist_edges_centers)+bin_size]-hist_edges_centers(end);
        imagesc(double_time, 1:size(psth,1),shifted_psth(sorted_by_r,:));
        axis xy
        xlim([-.2 .2])


        subplot(5,5,[18 23])
        pct50r = prctile(mvl,[40 50]);
        all_values_plot = [];
        osc_param(nn,1)      = max(mvl);
        osc_param(nn,2:3)    = pct50r;

        matrix2plot     = shifted_psth(mvl<pct50r(1), :);
        mean_response   = nanmean(matrix2plot);
        [~,~,ci] = ttest(matrix2plot);
        if size(ci,2)>2

            all_values_plot =[all_values_plot;[mean_response;ci]];
            fill([double_time(~any(isnan(ci))) fliplr(double_time(~any(isnan(ci))))],...
                [ci(1,~any(isnan(ci))) fliplr(ci(2,~any(isnan(ci))))],...
                'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
            hold on
            plot(double_time, mean_response, 'k')
            if numel(mean_response)>1
                mean_theta_lok(nn,1,:) = mean_response;
            end
            matrix2plot     = shifted_psth(mvl>pct50r(2), :);
            mean_response   = nanmean(matrix2plot);
            [~,~,ci] = ttest(matrix2plot);

            fill([double_time(~any(isnan(ci))) fliplr(double_time(~any(isnan(ci))))],...
                [ci(1,~any(isnan(ci))) fliplr(ci(2,~any(isnan(ci))))],...
                'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
            hold on
            plot(double_time, mean_response, 'r')
            all_values_plot =[all_values_plot;[mean_response;ci]];
            axis tight
            all_values_plot = all_values_plot(~isnan(all_values_plot));
            range2plot = prctile(all_values_plot, [1 99]);
            xlim([-.2 .2])
            ylim(range2plot)
        end

        subplot(5,5,4:5:14)
        temporal_shifts_norm = diag(1./ALL_TABLES.phase_entrainment(index,1))*ALL_TABLES.temporal_shifts(index,:);
        imagesc(temporal_shift_values, 1:size(temporal_shifts_norm,1), temporal_shifts_norm(sorted_by_r,:)   )
        axis xy

        subplot(5,5,[19 24])
        index_r = mvl<pct50r(1);
        matrix2plot     = temporal_shifts_norm(index_r,:);
        mean_response   = nanmean(matrix2plot);
        [~,~,ci] = ttest(matrix2plot);
        if size(ci,2)>2
            fill([temporal_shift_values(~any(isnan(ci))) fliplr(temporal_shift_values(~any(isnan(ci))))],...
                [ci(1,~any(isnan(ci))) fliplr(ci(2,~any(isnan(ci))))],...
                'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
            hold on
            plot(temporal_shift_values, mean_response, 'k')
            axis tight
            if numel(mean_response)>1
                mean_temporal_shifts(nn,1,:) = mean_response;
            end
            index_r         = mvl>=pct50r(2);
            matrix2plot     = temporal_shifts_norm(index_r,:);
            mean_response   = nanmean(matrix2plot);
            [~,~,ci] = ttest(matrix2plot);
            fill([temporal_shift_values(~any(isnan(ci))) fliplr(temporal_shift_values(~any(isnan(ci))))],...
                [ci(1,~any(isnan(ci))) fliplr(ci(2,~any(isnan(ci))))],...
                'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
            hold on
            plot(temporal_shift_values, mean_response, 'r')

            plot(temporal_shift_values, temporal_shift_values*0 +1 , ':k')
            xlim([-.2 .2])
            if numel(mean_response)>1
                mean_temporal_shifts(nn,2,:) = mean_response;
            end
        end

        subplot(5,5,5)
        histogram(mvl, 0:0.025:1)
        % xlabel('MVL')

        subplot(5,5,10)
        phase_preference  = ALL_TABLES.SinFit(index,3);
        os_angles       = phase_preference(mvl>=pct50r(2));
        nonos_angles    = phase_preference(mvl<pct50r(1));
        polarhistogram(os_angles, -pi:(pi/16):pi,'FaceColor', 'r', 'EdgeColor','none', 'Normalization','percentage')
        hold on
        polarhistogram(nonos_angles, -pi:(pi/16):pi,'FaceColor', 'b', 'EdgeColor','none', 'Normalization','percentage')

        osc_param(nn,3) = circ_mean(os_angles);
        osc_param(nn,4) = circ_mean(nonos_angles);

        polar_axis = gca;
        polarplot([osc_param(nn,3) osc_param(nn,3)], 1.1*polar_axis.RLim, 'r')
        polarplot([osc_param(nn,4) osc_param(nn,4)], 1.1*polar_axis.RLim, 'b')
        axis tight
        a = NaN;
        b = NaN;
        os_p = NaN;
        os_r = NaN;
        nonos_p = NaN;
        nonos_r = NaN;
        if numel(os_angles)>2 &&  numel(nonos_angles)>2
            os_p = circ_rtest(os_angles);
            os_r = circ_r(os_angles);
            
            nonos_p=circ_rtest(nonos_angles);
            nonos_r = circ_r(nonos_angles);
            
            if os_p<0.05 && nonos_p <0.05   
                [a,b] =circ_wwtest(os_angles,nonos_angles );
                b = b{2,5};
            end
            this_tilte = num2str(round([os_p os_r nonos_p  nonos_r a b],2));            
            title(strrep(this_tilte, '          ', '/')  )
        end
        osc_param(nn,5:10) = [os_p os_r nonos_p nonos_r a b];

        subplot(5,5,15)
        %
        amp_ratio  =abs(ALL_TABLES.SinFit(index,1))./abs(ALL_TABLES.SinFit(index,5));
        [~, ratio_ol] = rmoutliers(amp_ratio);

        %
        % amp_ratio  =abs(ALL_TABLES.SinFit(index,5));
        % [~, ratio_ol] = rmoutliers(amp_ratio);

        % plot(amp_ratio(~ratio_ol), mvl(~ratio_ol), 'k.')
        % [c,p]=corr(amp_ratio(~ratio_ol), mvl(~ratio_ol))
        % osc_param(nn,5:6) =[c,p];

        phase_preference  = ALL_TABLES.phase_entrainment(index,3);
         os_angles       = phase_preference(mvl>=pct50r(2));
        nonos_angles    = phase_preference(mvl<pct50r(1));
        polarhistogram(os_angles, -pi:(pi/16):pi,'FaceColor', 'r', 'EdgeColor','none', 'Normalization','percentage')
        hold on
        polarhistogram(nonos_angles, -pi:(pi/16):pi,'FaceColor', 'b', 'EdgeColor','none', 'Normalization','percentage')

        osc_param(nn,11) = circ_mean(phase_preference(mvl>=pct50r(2)));
        osc_param(nn,12) = circ_mean(phase_preference(mvl<pct50r(1)));

        polar_axis = gca;
        polarplot([osc_param(nn,11) osc_param(nn,11)], 1.1*polar_axis.RLim, 'r')
        polarplot([osc_param(nn,12) osc_param(nn,12)], 1.1*polar_axis.RLim, 'b')

        a = NaN;
        b = NaN;
        os_p = NaN;
        os_r = NaN;
        nonos_p = NaN;
        nonos_r = NaN;
        if numel(os_angles)>2 &&  numel(nonos_angles)>2
            os_p = circ_rtest(os_angles);
            os_r = circ_r(os_angles);
            
            nonos_p=circ_rtest(nonos_angles);
            nonos_r = circ_r(nonos_angles);
            
            if os_p<0.05 && nonos_p <0.05   
                [a,b] =circ_wwtest(os_angles,nonos_angles );
                b = b{2,5};
            end
             this_tilte = num2str(round([os_p os_r nonos_p  nonos_r a b],2));            
            title(strrep(this_tilte, '          ', '/')  )
        end
        osc_param(nn,13:18) = [os_p os_r nonos_p  nonos_r a b];

        axis tight


        subplot(5,5,20)

        lfp_amp= ALL_TABLES.lfp_param(index,1);
        os_amp = abs(ALL_TABLES.SinFit(index,1));
        [~, ol2r_os] = rmoutliers(os_amp);
        [~, ol2r_lfp] = rmoutliers(lfp_amp);
        plot(lfp_amp(~ol2r_os & ~ol2r_lfp), os_amp(~ol2r_os & ~ol2r_lfp), 'k.')
        hold on
        [c,p]=corr(lfp_amp(~ol2r_os & ~ol2r_lfp), os_amp(~ol2r_os & ~ol2r_lfp));
        title(num2str([c p]))
        osc_param(nn,19:20) =[c,p];
        xlabel('LFP AMP')
        ylabel('Os Am')

        subplot(5,5,25)

        lfp_amp= ALL_TABLES.lfp_param(index,1);
        os_amp = abs(ALL_TABLES.SinFit(index,5));
        [~, ol2r_os] = rmoutliers(os_amp);
        [~, ol2r_lfp] = rmoutliers(lfp_amp);
        plot(lfp_amp(~ol2r_os & ~ol2r_lfp), os_amp(~ol2r_os & ~ol2r_lfp), 'k.')
        hold on
        if sum(~ol2r_os & ~ol2r_lfp & mvl<pct50r(1))>0
            [c,p]=corr(lfp_amp(~ol2r_os & ~ol2r_lfp & mvl<pct50r(1)), os_amp(~ol2r_os & ~ol2r_lfp & mvl<pct50r(1)));
            osc_param(nn,21:22) =[c,p];

            title(num2str([c p]))
        end
        xlabel('LFP AMO')
        ylabel('Base Rate')
        pause(.1)
    end
    saveas(gcf,['ID =', num2str(id), ' CH = ' , num2str(ch), ' ',  char(c_type), '.jpg'] )
    saveas(gcf,['ID =', num2str(id), ' CH = ' , num2str(ch), ' ',  char(c_type), '.fig'] )
    close(gcf)

end
% freq = abs(ALL_TABLES.SinFit(index,2)/(2*pi));
% lfp_freq = ALL_TABLES.lfp_param(index,3);
%
% [~, freq2r] =rmoutliers(freq);
% [~, flp2r] =rmoutliers(lfp_freq);
%
% plot(freq(~freq2r & ~flp2r & mvl>.5 ),lfp_freq(~freq2r & ~flp2r & mvl>.5 ), '.r' )
% [c,p]=corr(freq(~freq2r & ~flp2r & mvl<.5 ),lfp_freq(~freq2r & ~flp2r & mvl<.5 ), 'type','Spearman')


% polarhistogram(ALL_TABLES.SinFit(index,3), -pi:(pi/16):pi, 'EdgeColor','none')

% bas_amp             = abs(SinFit(sorted_phases,5));
% sorted_bas_amp      = abs(bas_amp(sorted_phases));
% all_events          = 1:numel(bas_amp);
% all_events_woo      = all_events;
% [~, OtR]            = rmoutliers(sorted_bas_amp);
% all_events_woo(OtR) = [];
% sorted_bas_amp(OtR)   = [];
% sorted_bas_amp = interp1(all_events_woo,sorted_bas_amp,all_events);
% plot(abs(sorted_bas_amp),1:size(psth,1))
% plot(peak_values(:,2),1:size(psth,1))
%
% subplot(1,5,5)
% % plot(abs(sorted_osc_amplitud./sorted_bas_amp),1:size(psth,1))
% plot((peak_values(:,2)-peak_values(:,1))./(peak_values(:,2)+peak_values(:,1)),1:size(psth,1))
%
% figure
% plot(phases, (peak_values(:,2)-peak_values(:,1))./(peak_values(:,2)+peak_values(:,1)),'.')
% RI = (peak_values(:,2)-peak_values(:,1))./(peak_values(:,2)+peak_values(:,1)) ;
%
% [c,p]=circ_corrcl(phases(~isnan(RI)),RI(~isnan(RI)));
%% Session summary
figure

cells_actually_oss = osc_param(:,1)>.5;
subplot(2,2,3)
histogram(osc_param(cells_actually_oss,6), 0:0.05:1, 'FaceColor', 'r')
hold on
histogram(osc_param(cells_actually_oss,8), 0:0.05:1, 'FaceColor', 'k')



subplot(2,2,1)
histogram(osc_param(:,1), 0:0.05:1, 'FaceColor', 'k')


subplot(2,2,2)
polarhistogram(osc_param(cells_actually_oss,4), -pi:(pi/16):pi, 'FaceColor', 'k')
hold on
polarhistogram(osc_param(cells_actually_oss,3), -pi:(pi/16):pi, 'FaceColor', 'r')
ang_centers = (-(pi-pi/32):(pi/16):pi) ;
os_angles = histcounts(osc_param(cells_actually_oss,3), -pi:(pi/16):pi);
nonos_angles = histcounts(osc_param(cells_actually_oss,4), -pi:(pi/16):pi);
subplot(2,2,4)
plot(-pi:(pi/16):pi, sin((-pi:(pi/16):pi)  - pi/2), 'k')
yyaxis right
plot(ang_centers,os_angles)
hold on
plot(ang_centers,nonos_angles)
%%
shifted_psth = nan(size(psth,1), 2*size(psth,2));
all_closes_points = [];
for j=1:numel(phases)
    [~,closes_point] = min(abs(hist_edges_centers -opositephaseintime(j) ));
    all_closes_points = [all_closes_points;closes_point]   ;
    shifted_psth(sorted_phases(j),(numel(hist_edges_centers)/2-(closes_point-numel(hist_edges_centers)/2) + 1):...
        (numel(hist_edges_centers)/2-(closes_point-numel(hist_edges_centers)/2)  + size(psth,2))) = psth(sorted_phases(j),:);
end

[~,sorted_by_r] = sort(ALL_TABLES.phase_entrainment(index,1));
figure
subplot(1,2,1)
histogram(opositephaseintime, -.125:0.01:.125)
subplot(5,2,2:2:8)
double_time = [hist_edges_centers,hist_edges_centers + range(hist_edges_centers)+bin_size]-hist_edges_centers(end);
imagesc(double_time, 1:size(psth,1),shifted_psth(sorted_by_r,:));
axis xy
colormap(1-gray)

subplot(5,2,10)
matrix2plot        = shifted_psth;
mean_response   = nanmean(matrix2plot);
[~,~,ci] = ttest(matrix2plot);


fill([double_time(~any(isnan(ci))) fliplr(double_time(~any(isnan(ci))))],...
    [ci(1,~any(isnan(ci))) fliplr(ci(2,~any(isnan(ci))))],...
    'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
hold on
plot(double_time, mean_response, 'k')


%% next section are for cross correlation between cells
cross_corr_range        = [-2 2];
cross_corr_bin          = 0.01;
cross_correlogram_edges = cross_corr_range(1):cross_corr_bin:cross_corr_range(2);
cross_time              = (cross_correlogram_edges(1:end-1) + cross_correlogram_edges(2:end))/2;


figure
X = cell_firing_rate;
GMModel = fitgmdist(X,2);
y = [zeros(size(X,1),1);ones(size(X,1),1)];
rate_class1 = mean(GMModel.mu(1,:));
rate_class2 =  mean(GMModel.mu(2,:));
[~, max_col] = max([rate_class1 rate_class2]);

class                   = posterior(GMModel, X);
class                   = class(:,max_col)<.5;
h = gscatter(X(:,1),X(:,2),class);
hold on
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
g = gca;
fcontour(gmPDF,[0 35 -2 35 ],'MeshDensity',250, 'HandleVisibility','off')
% plot(GMModel.mu(1,:), 'r')
legend({'high firing','low firing'})
xlabel('Rate during firing episodes')
ylabel('Rate in entire session')



high_rate_neurons = cell_ids(~class);
low_rate_neurons  = cell_ids(class);
%% Estimating crosscorrelations during firing episodes
psth_indexes = [];
ALL_CROSS_CORR = [];
other_spikes            =  double(spike_times(ismember(spike_clusters,low_rate_neurons)))/30000;
other_spikes_clusters   = spike_clusters(ismember(spike_clusters,low_rate_neurons));

for id_index = 1:numel(high_rate_neurons)
    id = high_rate_neurons(id_index);

    this_cell_ranges = ALL_TABLES{ALL_TABLES.Id == id, {'RunStartTime','RunEndTime'}};
    % this_cell_ranges = this_cell_ranges(this_cell_ranges(:,2)-this_cell_ranges(:,1)>1, :);
    spike_times_sec =   double(spike_times(spike_clusters==id))/30000;
    all_psths = nan(size(this_cell_ranges,1), numel(cross_correlogram_edges)-1);
    active_spikes = false(size(this_cell_ranges,1),numel(low_rate_neurons));
    for rn=1:size(this_cell_ranges,1)
        current_ranges = this_cell_ranges(rn,:);

        this_spike_times        = spike_times_sec(spike_times_sec>=current_ranges(1) & spike_times_sec<=current_ranges(2));
        [counts,active_units] = groupcounts(other_spikes_clusters(other_spikes>=current_ranges(1) & other_spikes<=current_ranges(2)));
        active_spikes(rn, ismember(low_rate_neurons,active_units)) = counts;
        this_other_spike_times  =  other_spikes(other_spikes>=current_ranges(1)+cross_corr_range(1) & other_spikes<=current_ranges(2)+cross_corr_range(2));

        lags = crosscorrelogram(this_spike_times,this_other_spike_times, cross_corr_range);
        all_psths(rn,:) = histcounts(lags,cross_correlogram_edges)/sum(this_spike_times);
    end
    ALL_CROSS_CORR = [ALL_CROSS_CORR;all_psths];
    psth_indexes = [psth_indexes;[ones(size(all_psths,1),1)*id  this_cell_ranges active_spikes]];
end
%%
figure



matrix2plot = ALL_CROSS_CORR;
for j = 1:size(matrix2plot,1)
    matrix2plot(j,:) = zscore(matrix2plot(j,:));
end

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(matrix2plot);
for j=1:5
    subplot(1,5,j)
    plot(cross_time,COEFF(:,j))
    title(EXPLAINED(j))
end

%% positive and negative scores


figure
matrix2plot_neg = matrix2plot(SCORE(:,1)<0,:);
matrix2plot_pos = matrix2plot(SCORE(:,1)>0,:);

subplot(6,2,1:2:7)
imagesc(cross_time, 1:size(matrix2plot_neg,1), matrix2plot_neg)
clim([-2 2])
axis xy
subplot(6,2,9)
matrix2plot_neg = matrix2plot_neg(~any(isnan(matrix2plot_neg),2),:);
mean_response = mean(matrix2plot_neg);
ci =1.96*std(matrix2plot_neg)/sqrt(size(matrix2plot_neg,1));
fill([cross_time fliplr(cross_time)],[mean_response-ci fliplr(mean_response+ci)],'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
hold on
plot(cross_time, mean_response, 'r')



subplot(6,2,2:2:8 )
imagesc(cross_time, 1:size(matrix2plot_pos,1), matrix2plot_pos)
axis xy
clim([-2 2])
subplot(6,2,10  )
matrix2plot_pos = matrix2plot_pos(~any(isnan(matrix2plot_pos),2),:);
mean_response = mean(matrix2plot_pos);
ci =1.96*std(matrix2plot_pos)/sqrt(size(matrix2plot_pos,1));
fill([cross_time fliplr(cross_time)],[mean_response-ci fliplr(mean_response+ci)],'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
hold on
plot(cross_time, mean_response, 'r')



subplot(6,2, 11)

histogram(mean(psth_indexes(SCORE(:,1)<0, [2 3] ),2), 0:10:1500, 'Normalization','percentage', 'EdgeColor','none')
subplot(6,2, 12)
histogram(mean(psth_indexes(SCORE(:,1)>0, [2 3] ),2),0:10:1500, 'Normalization','percentage', 'EdgeColor','none')



%% distributin of scores per neuron

id_list = unique(psth_indexes(:,1))';
score_hist_range = [-2.1 2.1];
score_hist_bin = 0.1;
score_hist_edges = score_hist_range(1):score_hist_bin:score_hist_range(2);
score_hist_edges_centers = (score_hist_edges(1:end-1)+score_hist_edges(2:end))/2;
mean_score = nan(numel(id_list),1);

neurons_score_distr = nan(numel(id_list), numel(score_hist_edges)-1);
SCORE_ZSCORED = SCORE;

for j=1:size(SCORE,2)
    no_nan = ~isnan(SCORE(:,j));
    SCORE_ZSCORED(no_nan,j) = zscore(SCORE(no_nan,j));
end
z_scored_SCORE = zscore(SCORE(:,1));
neuron_contr2score = nan(numel(id_list),numel(low_rate_neurons));

for id_index = 1:numel(id_list)
    this_neuron_scores = SCORE_ZSCORED(psth_indexes(:,1)==id_list(id_index),2);
    neurons_pairs_this_neurons = psth_indexes(psth_indexes(:,1)==id_list(id_index), 4:end);
    figure
    [~, score_order] = sort(this_neuron_scores);
    imagesc(neurons_pairs_this_neurons(score_order,:))
    axis xy
    if size(neurons_pairs_this_neurons,1)>1
        anova_tsable = anova(neurons_pairs_this_neurons,this_neuron_scores);
        anova_tsable =stats(anova_tsable);
        neuron_contr2score(id_index,:) =anova_tsable.pValue(1:end-2);
    end
    mean_score(id_index) = mean(this_neuron_scores);
    % this_neuron_cross_corr = matrix2plot()
    neurons_score_distr(id_index,:) = histcounts(this_neuron_scores,score_hist_edges)/sum(psth_indexes(:,1)==id_list(id_index));
end
%% ploting distribution of score per neuron
figure
subplot(5,1,1:4)
matrix2plot = neurons_score_distr;
rows2include = sum(matrix2plot==0,2)<size(matrix2plot,2)*.1;
matrix2plot = matrix2plot(rows2include,:);


[mean_score_sorted, order] = sort(mean_score(rows2include));
imagesc(score_hist_edges_centers, 1:numel(order), matrix2plot(order,:))
hold on
plot(mean_score_sorted, 1:numel(order), 'w')
axis xy
clim([0 0.075])

subplot(5,1,5)

% matrix2plot = neurons_score_distr;

mean_response = mean(matrix2plot);
ci =1.96*std(matrix2plot)/sqrt(size(matrix2plot,1));
fill([score_hist_edges_centers fliplr(score_hist_edges_centers)],[mean_response-ci fliplr(mean_response+ci)],'r', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
hold on
plot(score_hist_edges_centers, mean_response, 'r')
%% relation between scores and theta osccliation
% id_index = 1:numel(high_rate_neurons)
factors2corr = 1:4;
oscillation_corr = [];
all_com_values = nan(numel(high_rate_neurons),6*numel(factors2corr)+1);

for id_index = 1:numel(high_rate_neurons)
    id       = high_rate_neurons(id_index);

    this_cell_table = ALL_TABLES(ALL_TABLES.Id == id,:);


    this_cells_oscillation_estate = abs(this_cell_table.SinFit(:,1)./this_cell_table.SinFit(:,5));
    this_cells_oscillation_phase = this_cell_table.SinFit(:,3);

    this_cell_cross_corr = ALL_CROSS_CORR(psth_indexes(:,1) == id,:);
    this_cell_score = SCORE_ZSCORED(psth_indexes(:,1) == id,2);
    this_cell_first_scores = SCORE_ZSCORED(psth_indexes(:,1) == id,factors2corr);
    oscillation_corr = [oscillation_corr;[this_cell_table.SinFit this_cell_first_scores]];

    all_com_values_this_cell =[];
    for fn = 1:size(this_cell_first_scores,2)

        positive_scores = this_cell_table.SinFit(this_cell_first_scores(:,fn) >0,3);
        negative_scores = this_cell_table.SinFit(this_cell_first_scores(:,fn)<0,3);

        if numel(positive_scores)>1 &&  numel(negative_scores)>1
            positive_mean = circ_mean(positive_scores);
            negative_mean = circ_mean(negative_scores);
            central_os     = circ_rtest(positive_scores);
            central_nonos  = circ_rtest(negative_scores);
            [pval, table] = circ_wwtest(positive_scores, negative_scores);
            all_com_values(id_index, (1:4) +(fn-1)*6 ) = [central_os central_nonos pval table{2,5}];
            all_com_values(id_index, (5:6) +(fn-1)*6 ) = [positive_mean negative_mean];
        else
            all_com_values(id_index, (1:6) +(fn-1)*6 ) = NaN;
        end
    end

    all_com_values(id_index,end) = sum(psth_indexes(:,1) == id);

    [this_cells_oscillation_phase,oscillation_order] = sort(this_cells_oscillation_phase);

    this_cell_cross_corr    =  this_cell_cross_corr(oscillation_order,:);
    this_cell_score         = this_cell_score(oscillation_order);

    this_cell_cross_corr_zscored = this_cell_cross_corr;
    for j = 1:size(this_cell_cross_corr,1)
        this_cell_cross_corr_zscored(j,:) = zscore(this_cell_cross_corr_zscored(j,:));
        this_cell_cross_corr_zscored(j,:)   = movmean ( this_cell_cross_corr_zscored(j,:) , .05/cross_corr_bin);
    end

    figure('units','normalized','outerposition',[0 0 .25 1])
    subplot(6,1,1:3)
    imagesc(cross_time,1:size(this_cell_cross_corr,1), this_cell_cross_corr_zscored)
    axis xy
    clim([-2 2])

    subplot(6,1,5:6)
    %
    % plot(this_cell_score,this_cells_oscillation_phase, 'k.')
    % hold on
    % plot(this_cell_score,this_cells_oscillation_phase+2*pi, 'r.')
    positive_scores = this_cells_oscillation_phase(this_cell_score>0);
    negative_scores = this_cells_oscillation_phase(this_cell_score<0);
    polarhistogram(positive_scores,0:(pi/16):(2*pi), 'FaceColor','r', 'FaceAlpha',.5, 'Normalization',  'percentage')
    hold on
    polarhistogram(negative_scores, 0:(pi/16):(2*pi),'FaceColor','b', 'FaceAlpha',.5, 'Normalization',  'percentage')

    if numel(positive_scores)>1 &&  numel(negative_scores)>1
        central_os     = circ_rtest(positive_scores);
        central_nonos  = circ_rtest(negative_scores);
        [pval, table] = circ_wwtest(positive_scores, negative_scores);
        title(num2str([central_os central_nonos pval table{2,5}]))
    end
    pause(.1)
end

%%
j = 2
figure
subplot(5,1,1:2)
% positive_scores     = oscillation_corr(oscillation_corr(:,6+j)>0,3);
% negative_scores     = oscillation_corr(oscillation_corr(:,6+j)<0,3);

positive_scores = all_com_values(:,17 );
positive_scores = positive_scores(~isnan(positive_scores));
negative_scores = all_com_values(:,18 );
negative_scores = negative_scores(~isnan(negative_scores));
positive_mean = circ_mean(positive_scores);
negative_mean = circ_mean(negative_scores);
central_os     = circ_rtest(positive_scores);
central_nonos  = circ_rtest(negative_scores);
[pval, table] = circ_wwtest(positive_scores, negative_scores);



polarhistogram(positive_scores,0:(pi/16):(2*pi), 'FaceColor','r', 'FaceAlpha',.5, 'Normalization',  'percentage')
hold on
polarhistogram(negative_scores, 0:(pi/16):(2*pi),'FaceColor','b', 'FaceAlpha',.5, 'Normalization',  'percentage')
title(pval)
polar_axis = gca;
r_range = polar_axis.RLim;
polarplot([positive_mean positive_mean],r_range, 'r', 'LineWidth',2 )
polarplot([negative_mean negative_mean],r_range, 'b', 'LineWidth',2 )


subplot(5,1,4:5)

histogram(log10(all_com_values(all_com_values(:,end)>=10, 15)), -10:.5:0, 'EdgeColor','none')
hold on
histogram(log10(all_com_values(all_com_values(:,end)>=10, 21)), -10:.5:0, 'EdgeColor','none')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% From now on behavior analysis
id_list = unique(ALL_TABLES.Id);

rate_data_bin_size = .025;
rate_data_time = 0:rate_data_bin_size:ceil(max(ALL_TABLES.RunEndTime));
rate_data_count = zeros(size(rate_data_time));
rate_data_id    = zeros(numel(id_list),numel(rate_data_time));
rate_data_phase  = nan(numel(id_list),numel(rate_data_time));
rate_data_amp    = nan(numel(id_list),numel(rate_data_time));

for j=1:size(ALL_TABLES,1)

    rate_data_count(rate_data_time>=ALL_TABLES.RunStartTime(j) & rate_data_time<=ALL_TABLES.RunEndTime(j)) = ...
        rate_data_count(rate_data_time>=ALL_TABLES.RunStartTime(j) & rate_data_time<=ALL_TABLES.RunEndTime(j)) +1;

    if ~isnan(ALL_TABLES.RunStartTime(j)) & ALL_TABLES.SinFit(j,6)>0.25

        rate_data_id(id_list==ALL_TABLES.Id(j),...
            rate_data_time>=ALL_TABLES.RunStartTime(j) & rate_data_time<=ALL_TABLES.RunEndTime(j))=1;

        phase_shift = mod(ALL_TABLES.SinFit(j,3),2*pi);

        amplitud    = ALL_TABLES.SinFit(j,1);
        phase_shift = phase_shift + (amplitud<0)*pi;
        base_rate   = ALL_TABLES.SinFit(j,5);
        phase_shift = mod(phase_shift, 2*pi);
        rate_data_phase(id_list==ALL_TABLES.Id(j),...
            rate_data_time>=ALL_TABLES.RunStartTime(j) & rate_data_time<=ALL_TABLES.RunEndTime(j))=phase_shift;
        rate_data_amp(id_list==ALL_TABLES.Id(j),...
            rate_data_time>=ALL_TABLES.RunStartTime(j) & rate_data_time<=ALL_TABLES.RunEndTime(j))= abs(amplitud./base_rate);
    end
end

raylight_time = nan(size(rate_data_phase,2),1);

for ti = 1:numel(raylight_time)
    current_angles = rate_data_phase(:,ti);
    current_angles = current_angles(~isnan(current_angles));

    if numel(current_angles)>1
        raylight_time(ti) = circ_rtest(current_angles);
    end
end

rate_data_time_matrix = repmat(rate_data_time, size(rate_data_phase, 1),1);
complex_angle_osc_cells = angle(exp(1i*rate_data_phase(rate_data_amp>.8)));
complex_angle_nonosc_cells = angle(exp(1i*rate_data_phase(rate_data_amp<.6)));

mean_angle_per_time_osc_cells = angle(complex_angle_osc_cells);
mean_angle_per_time_nonosc_cells = angle(complex_angle_nonosc_cells);
mean_angle_length   = abs(complex_angle);
mean_amp            = nanmean(rate_data_amp);

%%


session = 1;
behavior_names = {'CC'};
beh_index = ismember(Behavior.Animal,'Single') & ismember(Behavior.Type2, behavior_names);
behaviors2plot = [Behavior.Start(beh_index) Behavior.End(beh_index)];
behaviors2plot = behaviors2plot(behaviors2plot(:,1)>=animal_sessions(session,1) & ...
    behaviors2plot(:,2)<=animal_sessions(session,2),:);

behavior_state_vector = false(numel(rate_data_time),1);

for bn = 1:size(behaviors2plot,1)
    beh_start = behaviors2plot(bn,1);
    beh_end = behaviors2plot(bn,2);
    index = rate_data_time>=beh_start & rate_data_time<=beh_end;
    behavior_state_vector(index)=1;
end

[sorted_channels, channel_order] = sort(chanel_list);
figure
imagesc(rate_data_time,1:numel(chanel_list),  rate_data_amp(channel_order,:))
colormap(1-gray)
hold on
[~,pag_bottom] = min(abs(sorted_channels - 20));
[~, pag_upper] = min(abs(sorted_channels - 310));
for bn = 1:size(behaviors2plot,1)
    beh_start = behaviors2plot(bn,1);
    beh_end = behaviors2plot(bn,2);
    fill([beh_start beh_end beh_end beh_start], [0 0 numel(chanel_list) numel(chanel_list)], 'r', 'FaceAlpha',.5, 'EdgeColor','none')
end
hold on
plot([rate_data_time(1) rate_data_time(end)],[pag_bottom pag_bottom], 'b')
plot([rate_data_time(1) rate_data_time(end)],[pag_upper pag_upper], 'b')
plot([animal_sessions(session,1) animal_sessions(session,1)],[0 numel(chanel_list) ], 'b')
plot([animal_sessions(session,2) animal_sessions(session,2)],[0 numel(chanel_list) ], 'b')

xticks([animal_sessions(session,1) mean(animal_sessions(session,:)) mean(animal_sessions(session,:))+5 animal_sessions(session,2)])
xticklabels({'Beg Session','', '5 sec', ['End Session (', num2str(round(range(animal_sessions(session,:)))), ' s)']})
yticks([pag_bottom pag_upper ])
yticklabels({'PAG LOWER LIMIT', 'PAG UPPER LIMIT'})
clim([0 1])
axis xy

hold on
set(gca, 'FontSize', 18)
title(behavior_names)

%% for trainign SVM
data_range4training = rate_data_time>=animal_sessions(session,1) & ...
    rate_data_time<=animal_sessions(session,2);
trasposed_matrix_neuron = rate_data_amp';
trasposed_matrix_neuron(isnan(trasposed_matrix_neuron)) =0;
trasposed_matrix_neuron = trasposed_matrix_neuron(data_range4training,:);

% trasposed_matrix_neuron(isnan(trasposed_matrix_neuron)) =0;
% addedn_noise = randn(size(trasposed_matrix_neuron))*std(trasposed_matrix_neuron(:))/5;
%
% trasposed_matrix_neuron(trasposed_matrix_neuron==0) = addedn_noise(trasposed_matrix_neuron==0);


deleted_columns = find(sum(trasposed_matrix_neuron)==0);
trasposed_matrix_neuron(:,deleted_columns) = [];
trasposed_behavior      = behavior_state_vector(data_range4training);
%
% test_values= randsample(size(trasposed_matrix_neuron,1), round(.1*size(trasposed_matrix_neuron,1)), false);
% data_for_testing        = trasposed_matrix_neuron(test_values,:);
% behavior_for_testing    =  trasposed_behavior(test_values);
%
% trasposed_matrix_neuron(test_values,:) = [];
% trasposed_behavior(test_values) = [];
% SVMModel = fitcsvm(trasposed_matrix_neuron,trasposed_behavior,'Standardize',true,'KernelFunction','RBF',...
%     'KernelScale','auto');
fscnca_mdl = fscnca(trasposed_matrix_neuron,trasposed_behavior,'Solver','sgd');
GLMModel = fitglm(trasposed_matrix_neuron,trasposed_behavior, 'Distribution', 'binomial');
% c = cvpartition(n,"KFold",k)
%
% fcn = @(Xtr, Ytr, Xte) predict(...
%     GeneralizedLinearModel.fit(Xtr,Ytr,'linear','distr','normal'), ...
%     Xte);
%
% % perform cross-validation, and return average MSE across folds
% mse = crossval('mcr', trasposed_matrix_neuron, trasposed_behavior), 'Predfun',fcn, 'kfold',10);

% compute root mean squared error
% avrg_rmse = sqrt(mse)

%%
channels2plot = chanel_list;
channels2plot(deleted_columns) = [];
figure
subplot(1,3,1)
swarmchart(fscnca_mdl.FeatureWeights,25*round(channels2plot/25)  + 12.5, 'k.','YJitter','density')
ch2average = unique(25*round(channels2plot/25));
mean_values = ch2average*0;
for j=1:numel(ch2average)
    mean_values(j) = mean(fscnca_mdl.FeatureWeights(25*round(channels2plot/25) ==ch2average(j)));
end
hold on
plot(mean_values,ch2average  + 12.5, 'r')
hold on
x_lim = xlim;
plot(x_lim,[20 20], 'b')
plot(x_lim,[310 310], 'b')
ylim([0 350])
xlabel('Neuron Weight')
title('NCA')




subplot(1,3,2)
glm_coeff = -log10(GLMModel.Coefficients.SE(2:end));
swarmchart(glm_coeff, 25*round(channels2plot/25)  + 12.5,'k.','YJitter','density')
ch2average = unique(25*round(channels2plot/25));
mean_values = ch2average*0;
for j=1:numel(ch2average)
    mean_values(j) = mean(glm_coeff(25*round(channels2plot/25) ==ch2average(j)));
end
hold on
plot(mean_values ,ch2average+ 12.5, 'r')
x_lim = xlim;
plot(x_lim,[20 20], 'b')
plot(x_lim,[310 310], 'b')
ylim([0 350])
title('GLM')
xlabel('Neuron SE')

subplot(1,3,3)
swarmchart(abs(d_primes), 25*round(channels2plot/25)  + 12.5,'k.','YJitter','density')
ch2average = unique(25*round(channels2plot/25));
mean_values = ch2average*0;

for j=1:numel(ch2average)
    mean_values(j) = mean(abs(d_primes(25*round(channels2plot/25) ==ch2average(j))));
end
hold on
plot(mean_values ,ch2average+ 12.5, 'r')
x_lim = xlim;
plot(x_lim,[20 20], 'b')
plot(x_lim,[310 310], 'b')

ylim([0 350])
title('Naive Bayes')
xlabel("Neuron D'")



%%

figure
subplot(3,1,1)
plot(fscnca_mdl.FeatureWeights, -log10(GLMModel.Coefficients.pValue(2:end)), '.')
[c,p] = corr(fscnca_mdl.FeatureWeights, -log10(GLMModel.Coefficients.pValue(2:end)));
title([c, p])

subplot(3,1,2)
plot(fscnca_mdl.FeatureWeights, abs(d_primes), '.')
[c,p] = corr(fscnca_mdl.FeatureWeights, abs(d_primes));
title([c, p])

subplot(3,1,3)
plot(-log10(GLMModel.Coefficients.pValue(2:end)), abs(d_primes), '.')
[c,p] = corr(-log10(GLMModel.Coefficients.pValue(2:end)), abs(d_primes));
title([c, p])

%% svm classification

data_range4training = rate_data_time>=animal_sessions(session,1) & ...
    rate_data_time<=animal_sessions(session,2);
trasposed_matrix_neuron = rate_data_amp';
trasposed_matrix_neuron(isnan(trasposed_matrix_neuron)) =0;
trasposed_matrix_neuron = trasposed_matrix_neuron(data_range4training,:);

deleted_columns = find(sum(trasposed_matrix_neuron)==0);
trasposed_matrix_neuron(:,deleted_columns) = [];
trasposed_behavior      = behavior_state_vector(data_range4training);

K = 10;
cv = cvpartition(numel(trasposed_behavior), 'kfold',K);

mse_svm = zeros(K,4);
for k=1:K
    % training/testing indices for this fold
    trainIdx = cv.training(k);
    testIdx = cv.test(k);

    % train GLM model
    mdl =  fitcsvm(trasposed_matrix_neuron,trasposed_behavior,'Standardize',true,'KernelFunction','RBF',...
        'KernelScale','auto');
    % mdl = GeneralizedLinearModel.fit(trasposed_matrix_neuron(trainIdx,:), trasposed_behavior(trainIdx), 'Distribution', 'binomial');

    % predict regression output
    Y_hat = predict(mdl, trasposed_matrix_neuron(testIdx,:));

    % Y_hat = Y_hat>.5;
    % compute mean squared error
    mse_svm(k,1) = sum(trasposed_behavior(testIdx)==1 &  Y_hat>.5)/sum(trasposed_behavior(testIdx)==1);
    mse_svm(k,2) = sum(trasposed_behavior(testIdx)==0 &  Y_hat>.5)/sum(trasposed_behavior(testIdx)==0);
    mse_svm(k,3) = sum(trasposed_behavior(testIdx)==1 &  Y_hat==0)/sum(trasposed_behavior(testIdx)==1);
    mse_svm(k,4) = sum(trasposed_behavior(testIdx)==0 &  Y_hat==0)/sum(trasposed_behavior(testIdx)==0);
end

% average RMSE across k-folds
figure
boxplot(100*mse_svm)
xticklabels({'TP','FP','FN','TN'})
set(gca, 'FontSize', 18)
ylabel('%')
title('SVM classification')
%% for training Naiv Bayes

data_range4training = rate_data_time>=animal_sessions(session,1) & ...
    rate_data_time<=animal_sessions(session,2);
trasposed_matrix_neuron = rate_data_amp';
trasposed_matrix_neuron(isnan(trasposed_matrix_neuron)) =0;
addedn_noise = randn(size(trasposed_matrix_neuron))*std(trasposed_matrix_neuron(:))/5;

trasposed_matrix_neuron(trasposed_matrix_neuron==0) = addedn_noise(trasposed_matrix_neuron==0);
trasposed_matrix_neuron = trasposed_matrix_neuron(data_range4training,:);
trasposed_behavior      = behavior_state_vector(data_range4training);


deleted_columns = find(std(trasposed_matrix_neuron(trasposed_behavior==1,:))==0 | std(trasposed_matrix_neuron(trasposed_behavior==0,:))==0 );
trasposed_matrix_neuron(:,deleted_columns) = [];

Mdl = fitcnb(trasposed_matrix_neuron,trasposed_behavior);

K = 25;
cv = cvpartition(numel(trasposed_behavior), 'kfold',K);

mse_nb = nan(K,4);
for k=1:K
    % training/testing indices for this fold
    trainIdx = cv.training(k);
    testIdx = cv.test(k);

    % train GLM model
    X_train = trasposed_matrix_neuron(trainIdx,:);
    y_train = trasposed_behavior(trainIdx);
    if ~any(std(X_train(y_train==1,:))==0) && ~any(std(X_train(y_train==0,:))==0)
        mdl =  fitcnb(X_train,y_train);
        % mdl = GeneralizedLinearModel.fit(trasposed_matrix_neuron(trainIdx,:), trasposed_behavior(trainIdx), 'Distribution', 'binomial');

        % predict regression output
        Y_hat = predict(mdl, trasposed_matrix_neuron(testIdx,:));

        % Y_hat = Y_hat>.5;
        % compute mean squared error
        mse_nb(k,1) = sum(trasposed_behavior(testIdx)==1 &  Y_hat>.5)/sum(trasposed_behavior(testIdx)==1);
        mse_nb(k,2) = sum(trasposed_behavior(testIdx)==0 &  Y_hat>.5)/sum(trasposed_behavior(testIdx)==0);
        mse_nb(k,3) = sum(trasposed_behavior(testIdx)==1 &  Y_hat==0)/sum(trasposed_behavior(testIdx)==1);
        mse_nb(k,4) = sum(trasposed_behavior(testIdx)==0 &  Y_hat==0)/sum(trasposed_behavior(testIdx)==0);
    end
end

% average RMSE across k-folds
figure
boxplot(100*mse_nb)
xticklabels({'TP','FP','FN','TN'})
set(gca, 'FontSize', 18)
ylabel('%')
title('NB classification')


%
% test_values= randsample(size(trasposed_matrix_neuron,1), round(.1*size(trasposed_matrix_neuron,1)), false);
% data_for_testing        = trasposed_matrix_neuron(test_values,:);
% behavior_for_testing    =  trasposed_behavior(test_values);
%
% trasposed_matrix_neuron(test_values,:) = [];
% trasposed_behavior(test_values) = [];

%%
animal_types
session = 1;
behavior_names = {'CC', 'CB'};
beh_index = ismember(Behavior.Animal,'Single') & ismember(Behavior.Type2, behavior_names);
behaviors2plot = [Behavior.Start(beh_index) Behavior.End(beh_index)];
behaviors2plot = behaviors2plot(behaviors2plot(:,1)>=animal_sessions(session,1) & ...
    behaviors2plot(:,2)<=animal_sessions(session,2),:);


figure
subplot(3,1,1)
plot(rate_data_time,log10(raylight_time), 'k')
hold on
for bn = 1:size(behaviors2plot,1)
    beh_start = behaviors2plot(bn,1);
    beh_end = behaviors2plot(bn,2);
    fill([beh_start beh_end beh_end beh_start], [0 0 -10 -10], 'r', 'FaceAlpha',.25, 'EdgeColor','none')
end
plot([rate_data_time(1) rate_data_time(end)], log10(0.05)*[1 1], 'g')
% plot(rate_data_time, sum(rate_data_id))

subplot(3,1,2)
plot(rate_data_time, mean_angle_per_time, '.k')
hold on
for bn = 1:size(behaviors2plot,1)
    beh_start = behaviors2plot(bn,1);
    beh_end = behaviors2plot(bn,2);
    fill([beh_start beh_end beh_end beh_start], [-pi -pi pi pi], 'r', 'FaceAlpha',.25, 'EdgeColor','none')
end
% plot([behaviors2plot(:,1)';behaviors2plot(:,1)'], [-pi pi], 'r')

subplot(3,1,3)
plot(rate_data_time, mean_angle_length , 'k')
hold on
hold on
for bn = 1:size(behaviors2plot,1)
    beh_start = behaviors2plot(bn,1);
    beh_end = behaviors2plot(bn,2);
    fill([beh_start beh_end beh_end beh_start], [0 0 1 1], 'r', 'FaceAlpha',.25, 'EdgeColor','none')
end




evoked_response_range = [-2 2];
evoked_bin_size = rate_data_bin_size;
evoked_bin_deges = evoked_response_range(1):evoked_bin_size:evoked_response_range(2);

recruitment_response    = nan(size(behaviors2plot,1), range(evoked_response_range)/evoked_bin_size + 1);
angular_response        = nan(size(behaviors2plot,1), range(evoked_response_range)/evoked_bin_size + 1);
length_response         = nan(size(behaviors2plot,1), range(evoked_response_range)/evoked_bin_size + 1);
amplitud_response       = nan(size(behaviors2plot,1), range(evoked_response_range)/evoked_bin_size + 1);
for bn = 1:size(behaviors2plot,1)
    beh_start = behaviors2plot(bn,1);
    [~,angle_index] = min(abs(rate_data_time-beh_start));
    index2select = (angle_index + evoked_response_range(1)/evoked_bin_size):(angle_index + evoked_response_range(2)/evoked_bin_size);
    if min(index2select)>0 && max(index2select)<=numel(mean_angle_per_time)
        angular_response(bn,:) = mean_angle_per_time(index2select);
        length_response(bn,:) = mean_angle_length(index2select);
        recruitment_response(bn,:) = rate_data_count(index2select);
        amplitud_response(bn,:) = mean_amp(index2select);
    end
end

baseline_correction = amplitud_response;
baseline_index = evoked_bin_deges<-.5;
for j =1:size(baseline_correction,1)

    baseline_correction(j,:) = (baseline_correction(j,:) - mean(baseline_correction(j,baseline_index)))/std(baseline_correction(j,baseline_index));
end


figure
subplot(5,3,1:3:10)
matrix2plot = recruitment_response(~any(isnan(recruitment_response),2),:);
imagesc(evoked_bin_deges, 1:size(matrix2plot,1), matrix2plot)
title(behavior_names)
mean_response = mean(matrix2plot);

subplot(5,3,13)
ci =1.96*std(matrix2plot)/sqrt(size(matrix2plot,1));
fill([evoked_bin_deges fliplr(evoked_bin_deges)],[mean_response-ci fliplr(mean_response+ci)],'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
hold on
plot(evoked_bin_deges, mean_response, 'k')


subplot(5,3,2:3:11)
matrix2plot = length_response(~any(isnan(length_response),2),:);
imagesc(evoked_bin_deges, 1:size(matrix2plot,1), matrix2plot)

subplot(5,3,14)
mean_response = mean(matrix2plot);
ci =1.96*std(matrix2plot)/sqrt(size(matrix2plot,1));
fill([evoked_bin_deges fliplr(evoked_bin_deges)],[mean_response-ci fliplr(mean_response+ci)],'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
hold on
plot(evoked_bin_deges, mean_response, 'k')



subplot(5,3,3:3:12)
matrix2plot = baseline_correction(~any(isnan(baseline_correction),2),:);


imagesc(evoked_bin_deges, 1:size(matrix2plot,1), matrix2plot)
title(behavior_names)
mean_response = mean(matrix2plot);
subplot(5,3,15)


ci =1.96*std(matrix2plot)/sqrt(size(matrix2plot,1));
fill([evoked_bin_deges fliplr(evoked_bin_deges)],[mean_response-ci fliplr(mean_response+ci)],'k', 'FaceAlpha',.5, 'EdgeColor','none', 'HandleVisibility','off')
hold on
plot(evoked_bin_deges, mean_response, 'k')

% plot([behaviors2plot';behaviors2plot'], [0 1], 'r')

% t = rate_data_time(50);
%
% relevant_ranges = ALL_TABLES.RunStartTime<=t & ALL_TABLES.RunEndTime>=t;

%%
treshold = 2;
find_events
[L,n] = bwlabeln(all_rate);

nn = 1;

group_indexes = find(L==nn);
max_rate = max(all_rate(group_indexes));

% selected_indexes = all_rate(group_indexes)>.

figure

plot(this_neuron_spikes,this_neuron_phases, 'k.')
hold on
% plot(this_neuron_spikes,this_neuron_phases+2*pi, 'k.')


yyaxis right
plot(all_rate_time, all_rate)
%%
% function pred = classf(X1train,ytrain,X1test)
% Xtrain = table(X1train,ytrain);
% Xtest = table(X1test);
% mdl = fitglm(Xtrain,ytrain);
% yfit = predict(mdl,Xtest);
% pred = (yfit > 0.5);
% end