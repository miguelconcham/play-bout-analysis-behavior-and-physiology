%% super synch check



mean_pctl =1- squeeze(mean(pctl_spikes_histogram));
mean_pctl(mean_pctl==.5)= NaN;

figure
histogram(mean_pctl(:), 0:0.01:1, 'Normalization','percentage')
super_synch_pairs = find(sum(mean_pctl,2)>.1*size(mean_pctl,2));
numel(super_synch_pairs)
coincidence_ssp = squeeze(mean(synch_spikes_histogram(:, super_synch_pairs,:)));


figure

plot(time_centers,mean(coincidence_ssp, 'omitmissing'), 'r')
hold on
plot(time_centers,squeeze(mean(mean(synch_spikes_histogram), 'omitmissing')), 'k')






 all_area_comb_ssp  = all_area_comb(super_synch_pairs,:);

 [Mu,ia,ic] = unique(cell2table(all_area_comb_ssp), 'rows', 'stable');           % Unique Values By Row, Retaining Original Order
h = accumarray(ic, 1);                              % Count Occurrences
Mu.count = h;

[all_Mu, all_ia, all_ic] = unique(cell2table(all_area_comb), 'rows', 'stable'); 
h = accumarray(all_ic, 1); 
Mu.TotalCount =- nan(size(Mu,1),1);
for j=1:size(Mu,1)
    Mu.TotalCount(j) = h(ismember(all_Mu.all_area_comb1,Mu.all_area_comb_ssp1(j)) & ismember(all_Mu.all_area_comb2,Mu.all_area_comb_ssp2(j)));
end
h = accumarray(all_ia, 1);     

data = string(table2array(Mu(:, {'all_area_comb_ssp1','all_area_comb_ssp2'})));
sortedRows = sort(data, 2);

% 2. Get a unique ID for every unique combination
[~, ~, ic] = unique(sortedRows, 'rows');

% 3. Group the original row indexes based on those unique IDs
% 'rowGroups' will be a cell array where each cell contains the indexes for a combination
rowGroups = accumarray(ic, (1:size(data,1))', [], @(x) {sort(x')});

% 4. Filter to keep only the groups with more than 1 index (the mirrored pairs)
mirroredPairs = rowGroups(cellfun(@(x) length(x) > 1, rowGroups));

% Display the pairs
disp('Mirrored row index pairs [m, n]:');
celldisp(mirroredPairs);

for j=1:size(mirroredPairs,1)

    Mu.count(mirroredPairs{j}) = repmat(sum(Mu.count(mirroredPairs{j}) ),2,1);

     Mu.TotalCount(mirroredPairs{j}) = repmat(sum(Mu.TotalCount(mirroredPairs{j}) ),2,1);
end

Mu.Percentage = 100*Mu.count./Mu.TotalCount;
session_id_ssp     = session_id(super_synch_pairs);
 clusters_id_ssp    = clusters_id(super_synch_pairs,:);

 %%
cc_indexes = nan(size(session_id_ssp,1),1);
for j=1:size(session_id_ssp,1)

   putative_index =  find(ismember(session_id_cc,session_id_ssp(j)) & all(ismember(clusters_id_cc,clusters_id_ssp(j,:)),2));

   if ~isempty(putative_index)
       cc_indexes(j) = putative_index;
   end
end

%%
smoth_wind = 1;
cc_time= (psth_edges_cc(1:end-1)+psth_edges_cc(2:end))/2;
this_cc = squeeze(play_cc(:,:));
for k = 1:size(this_cc,1)
    this_cc(k,:) = movmean( this_cc(k,:),smoth_wind);
end

this_cc_control = all_cross_corr_play_pctl(:,:);
for k = 1:size(this_cc_control,1)
    this_cc_control(k,:) = movmean( this_cc_control(k,:),smoth_wind);
end

this_cc_zs          = this_cc;
this_cc_zs_within = this_cc_zs;
this_cc_control_zs  = this_cc_control;
for k = 1:size(this_cc,1)
    this_cc_zs_within(k,:) =(this_cc(k,:)-mean(this_cc_zs_within(k,abs(cc_time)>.1 & abs(cc_time)<.5)))/std(this_cc_zs_within(k,abs(cc_time)>.1 & abs(cc_time)<.5));
    this_cc_zs(k,:) =(this_cc(k,:)-mean(this_cc_control(k,:)))/std(this_cc_control(k,:));
    this_cc_control_zs(k,:) =(this_cc_control_zs(k,:)-mean(this_cc_control(k,:)))/std(this_cc_control(k,:));
end


this_cc = this_cc_zs-this_cc_control_zs;



%%

sharpness = mean(this_cc(:,abs(cc_time)<.02),2)./ mean(this_cc(:,abs(cc_time)>.1 & abs(cc_time)<.5),2);
prominenc = (mean(this_cc(:,abs(cc_time)<.02),2) - mean(this_cc(:,abs(cc_time)>.1 & abs(cc_time)<.5),2))./std(this_cc(:,abs(cc_time)>.1 & abs(cc_time)<.5),[],2);
%%
figure
sharp_tresh = 6;

subplot(2,2,1)
plot(cc_time,median(this_cc_zs_within(sharpness>sharp_tresh,:), 'omitmissing'), 'b')
hold on
plot(cc_time,median(this_cc_zs_within(sharpness<sharp_tresh,:), 'omitmissing'), ':b')
hold on

ssp_this_cc = this_cc_zs_within(cc_indexes,:);
plot(cc_time,median(ssp_this_cc(sharpness(cc_indexes)>sharp_tresh,:), 'omitmissing'), 'r')
plot(cc_time,median(ssp_this_cc(sharpness(cc_indexes)<sharp_tresh,:), 'omitmissing'), ':r')
xlim([-.5 .5])

subplot(2,2,2)

hold on
plot(cc_time,median(this_cc, 'omitmissing'), ':b')
hold on

plot(cc_time,median(this_cc(cc_indexes,:), 'omitmissing'), 'r')
xlim([-.5 .5])

subplot(2,2,3)
xi = -10:0.1:10;
% 1. Calculate the empirical CDF
[f, x] = ecdf(sharpness(cc_indexes));
% 2. Remove duplicate x-values (keeping the highest f for each)
[x_unique, idx] = unique(x, 'last');
f_unique = f(idx);
% 3. Interpolate safely
f_new = interp1(x_unique, f_unique, xi, 'previous', 'extrap');
plot(xi,f_new, 'r' )
hold on

all_sharpness =sharpness;
all_sharpness(cc_indexes) = [];
all_sharpness(isinf(all_sharpness)) = [];

% 1. Calculate the empirical CDF
[f, x] = ecdf(all_sharpness);
% 2. Remove duplicate x-values (keeping the highest f for each)
[x_unique, idx] = unique(x, 'last');
f_unique = f(idx);
% 3. Interpolate safely

f_new = interp1(x_unique, f_unique, xi, 'previous', 'extrap');
plot(xi,f_new, 'k' )





subplot(2,2,4)
xi = -10:0.1:10;
% 1. Calculate the empirical CDF
[f, x] = ecdf(prominenc(cc_indexes));
% 2. Remove duplicate x-values (keeping the highest f for each)
[x_unique, idx] = unique(x, 'last');
f_unique = f(idx);
% 3. Interpolate safely
f_new = interp1(x_unique, f_unique, xi, 'previous', 'extrap');
plot(xi,f_new, 'r' )
hold on

all_prominenc =prominenc;
all_prominenc(cc_indexes) = [];
all_prominenc(isinf(all_prominenc)) = [];

% 1. Calculate the empirical CDF
[f, x] = ecdf(all_prominenc);
% 2. Remove duplicate x-values (keeping the highest f for each)
[x_unique, idx] = unique(x, 'last');
f_unique = f(idx);
% 3. Interpolate safely

f_new = interp1(x_unique, f_unique, xi, 'previous', 'extrap');
plot(xi,f_new, 'k' )



%%
sharp_tresh = 5;
j =14;
 raw_n = 500;
 raw_n = 10;
x_lim = [-.5 .5];
entrained_group1 = trough;
entrained_group2 = trough;
bin_size_cc = mean(diff(time_centers));
area_combination_indexes       = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
                                      (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)));
     cell_type_comb_trough          = find(entrained_group1(idx_pairs(:,1)) &  entrained_group2(idx_pairs(:,2)) &  area_combination_indexes );








figure
trough_this_cc = this_cc_zs_within(cell_type_comb_trough,:);
subplot(1,4,1)
imagesc(cc_time, 1:size(ssp_this_cc,1), ssp_this_cc)
clim([-2 2])


subplot(1,4,2)
no_trough_this_cc = this_cc_zs_within;
[~, order] = sort(sharpness);

imagesc(cc_time, 1:size(no_trough_this_cc,1), no_trough_this_cc(order,:))


clim([-2 2])
subplot(1,4,3)
hold on


plot(cc_time,median(ssp_this_cc(sharpness(cell_type_comb_trough)>sharp_tresh,:), 'omitmissing'), 'r')
plot(cc_time,median(ssp_this_cc(sharpness(cell_type_comb_trough)<sharp_tresh,:), 'omitmissing'), ':r')

subplot(1,4,4)

xi = -10:0.1:10;
% 1. Calculate the empirical CDF
[f, x] = ecdf(sharpness(cell_type_comb_trough));
% 2. Remove duplicate x-values (keeping the highest f for each)
[x_unique, idx] = unique(x, 'last');
f_unique = f(idx);
% 3. Interpolate safely
f_new = interp1(x_unique, f_unique, xi, 'previous', 'extrap');
plot(xi,f_new, 'r' )
hold on

all_sharpness =sharpness;
all_sharpness(cell_type_comb_trough) = [];
all_sharpness(isinf(all_sharpness)) = [];

% 1. Calculate the empirical CDF
[f, x] = ecdf(all_sharpness);
% 2. Remove duplicate x-values (keeping the highest f for each)
[x_unique, idx] = unique(x, 'last');
f_unique = f(idx);
% 3. Interpolate safely

f_new = interp1(x_unique, f_unique, xi, 'previous', 'extrap');
plot(xi,f_new, 'k' )