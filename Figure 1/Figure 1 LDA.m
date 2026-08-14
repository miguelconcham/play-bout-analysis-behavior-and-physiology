%%  define folders and load main data
data_root                       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\Codes repository\Data';
hmm_data_folder                 = [data_root, '\HMM data\HMM raw data'];
labeled_data_folder             = [data_root, '\HMM data\locomotive behaviors'];
segmented_data_folder           = [data_root, '\Analysis results\locomotive behaviors 2 partners'];
not_labeled_data_folder         = [data_root, '\Analysis results\non labeled behavior 2 partners'];
folder_with_anima_names         = [data_root, '\Behavior backups'];
behavior_classification_folder  = [data_root, '\Analysis results\Behavior classification'];
dir_list = dir(fullfile(folder_with_anima_names, '*.txt'));
spatial_property_names      = {  'Speed','AngleSpeed','AngleAcc','Acc','Wall2CenterPos'...
    'RelativeDistance','RelativeSpeed','RelativeAngleSpeed','RelativeAngleAcc', 'RelativeAcc'};
call_prop_list = {'PrincipalFrequencykHz', 'SlopekHzs', 'Sinuosity', 'DeltaFreqkHz', 'FrequencyStandardDeviationkHz'};


behaviors2check = {'Pin', 'Boxing', 'Evasion', 'Pounce_A','Pounce_B','CD','Escape','CC','CB','Pounce_Ai','Pounce_Bi','Rearing','Sniffing', 'Bite', 'Scratch', 'Grooming'}; %% here you decide what behavior to extract


load([segmented_data_folder,'/all_behavior.mat'],'all_behavior')

load([segmented_data_folder,'/all_not_labeled_behavior.mat'],'all_not_labeled_behavior')
load([segmented_data_folder,'/all_var_names.mat'],'all_var_names')
%% Fig 1C estimate all properties, estiamte umap, color events and behavior and plot

disp('Creating properties matrix, where each raw is a behavior or a random event and the columns ar the motion-USV space')
all_properties = [ all_behavior(:,1);all_not_labeled_behavior(:,1)];
all_mean_prop = [];
for j=1:size(all_properties,1)
    this_var = all_properties{j,1};
    if size(this_var,1)==1
        all_mean_prop = [all_mean_prop;[this_var size(this_var,1) sum(this_var(:,1))]];
    else
        all_mean_prop = [all_mean_prop;[mean(this_var,'omitmissing') size(this_var,1) sum(this_var(:,1))]];
    end
    behavior_length(j) =size(this_var,1);
end
all_mean_prop(:, end-1) = zscore(all_mean_prop(:, end-1));
all_mean_prop(:, end) = (all_mean_prop(:, end) - mean(all_mean_prop(:, end), 'omitmissing'))/std(all_mean_prop(:, end), 'omitmissing');
all_spatial_prop = all_var_names(1:numel(all_var_names)-numel(call_prop_list)-1);
call_prop    = all_mean_prop(:, numel(all_spatial_prop)+1: numel(all_spatial_prop)+numel(call_prop_list)+1);
spatial_prop = all_mean_prop(:, 1:numel(all_spatial_prop));



[coeff_call,score_call,latent_call,tsquared_call,explained_call,mu_call]  = pca(call_prop);
[coeff_spatial,score_spatial,latent_spatial,tsquared_spatial,explained_spatial,mu_spatial]  = pca(spatial_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% create behavioral labels array %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('asigning a behavior label each raw of the matrix, it stay emtpy for random events')

behavior_labels = [all_behavior(:,2);repmat({''},size(all_not_labeled_behavior,1),1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% estimating  umap (using the mean of varaibles as features) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Estimating UMAP')

data2umap = all_mean_prop;
data2umap(:, ismember(all_var_names, 'WallPos')) = [];
nan_values = any(isnan(all_mean_prop),2);

data2umap = data2umap(~nan_values,:);
matched_labels = behavior_labels(~nan_values,:);



n_dim = 2;
% [reduction, umap] = run_umap(data2classify, 'n_components',2 );
[embedding_umaps, umapStruct] = run_umap(data2umap, ...
    'n_components',  n_dim, ...
    'n_neighbors', 25, ...
    'min_dist', 0.1, ...
    'metric', 'euclidean',...
    'verbose', 'none');



% play_var        = {'CC','Escape', 'CB'};
% aggresive_var   = {'Bite'};
% rearing_var     = {'Grooming','Scratch','Rearing'};
% 
% behavior_labels = (all_behavior(:,2));
% behavior_labels(ismember(behavior_labels, 'Pounce_B')) = {'Pounce_A'};
% behavior_labels(ismember(behavior_labels, 'Pounce_Bi')) = {'Pounce_Ai'};
% all_behavior_labels = behavior_labels;
% 
% 
% behavior_labels(ismember(behavior_labels,play_var))         = {'Play'};
% behavior_labels(ismember(behavior_labels,aggresive_var))    = {'Aggression'};
% behavior_labels(ismember(behavior_labels,rearing_var))      = {'Regulation'};
% 
% 
% behavior_motionUSV_features        = all_mean_prop(ismember(behavior_labels, {'Play','Aggression', 'Regulation'}),:);
% behavior_labels                = behavior_labels(ismember(behavior_labels, {'Play','Aggression', 'Regulation'}));
% 
% behavior_motionUSV_features = array2table(behavior_motionUSV_features);
% 
% n = size(behavior_motionUSV_features,2);
% colNames = arrayfun(@(i) sprintf('prop_%d', i), 1:n, 'UniformOutput', false);
% behavior_motionUSV_features.Properties.VariableNames = colNames;
% behavior_motionUSV_features(:,5) = []; % varaible 5 is wall distance that we havent used

% LDA_PLAY_AGG_REG classifier is a LDA trained on X = behavior_motionUSV_features and 
% y = behavior_labels 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% here i loaded the LDA that classify behavior as belongin to play   %%%%
%%%%%% aggression or regulation, as defined in the commented lines before %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading LDA classifier')
load([behavior_classification_folder,'\LDA_PLAY_AGG_REG.mat'])

disp('Assigning color to each event')
trained_model   = LDA_PLAY_AGG_REG.ClassificationDiscriminant;
Features        = all_mean_prop;   
Features(:,5)   = []; %removing wall pos
n               = size(Features,2);
Features        = array2table(Features);
colNames        = arrayfun(@(i) sprintf('mean_prop_%d', i), 1:n, 'UniformOutput', false);

Features.Properties.VariableNames = colNames;

Features        = Features(~nan_values,:);
umap_pca        = embedding_umaps;
% Assume you already have:
% trained_model: your trained LDA model
% Features: your feature matrix (N x d)
% umap_pca: your UMAP/PCA embedding (N x 2)
% trained_model.ClassNames: cell array of class names

% 1. Get LDA predictions and posterior probabilities
[predictedLabels, posteriorProbs] = predict(trained_model, Features);

% 2. Get max posterior probability for each observation (proximity to closest centroid)
[maxProb, ~] = max(posteriorProbs, [], 2);

% 3. Map categorical predicted class to numeric for colormap
[~, classIdx] = ismember(predictedLabels, trained_model.ClassNames);

% 4. Create a custom colormap for the 3 classes
colors = lines(numel(trained_model.ClassNames)); % distinct colors for each class\
% colors = hsv(numel(trained_model.ClassNames));
pointColors = zeros(size(classIdx,1),3);
for i = 1:numel(trained_model.ClassNames)
    % Blend class color with white based on confidence (maxProb)
    pointColors(classIdx==i,:) = (1-maxProb(classIdx==i))*[1 1 1] + maxProb(classIdx==i)*colors(i,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% and finally plot %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Ploting figure 1C, notice that it can be fliped in the x or y axis')

figure;

[~,re_order] = sort(posteriorProbs(:,3), 'descend');
scatter(umap_pca(re_order,1), umap_pca(re_order,2), 40, posteriorProbs(re_order,:)*colors, 'filled', 'MarkerEdgeColor','none','MarkerFaceAlpha',0.25);
title('Smoothed Posterior Map (Masked to Data Region)');
xlabel('UMAP/PCA 1'); ylabel('UMAP/PCA 2');

%% plot clusters of behavior sun embedign space (using umap or pca FIGURE 1 D)

behaviors2check = {'Pin' ,'Boxing' ,'Evasion' ,'Pounce_A' ,'Pounce_B','CD'  ,'Escape','CC'  ,'CB'  ,'Pounce_Ai' ,'Pounce_Bi' ,'Rearing' ,'Sniffing' ,'Bite' , 'Scratch' ,'Grooming',''}; %% here you decide what behavior to extract


% varaibles2cluster = num2cell(1:16); % tjis woul dbe ploting each behavior
% independetly
varaibles2cluster = {[10 11]}  % if you put two behavior indexes together
% in the same cell as in this line, then they will be plot together
cluster_cloors = repmat({[0 0 1]},1,numel(varaibles2cluster));
% dimension2project = [score_call(:,1), score_spatial(:,1)];
dimension2project= embedding_umaps;

pcts = [0.1:0.1:0.8];
  axiis_lim = [-13    7.75   -3.5    7.3]
for nc = 1:numel(varaibles2cluster)
  figure
  
    scatter(dimension2project(:,1), dimension2project(:,2), 20, '.k');
    hold on
    index = ismember(matched_labels,behaviors2check(varaibles2cluster{nc})) ;
    X = dimension2project(index,:);
    % Define grid over which to evaluate the density
    x1 = linspace(min(dimension2project(:,1))-1, max(dimension2project(:,1))+1, 100);
    x2 = linspace(min(X(:,2))-1, max(X(:,2))+1, 100);
    [xg, yg] = meshgrid(x1, x2);
    grid_points = [xg(:), yg(:)];

    % Perform 2D Kernel Density Estimation
    [f, xi] = ksdensity(X, grid_points);  % f is the density at each grid point
    f_grid = reshape(f, length(x2), length(x1)); % reshape to grid

    % Normalize to get cumulative density
    f_sorted = sort(f(:), 'descend');
    cdf = cumsum(f_sorted) / sum(f_sorted);

    % Define percentile levels (in density units)
    % desired percentiles
    levels = zeros(size(pcts));
    for i = 1:length(pcts)
        idx = find(cdf >= pcts(i), 1, 'first');
        levels(i) = f_sorted(idx);
    end
    % contour(x1, x2, f_grid, sort(levels), 'LineWidth', .1, 'EdgeColor',cluster_cloors{nc});

    sorted_levels = sort(levels);
    for j=fliplr(1:numel(levels)-1)
        contourf(x1, x2, f_grid, sorted_levels([j j+1]), 'LineWidth', .5, 'FaceColor',cluster_cloors{nc}, 'FaceAlpha',.4*j/numel(levels),  'EdgeColor','None');
    end
    % axis(axiis_lim)

    title(behaviors2check(varaibles2cluster{nc}))

        pause(.1)
end

%% Fig S1D (ploting behavior according to their score in the play-bite sapce)

if exist('LDA_PLAY_AGG_REG', 'var')~=1
    load([behavior_classification_folder,'\LDA_PLAY_AGG_REG.mat'])
end

trained_model = LDA_PLAY_AGG_REG.ClassificationDiscriminant;
Features = all_mean_prop(1:size(all_behavior,1),:);   
Features(:,5) = [];
n = size(Features,2);
Features = array2table(Features);
colNames = arrayfun(@(i) sprintf('mean_prop_%d', i), 1:n, 'UniformOutput', false);
Features.Properties.VariableNames = colNames;

% Features = Features(~nan_pca_variables,:);

% Assume you already have:
% trained_model: your trained LDA model
% Features: your feature matrix (N x d)
% umap_pca: your UMAP/PCA embedding (N x 2)
% trained_model.ClassNames: cell array of class names

% 1. Get LDA predictions and posterior probabilities
[predictedLabels, posteriorProbs] = predict(trained_model, Features);

classMeans = trained_model.Mu;          % C x D
Sigma = trained_model.Sigma;            % D x D (pooled covariance)
classes = trained_model.ClassNames;
nClasses = size(classMeans,1);

% Compute between-class scatter
overallMean = mean(classMeans,1);
B = zeros(size(Sigma));
for k = 1:nClasses
    diff = classMeans(k,:) - overallMean;
    B = B + (diff' * diff);
end

% Solve for discriminant directions (eigenvectors of inv(Sigma)*B)
[W,~] = eig(pinv(Sigma)*B);
W = real(W);   % remove tiny imaginary parts

% Sort by eigenvalues (descending)
[~, idx] = sort(diag(eig(pinv(Sigma)*B)), 'descend');
W = W(:, idx);

% Project features into this space
Xproj = (Features{:,:} - overallMean) * W;

% Get class centroids in projected space
projMeans = (classMeans - overallMean) * W;


i = 1; j = 2; % classes

v = projMeans(j,:) - projMeans(i,:);
scores = ((Xproj - projMeans(i,:)) * v') ./ (norm(v)^2);
scores = max(0,min(1,scores)); % clamp to [0,1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  now sort behavior according to the value in the projected score %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_behavior_relabeled = all_behavior(:,2);
all_behavior_relabeled(ismember(all_behavior_relabeled, 'Pounce_B')) = {'Pounce_A'};
all_behavior_relabeled(ismember(all_behavior_relabeled, 'Pounce_Bi')) = {'Pounce_Ai'};
behavior_list = unique(all_behavior(:,2))';
mean_val_per_category = nan(numel(behavior_list), 1);
behavior_class = nan(size(scores,1),1);

for bn = 1:numel(behavior_list)
    mean_val_per_category(bn) = mean(scores(ismember(all_behavior_relabeled, behavior_list{bn})));
    behavior_class(ismember(all_behavior_relabeled, behavior_list{bn})) = bn;
end
[sorted_means, scores_sorted] = sort(mean_val_per_category, 'desc');

new_order = zeros(size(scores_sorted));
new_order(scores_sorted) = 1:numel(behavior_list);
behavior_class_sorted = new_order(behavior_class);

% 4. Re-map behavior_class to this new order


behavior_labels = all_behavior_relabeled;
x_pos = behavior_class_sorted + (rand(size(behavior_class_sorted))-.5)*.5;

%%%%%%%%%%%%%%%%%%
%%%% now plot %%%%
%%%%%%%%%%%%%%%%%%

figure
subplot(1,2,1)
swarmchart(behavior_class_sorted,scores, '.k')
hold on
plot(1:numel(sorted_means), sorted_means, '_r', 'MarkerSize', 10)
xticks(1:numel(behavior_list))
xticklabels(behavior_list(scores_sorted))



[~,re_order] = sort(scores);
scatter(x_pos(re_order), scores(re_order), 15, scores(re_order)*[1 0 0] + (1-scores(re_order))*[0 0 1], 'filled', 'MarkerFaceAlpha',1);
hold on
scatter(1:numel(sorted_means), sorted_means, 'ko', 'filled')
plot(1:numel(sorted_means), sorted_means, 'k', 'LineWidth',2)
xticks(1:numel(behavior_list))
xticklabels(behavior_list(scores_sorted))
axis tight

subplot(1,2,2)
[~,re_order] = sort(posteriorProbs(:,3), 'descend');
scatter(x_pos(re_order), scores(re_order), 15, posteriorProbs(re_order,:)*colors, 'filled', 'MarkerFaceAlpha',.5);
hold on
% plot(sort([(1:(numel(sorted_means)))-.5,(1:numel(sorted_means))+.5 ]), sorted_means(sort([1:numel(sorted_means),1:numel(sorted_means) ])), 'k');
% plot([(1:(numel(sorted_means)))-.5;(1:numel(sorted_means))-.5 ], [0*sorted_means'; sorted_means'], 'k')
% plot([(1:(numel(sorted_means)))+.5;(1:numel(sorted_means))+.5 ], [0*sorted_means'; sorted_means'], 'k')
scatter(1:numel(sorted_means), sorted_means, 'ko', 'filled')
plot(1:numel(sorted_means), sorted_means, 'k', 'LineWidth',2)
xticks(1:numel(behavior_list))
xticklabels(behavior_list(scores_sorted))
xlim tight
%% Fig S1B
load([behavior_classification_folder,'\trainedModel_LDA_ALL_BEHAVIORS.mat'],'trainedModel_LDA_ALL_BEHAVIORS') 

trained_model = trainedModel_LDA_ALL_BEHAVIORS.ClassificationDiscriminant;

cvMdl = trained_model.crossval;

k = cvMdl.KFold;
class_names =  trained_model .ClassNames;
cm_folds = nan(k,numel(class_names),numel(class_names));
sim    = cm_folds;
 n = numel(class_names);
for i = 1:k
    testIdx = test(cvMdl.Partition, i);
    preds = predict(cvMdl.Trained{i}, cvMdl.X(testIdx,:));
    CM = confusionmat(cvMdl.Y(testIdx), preds);    
    cm_folds(i,:,:) = CM;  


end
mean_cm = squeeze(mean(cm_folds));

CM_row = mean_cm ./ sum(mean_cm,2);
CM_col = mean_cm ./ sum(mean_cm,1);

S = (CM_col + CM_row') / 2;
S = (S + S') / 2;
S(1:n+1:end) = 0;

D = 1 - S/100;

D = (D + D') / 2;           % enforce symmetry (important!)
D(1:n+1:end) = 0;           % enforce zero diagonal
D(isnan(D)) = 1;            % handle any NaNs (max distance)

Z = linkage(D, 'ward');  % or 'complete', 'single', 'ward', etc.
n_clust=2;
clusters = cluster(Z, 'maxclust', n_clust);
labels = cellstr(class_names);

clust2 = find(clusters==2);
clust1 = find(clusters==1);
re_order = [clust1;clust2]';

row_sums = sum(mean_cm, 2);
CM_perc = mean_cm ./ row_sums;   % fraction per column
CM_perc = CM_perc * 100;

figure
colormap(1-gray)
CM_perc(1:n+1:end) = 0;
imagesc(1:numel(re_order),1:numel(re_order),CM_perc(re_order,re_order))
xticks(1:numel(re_order))
yticks(1:numel(re_order))

xticklabels(class_names(re_order))
yticklabels(class_names(re_order))
hold on
for j=1:n

    this_class_conf = CM_perc(re_order,re_order(j));
    [top_pctg, top_2] = sort(this_class_conf, 'descend');
    top_2 = top_2(1:2);
    top_pctg = top_pctg(1:2);
    for k=1:2
        text(j-.3,top_2(k),num2str(round(top_pctg(k),0)),'Color','r')
    end
end
axis square

plot([numel(clust1) numel(clust1)]+.5, [0 numel(clust1)+.5], 'r')
plot([0 numel(clust1)+.5], [numel(clust1) numel(clust1)]+.5, 'r')

plot([numel(clust1) numel(clust1)]+.5, [0 numel(clust2)+.5]+numel(clust1), ':b')
plot([0 numel(clust2)+.5]+numel(clust1), [numel(clust1) numel(clust1)]+.5, ':b')


colorbar

%% Fig S1C
% trained LDA classificator to classify each behavior using mean properties
% of the behaviors, EXCLUDING tigmotaxis.
trained_model = trainedModel_LDA_ALL_BEHAVIORS.ClassificationDiscriminant;
centroids = trained_model.Mu;  % size: [numClasses x numFeatures]
covMat = trained_model.Sigma;  % shared covariance matrix

% Compute pairwise Mahalanobis distances between centroids
numClasses = size(centroids,1);
D = zeros(numClasses);
for i = 1:numClasses
    for j = 1:numClasses
        diff = centroids(i,:) - centroids(j,:);
        D(i,j) = sqrt(diff / covMat * diff');  % Mahalanobis distance
    end
end

% Use MDS to project these distances into 2D
Ymds = mdscale(D, 2, 'Criterion', 'metricstress');

% Plot the behaviors in this 2D space
figure;
scatter(Ymds(:,1), Ymds(:,2), 80, 'filled');
text(Ymds(:,1)+0.02, Ymds(:,2), trained_model.ClassNames, 'FontSize', 12);
title('Behaviors in LDA-centroid space (MDS of Mahalanobis distances)');
xlabel('Dimension 1');
ylabel('Dimension 2');
axis equal;
axis square