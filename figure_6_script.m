%% NMI of real data

%% EPI of real data

%% Parameters
min_perc_samples = 0.1; % Min num of samples per OTU
N = 1000; % Total number of top OTU
threshold_min_pres = 0.2; % Minimum num of samples per OTU - FOR ANALYSIS
threshold_max_pres = 0.8; % Maximum num of samples per OTU - FOR ANALYSIS
threshold_net = 0.25; % Percentile threshold for network construction of Q
fontsize = 15;
num_reshuffling = 20; % Number of reshuffling for pvalue calculations D/Q


% Plot parameters
bcolor = [0  113.9850  188.9550]/255;
rcolor = [216.7500   82.8750   24.9900]/255;
kcolor = [0 0 0];
node_color_parm = 0.8;
linewidth = 1.5;
markersize = 3;
num_bins = 100;
threshold_MI = 0.01;
MI_marker_size = 5;
edgealpha = 0.5;
hist_color = [0.4, 0.4, 0.4];

%% Loading data
load HMP_stool_2.mat;
load hmp_plaque.mat
load HMP_throat.mat;
load HMP_stool_only1.mat
hmp.data = HMP_stool_2.data;
hmp.OTU = HMP_stool_2.OTU;
% hmp.data = HMP_throat.data;
% hmp.OTU = HMP_throat.OTU;
% hmp.data = hmp_plaque.data;
% hmp.OTU = hmp_plaque.OTU;
% hmp.data = HMP_stool_only1.data;
% hmp.OTU = HMP_stool_only1.OTU;

% hmp.data = hmp.data(:, 1:20);

%% Removing constant rows
[data,PS] = removeconstantrows(hmp.data);
hmp.data = data;
hmp.OTU(PS.remove,:) = [];

%% Normalizing the data
hmp.data = hmp.data./(sum(hmp.data));

%% Removing OTUs that are not in min_perc_samples of the samples
% [~, num_samples] = size(hmp.data);
% temp = find(sum(hmp.data>0,2)<round(num_samples*min_perc_samples));
% hmp.data(temp,:) = [];
% hmp.OTU(temp,:) = [];

%% Ordering by mean abundance
[~,temp] = sort(mean(hmp.data,2),'descend');
hmp.data = hmp.data(temp,:);
hmp.OTU = hmp.OTU(temp,:);

%% Choosing top N OTUs
hmp.data = hmp.data(1:N,:);
hmp.OTU = hmp.OTU(1:N,:);

%% Converting to presence/absence
hmp.datad = double(hmp.data>0);

%% Removing constant rows of presence/absence
% [datad,PS] = removeconstantrows(hmp.datad);
% hmp.datad = datad;
% hmp.data(PS.remove,:) = [];
% hmp.OTU(PS.remove,:) = [];

%% Removing repeated rows of presence/absance (NO NEED)
% [hmp.datad,ia] = unique(hmp.datad,'rows','stable');
% hmp.data = hmp.data(ia,:);
% hmp.OTU = hmp.OTU(ia,:);

%% Removing species with less then 30% presence/absance
% temp = sum(hmp.datad,2)>threshold_min_pres*num_samples & sum(hmp.datad,2)<threshold_max_pres*num_samples;
% hmp.datad = hmp.datad(temp,:);
% hmp.data = hmp.data(temp,:);
% hmp.OTU = hmp.OTU(temp,:);

% Updating N/num_samples
N = size(hmp.datad,1);
num_samples = size(hmp.datad,2);

%% Randomly reordering species
ind = randperm(N);
hmp.data = hmp.data(ind, :);
hmp.datad = hmp.datad(ind, :);
hmp.OTU = hmp.OTU(ind, :);

hmp.data_shuffled = shuffle_data(hmp.data);
hmp.datad_shuffled = hmp.data_shuffled > 0;

%% Calculating relative freq of all speices
hmp.rel_freq = sum(hmp.data>0, 2)/num_samples;

%% Calculating modularity/D1/D2/degree for all species within range of rel_freq
% Init
Q = nan(N, 1);
D1 = nan(N, 1);
D2 = nan(N, 1);
degree = nan(N, 1);
num_reshuffle_degree = 1000;

M = num_samples;
data = hmp.data;
datad = hmp.datad;

hmp.B_i_graph = cell(N, 1);
hmp.s_i = cell(N, 1);
hmp.dist_i = cell(N, 1);

% Indecies of rel_freq
ind_rel = find(hmp.rel_freq > threshold_min_pres & hmp.rel_freq < threshold_max_pres);

% Start
for ii = 1:length(ind_rel)
    i = ind_rel(ii);
    i
    % Real data
    data_i = data; data_i(i, :) = [];
    
    % Renormalizing
    data_i = data_i./sum(data_i);
    
    distances_i = pdist(data_i', @BC_pdist2);
    [cdf_dist, dist] = ecdf(distances_i);
    dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last')); % Finding distance threshold
    B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
    B_i(1:M+1:end) = 0;
    B_i_graph = graph(B_i);
    
    % Modularity 
    s_i = double(datad(i, :)); s_i(s_i==0) = -1; s_i = s_i';
    Q(i) = modularity_guy(B_i, s_i);
    
    hmp.B_i_graph{i} = B_i_graph;
    hmp.s_i{i} = s_i;
    
    % D1 & D2
    ind_pres = data(i, :) ~= 0;
    data_pres = data(:, ind_pres);
    data_abs = data(:, not(ind_pres));
    data_pres(i, :) = [];
    data_abs(i, :) = [];
    
    % Renormalizing
    data_pres = data_pres./sum(data_pres);
    data_abs = data_abs./sum(data_abs);
    
%     dist_i = pdist2(data_pres', data_abs', @BC_pdist2);
%     hmp.dist_i{i} = dist_i;
    hmp.dist_i{i} = pdist2(data_i', data_i', @BC_pdist2);
    
    D1(i) = sum(sum(pdist2(data_pres', data_abs', @BC_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));
    D2(i) = pdist2(mean(data_pres'), mean(data_abs'), @BC_pdist2);    
    
end

% Degree
% % Creating co-occurence network -------------------------------------
D = 1 - pdist(data, 'spearman');
D(isnan(D)) = 0;

% Calculating p-values with reshuffling -------------------------------
D_shuffle = zeros(num_reshuffle_degree, N*(N-1)/2);
for m = 1:num_reshuffle_degree
    m
    data_shuffle = shuffle_data(hmp.data);
    D_shuffle(m, :) = 1 - pdist(data_shuffle, 'spearman');
end
z = (D - mean(D_shuffle))./std(D_shuffle);
pvals = 2*normcdf(-abs(z));

% Multiple comparison using Benjamini & Hochberg ----------------------
q = 0.05;
method = 'pdep';
report = 'yes';
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report);
D = D.*h;
D = squareform(D);
D(1:N+1:end) = 1;

% Converting nan to 0 SHOULD I DO THIS????? ---------------------------
D(isnan(D)) = 0;

% Calculating topology -------------------------------------------------
G_D = graph(D);
degree = G_D.degree;
%     betweenness_centrality = centrality(G_D, 'betweenness');
%     closenness_centrality = centrality(G_D, 'closeness');

hmp.degree = degree;

hmp.Q = Q;
hmp.D1 = D1;
hmp.D2 = D2;
hmp.sp = sum(hmp.datad, 2)/num_samples;


%% Correlation between different measures
figure_6 = figure; 
set(gcf, 'Position', [20 100 1000 500], 'Units', 'centimeters');

fontsize = fontsize*2;
sz = 20;

% Q vs D1 -----------------------------------------------------------------
subplot(1, 3, 2);
scatter(hmp.D1, hmp.Q, sz, [0 0 0], 'filled');
ylabel('Q'); xlabel('D_1');
set(gca, 'FontSize', fontsize);
text(0, 1.1, 'b', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
% legend('Present', 'Absent', 'FontSize', 13);
axis square;
[r,p] = corr(hmp.D1, hmp.Q, 'rows', 'complete');
text(0.6, 0.15, "r = " + string(r) + newline + "p = " + string(p),...
    'Units', 'normalized', 'FontSize', 15,...
    'HorizontalAlignment', 'left',...
    'BackgroundColor', [1 1 1], 'EdgeColor', 'k');
set(gca, 'FontSize', fontsize);

% Q vs D2 -----------------------------------------------------------------
subplot(1, 3, 3);
scatter(hmp.D2, hmp.Q, sz, [0 0 0], 'filled');
ylabel('Q'); xlabel('D_2');
set(gca, 'FontSize', fontsize);
text(0, 1.1, 'c', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
% legend('Present', 'Absent', 'FontSize', 13);
axis square;
[r,p] = corr(hmp.D2, hmp.Q, 'rows', 'complete');
text(0.6, 0.15, "r = " + string(r) + newline + "p = " + string(p),...
    'Units', 'normalized', 'FontSize', 15,...
    'HorizontalAlignment', 'left',...
    'BackgroundColor', [1 1 1], 'EdgeColor', 'k');
set(gca, 'FontSize', fontsize);

% D1 vs D2 -----------------------------------------------------------------
subplot(1, 3, 1);
scatter(hmp.D2, hmp.D1, sz, [0 0 0], 'filled');
ylabel('D_1'); xlabel('D_2');
set(gca, 'FontSize', fontsize);
text(0, 1.1, 'a', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
% legend('Present', 'Absent', 'FontSize', 13);
axis square;
[r,p] = corr(hmp.D1, hmp.D2, 'rows', 'complete');
text(0.6, 0.15, "r = " + string(r) + newline + "p = " + string(p),...
    'Units', 'normalized', 'FontSize', 15,...
    'HorizontalAlignment', 'left',...
    'BackgroundColor', [1 1 1], 'EdgeColor', 'k');
set(gca, 'FontSize', fontsize);
