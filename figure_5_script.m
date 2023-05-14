%% EPI of real data

%% Parameters
min_perc_samples = 0.1; % Min num of samples per OTU
N = 1000; % Total number of top OTU
threshold_min_pres = 0.2; % Minimum num of samples per OTU - FOR ANALYSIS
threshold_max_pres = 0.8; % Maximum num of samples per OTU - FOR ANALYSIS
threshold_net = 0.2; % Percentile threshold for network construction of Q
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

%% Shuffling data
% hmp.data = shuffle_data(hmp.data);

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
    
    % Renomralize:
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
    
    % Renormalize:
    data_pres = data_pres./sum(data_pres);
    data_abs = data_abs./sum(data_abs);
    
%     dist_i = pdist2(data_pres', data_abs', @BC_pdist2);
%     hmp.dist_i{i} = dist_i;
    hmp.dist_i{i} = pdist2(data_i', data_i', @BC_pdist2);
    
    D1(i) = sum(sum(pdist2(data_pres', data_abs', @BC_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));
    D2(i) = pdist2(mean(data_pres'), mean(data_abs'), @BC_pdist2);    
    
end

%% Degree
% % % Creating co-occurence network -------------------------------------
% D = 1 - pdist(data, 'spearman');
% D(isnan(D)) = 0;
% 
% % Calculating p-values with reshuffling -------------------------------
% D_shuffle = zeros(num_reshuffle_degree, N*(N-1)/2);
% for m = 1:num_reshuffle_degree
%     m
%     data_shuffle = shuffle_data(hmp.data);
%     D_shuffle(m, :) = 1 - pdist(data_shuffle, 'spearman');
% end
% z = (D - mean(D_shuffle))./std(D_shuffle);
% pvals = 2*normcdf(-abs(z));
% 
% % Multiple comparison using Benjamini & Hochberg ----------------------
% q = 0.05;
% method = 'pdep';
% report = 'yes';
% [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report);
% D = D.*h;
% D = squareform(D);
% D(1:N+1:end) = 1;
% 
% % Converting nan to 0 SHOULD I DO THIS????? ---------------------------
% D(isnan(D)) = 0;
% 
% % Calculating topology -------------------------------------------------
% G_D = graph(D);
% degree = G_D.degree;
% %     betweenness_centrality = centrality(G_D, 'betweenness');
% %     closenness_centrality = centrality(G_D, 'closeness');
% 
% hmp.degree = degree;

hmp.Q = Q;
hmp.D1 = D1;
hmp.D2 = D2;
hmp.sp = sum(hmp.datad, 2)/num_samples;


%% Network based NMI
% hmp.MI = pdist_MI(hmp.datad); % Similarity
% temp = hmp.MI; temp = 1 - temp; temp(1:N+1:end) = 0; % Distance
% hmp.MI_dist = temp;
% 
% % Calculating PCoA based on MI
% % MI_y = cmdscale(hmp.MI_dist);
% MI_y = tsne(hmp.datad, 'Distance', @tsne_MI);
% hmp.MI_y = MI_y;
% 
% % Setting MI threshold for network
% MI_vector = squareform(hmp.MI_dist, 'tovector');
% [cdf_dist, dist] = ecdf(MI_vector);
% MI_thresh = dist(find(cdf_dist<threshold_MI, 1, 'last')); 
% hmp.MI_vector = MI_vector;
% hmp.MI_thresh = MI_thresh;
% 
% % Building MI network
% MI_mat = squareform(MI_vector, 'tomatrix')<=MI_thresh;
% MI_mat(1:N+1:end) = 0;
% MI_graph = graph(MI_mat);
% hmp.MI_mat = MI_mat;
% hmp.MI_graph = MI_graph;
% 
% % Correlation
% hmp.MI_corr = corr(hmp.datad', hmp.datad');
% 
% % Positive indecies
% temp = hmp.MI_corr.*(hmp.MI_dist<MI_thresh);
% temp(1:N+1:end) = 0;
% temp = temp>0;
% [node_start,node_end] = find(temp);
% hmp.node_start_pos = node_start;
% hmp.node_end_pos = node_end;
% % highlight(h_MI,node_start,node_end,'EdgeColor',bcolor);
% 
% % Negative indecies
% temp = hmp.MI_corr.*(hmp.MI_dist<MI_thresh);
% temp(1:N+1:end) = 0;
% temp = temp<0;
% [node_start,node_end] = find(temp);
% hmp.node_start_neg = node_start;
% hmp.node_end_neg = node_end;
% % highlight(h_MI,node_start,node_end,'EdgeColor',rcolor);
% % highlight(s,node_start,node_end,'EdgeColor',[1,0,0]);

%% PCoA of top modularity/D1/D2 within rel_freq @@@@@@@@@@@@@@@@@@@@@@@@@@@
figure_5 = figure;
letterposx = -0.25;
letterposy = 1.2;
set(gcf, 'Position', [20 100 1000 1000], 'Units', 'centimeters');
% Top Q -------------------------------------------------------------------
subplot(3, 3, 8); box on;
[~, ind] = nanmax(hmp.Q);
ind_max_Q = ind;
B_i_graph = hmp.B_i_graph{ind};
s_i = hmp.s_i{ind};

% Plotting network
h1 = plot(B_i_graph, 'NodeLabel', [], 'EdgeColor', [0 0 0], 'EdgeAlpha', 1, 'Layout', 'force');
hold on;
axis square;

xvalues = h1.XData;
yvalues = h1.YData;

% Removing nodes
temp = find(s_i==1);
highlight(h1, temp, 'NodeColor', bcolor, 'MarkerSize', 30, 'Marker', 'none');
temp = find(s_i==-1);
highlight(h1, temp, 'NodeColor', bcolor, 'MarkerSize', 30, 'Marker', 'none');

% Redrawing nodes
sz = 30;
temp = find(s_i==-1);
scatter(xvalues(temp), yvalues(temp), sz, bcolor, 'MarkerFaceColor', [1 1 1]);
temp = find(s_i==1);
scatter(xvalues(temp), yvalues(temp), sz, bcolor, 'filled');


% Highlighting edges
% edges = table2array(B_i_graph.Edges);
% for i = 1:size(edges, 1)
%     temp1 = s_i(edges(i, 1));
%     temp2 = s_i(edges(i, 2));
%     if temp1==1 && temp2 == 1 % Blue-Blue
%         highlight(h1, edges(i, 1), edges(i, 2), 'EdgeColor', bcolor, 'LineWidth', linewidth);
%     elseif temp1==-1 && temp2==-1 % Red-Red
%         highlight(h1, edges(i, 1), edges(i, 2), 'EdgeColor', rcolor, 'LineWidth', linewidth);
%     else % Blue-Red/Red-Blue
%         highlight(h1, edges(i, 1), edges(i, 2), 'EdgeColor', kcolor, 'LineWidth', linewidth);
%     end
% end
text(letterposx, letterposy, 'h', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');

% Adding legend
% h = zeros(2, 1);
% h(1) = plot(NaN,NaN,'.','Color',node_color_parm.*bcolor);
% h(2) = plot(NaN,NaN,'.','Color',node_color_parm.*rcolor);
% legend(h, 'Present', 'Absent');

set(gca,'XTick',[]);
set(gca,'YTick',[]);
text(0.5, -0.1 ,'Sample-to-sample network, B^i', 'Units', 'normalized', 'FontSize', 13, 'HorizontalAlignment', 'center');
% box on;
axis off;
set(gca, 'FontSize', fontsize);


% Top D1 ------------------------------------------------------------------
subplot(3, 3, 2); box off;
[~, ind] = nanmax(hmp.D1);
ind_max_D1 = ind;
ind_pres = hmp.data(ind, :) ~= 0;
ind_abs = not(ind_pres);
data_pres = hmp.data(:, ind_pres);
data_abs = hmp.data(:, ind_abs);
% Deleting species
data_pres(ind, :) = [];
data_abs(ind, :) = [];
% Calculating distances between all samples
dist_i = hmp.dist_i{ind};
% dist_i(isnan(dist_i)) = 0;
y = mdscale(dist_i, 2);
hold on;
% plot(y(ind_pres, 1), y(ind_pres, 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*bcolor);
% plot(y(not(ind_pres), 1), y(not(ind_pres), 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*rcolor);
sz = 30;
scatter(y(ind_pres, 1), y(ind_pres, 2), sz, bcolor, 'filled');
scatter(y(not(ind_pres), 1), y(not(ind_pres), 2), sz, bcolor);
xlabel('PC1'); ylabel('PC2');
set(gca, 'FontSize', fontsize);
set(gca, 'XTick', []); set(gca, 'YTick', []);
text(letterposx, letterposy, 'b', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');

l(1) = scatter(nan,nan, sz, [0 0 0], 'filled');
l(2) = scatter(nan,nan, sz, [0 0 0]);
legend(l, 'Present', 'Absent');

text(0.5, 1.1, 'Candidate keystone', 'Units', 'normalized', 'FontSize', 18, 'HorizontalAlignment', 'center');

% legend('Present', 'Absent', 'FontSize', 13);
axis square;


% Top D2 ------------------------------------------------------------------
subplot(3, 3, 5); box off;
[~, ind] = nanmax(hmp.D2);
ind_max_D2 = ind;
ind_pres = hmp.data(ind, :) ~= 0;
ind_abs = not(ind_pres);
data_pres = hmp.data(:, ind_pres);
data_abs = hmp.data(:, ind_abs);
% Deleting species
data_pres(ind, :) = [];
data_abs(ind, :) = [];
% Calculating distances between all samples
dist_i = hmp.dist_i{ind};
y = mdscale(dist_i, 2);
hold on;
% plot(y(ind_pres, 1), y(ind_pres, 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*bcolor);
% plot(y(not(ind_pres), 1), y(not(ind_pres), 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*rcolor);
sz = 30;
scatter(y(ind_pres, 1), y(ind_pres, 2), sz, bcolor, 'filled');
scatter(y(not(ind_pres), 1), y(not(ind_pres), 2), sz, bcolor);

y_pres_avg = [mean(y(ind_pres, 1)), mean(y(ind_pres, 2))];
y_abs_avg = [mean(y(not(ind_pres), 1)), mean(y(not(ind_pres), 2))];
plot(y_pres_avg(1), y_pres_avg(2), 'X', 'MarkerSize', 15, 'Color', 'k', 'LineWidth', 2);
plot(y_abs_avg(1), y_abs_avg(2), 'X', 'MarkerSize', 15, 'Color', 'k', 'LineWidth', 2)
plot([y_pres_avg(1), y_abs_avg(1)], [y_pres_avg(2), y_abs_avg(2)], 'k', 'LineWidth', 1.5)

xlabel('PC1'); ylabel('PC2');
set(gca, 'FontSize', fontsize);
set(gca, 'XTick', []); set(gca, 'YTick', []);
text(letterposx, letterposy, 'e', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
% legend('Present', 'Absent', 'FontSize', 13);
axis square;








% Rand Q ------------------------------------------------------------------
subplot(3, 3, 9); box on;
ind = randi([1, length(ind_rel)], 1);
ind = ind_rel(ind);
ind_rand_Q = ind;
B_i_graph = hmp.B_i_graph{ind};
s_i = hmp.s_i{ind};

% Plotting network
h1 = plot(B_i_graph, 'NodeLabel', [], 'EdgeColor', [0 0 0], 'EdgeAlpha', 1, 'Layout', 'force');
hold on;
axis square;

xvalues = h1.XData;
yvalues = h1.YData;

% Removing nodes
temp = find(s_i==1);
highlight(h1, temp, 'NodeColor', rcolor, 'MarkerSize', 30, 'Marker', 'none');
temp = find(s_i==-1);
highlight(h1, temp, 'NodeColor', rcolor, 'MarkerSize', 30, 'Marker', 'none');


% Highlighting edges
% edges = table2array(B_i_graph.Edges);
% for i = 1:size(edges, 1)
%     temp1 = s_i(edges(i, 1));
%     temp2 = s_i(edges(i, 2));
%     if temp1==1 && temp2 == 1 % Blue-Blue
%         highlight(h1, edges(i, 1), edges(i, 2), 'EdgeColor', bcolor, 'LineWidth', linewidth);
%     elseif temp1==-1 && temp2==-1 % Red-Red
%         highlight(h1, edges(i, 1), edges(i, 2), 'EdgeColor', rcolor, 'LineWidth', linewidth);
%     else % Blue-Red/Red-Blue
%         highlight(h1, edges(i, 1), edges(i, 2), 'EdgeColor', kcolor, 'LineWidth', linewidth);
%     end
% end
text(0, 1.1, 'i', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');

% Redrawing nodes
sz = 30;
temp = find(s_i==-1);
scatter(xvalues(temp), yvalues(temp), sz, rcolor, 'MarkerFaceColor', [1 1 1]);
temp = find(s_i==1);
scatter(xvalues(temp), yvalues(temp), sz, rcolor, 'filled');

% Adding legend
% h = zeros(2, 1);
% h(1) = plot(NaN,NaN,'.','Color',node_color_parm.*bcolor);
% h(2) = plot(NaN,NaN,'.','Color',node_color_parm.*rcolor);
% legend(h, 'Present', 'Absent');

set(gca,'XTick',[]);
set(gca,'YTick',[]);
text(0.5, -0.1 ,'Sample-to-sample network, B^j', 'Units', 'normalized', 'FontSize', 13, 'HorizontalAlignment', 'center');
box off;
axis off;
set(gca, 'FontSize', fontsize);


% Rand D1 -----------------------------------------------------------------
subplot(3, 3, 3); box off;
ind = randi([1, length(ind_rel)], 1);
ind = ind_rel(ind);
ind_rand_D1 = ind;
ind_pres = hmp.data(ind, :) ~= 0;
ind_abs = not(ind_pres);
data_pres = hmp.data(:, ind_pres);
data_abs = hmp.data(:, ind_abs);
% Deleting species
data_pres(ind, :) = [];
data_abs(ind, :) = [];
% Calculating distances between all samples
dist_i = hmp.dist_i{ind};
y = mdscale(dist_i, 2);
hold on;
% plot(y(ind_pres, 1), y(ind_pres, 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*bcolor);
% plot(y(not(ind_pres), 1), y(not(ind_pres), 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*rcolor);
sz = 30;
scatter(y(ind_pres, 1), y(ind_pres, 2), sz, rcolor, 'filled');
scatter(y(not(ind_pres), 1), y(not(ind_pres), 2), sz, rcolor);
xlabel('PC1'); ylabel('PC2');
set(gca, 'FontSize', fontsize);
set(gca, 'XTick', []); set(gca, 'YTick', []);
text(letterposx, letterposy, 'c', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
% legend('Present', 'Absent', 'FontSize', 13);
text(0.5, 1.1, 'Random species', 'Units', 'normalized', 'FontSize', 18, 'HorizontalAlignment', 'center');
axis square;


% Rand D2 -----------------------------------------------------------------
subplot(3, 3, 6); box off;
ind = randi([1, length(ind_rel)], 1);
ind = ind_rel(ind);
ind_rand_D2 = ind;
ind_pres = hmp.data(ind, :) ~= 0;
ind_abs = not(ind_pres);
data_pres = hmp.data(:, ind_pres);
data_abs = hmp.data(:, ind_abs);
% Deleting species
data_pres(ind, :) = [];
data_abs(ind, :) = [];
% Calculating distances between all samples
dist_i = hmp.dist_i{ind};
y = mdscale(dist_i, 2);
hold on;
% plot(y(ind_pres, 1), y(ind_pres, 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*bcolor);
% plot(y(not(ind_pres), 1), y(not(ind_pres), 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*rcolor);
sz = 30;
scatter(y(ind_pres, 1), y(ind_pres, 2), sz, rcolor, 'filled');
scatter(y(not(ind_pres), 1), y(not(ind_pres), 2), sz, rcolor);

y_pres_avg = [mean(y(ind_pres, 1)), mean(y(ind_pres, 2))];
y_abs_avg = [mean(y(not(ind_pres), 1)), mean(y(not(ind_pres), 2))];
plot(y_pres_avg(1), y_pres_avg(2), 'X', 'MarkerSize', 15, 'Color', 'k', 'LineWidth', 2);
plot(y_abs_avg(1), y_abs_avg(2), 'X', 'MarkerSize', 15, 'Color', 'k', 'LineWidth', 2)
plot([y_pres_avg(1), y_abs_avg(1)], [y_pres_avg(2), y_abs_avg(2)], 'k', 'LineWidth', 1.5)

xlabel('PC1'); ylabel('PC2');
set(gca, 'FontSize', fontsize);
set(gca, 'XTick', []); set(gca, 'YTick', []);
text(letterposx, letterposy, 'f', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
% legend('Present', 'Absent', 'FontSize', 13);
axis square;




% Histogram Q -------------------------------------------------------------
% Should I do bar plot instead? @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subplot(3, 3, 7);
% Calculating Qi for all species
h1 = histogram(hmp.Q, num_bins, 'EdgeColor', 'none', 'FaceColor', hist_color);
hold on;
ylimits = get(gca, 'Ylim'); 
ymin = ylimits(1); ymax = ylimits(2);
xlimits = get(gca, 'Xlim');
xmin = xlimits(1); xmax = xlimits(2);
% plot([D_key, D_key], [0, ymax/2]);
Annotate(gca, 'arrow', [xmin + (xmax-xmin)*(2/3), hmp.Q(ind_max_Q)], [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', bcolor);
% text(xmin + (xmax-xmin)*(2/3), ymin + (ymax-ymin)/2 + 1, 'Keystone species',...
%     'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.8], 'FontWeight', 'bold')
Annotate(gca, 'arrow', [xmin + (xmax-xmin)/3, hmp.Q(ind_rand_Q)],  [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', rcolor);
% text(xmin + (xmax-xmin)/3, ymin + (ymax-ymin)/2 + 5, 'Random species',...
%     'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.8], 'FontWeight', 'bold')
% annotation('arrow', [D_key, D_key], [0, ymax/2]);
xlabel('Q');
ylabel('Number of species');
set(gca, 'FontSize', fontsize);
text(letterposx, letterposy, 'g', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
% Adding color background
area([nanmean(hmp.Q) + 2*nanstd(hmp.Q) ,xmax+1], [ymax, ymax], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlim(xlimits);
box off;
axis square;
num_Q  = sum(hmp.Q(not(isnan(hmp.Q))) > (nanmean(hmp.Q) + 2*nanstd(hmp.Q)))./sum(not(isnan(hmp.Q)))

% Bar plot
% Ordering by Q
% [~, ind] = sort(hmp.Q, 'descend');
% bar(hmp.Q(ind), 'FaceColor', [0 0 0], 'FaceAlpha', 0.5);
% ylabel('$Q^i$', 'Interpreter', 'latex', 'Rotation', 0, 'HorizontalAlignment', 'right');
% xlabel('Species ranking', 'Interpreter', 'latex');
% set(gca, 'FontSize', fontsize);
% text(0, 1.1, '(a)', 'Units', 'normalized', 'FontSize', 18, 'Interpreter', 'latex');
% y_line_value = 2*std(hmp.Q) + mean(hmp.Q);
% hold on; plot([0, N], [y_line_value, y_line_value], 'Color', 'k', 'LineWidth', 2);
% xmin = 0; xmax = N; 
% ylimits = get(gca, 'Ylim'); 
% ymin = ylimits(1); ymax = ylimits(2);
% ind_max_Q_annot = find(ind==ind_max_Q);
% ind_rand_Q_annot = find(ind==ind_rand_Q);
% annot_const = 300;
% % Annotate(gca, 'arrow', [ind_max_Q_annot + annot_const, ind_max_Q_annot], [hmp.Q(ind_max_Q), hmp.Q(ind_max_Q)]);
% % Annotate(gca, 'arrow', [ind_rand_Q_annot + annot_const, ind_rand_Q_annot],  [hmp.Q(ind_rand_Q), hmp.Q(ind_rand_Q)], 'Linestyle', '--');
% axis square;

% Histogram D1 ------------------------------------------------------------
subplot(3, 3, 1);
% Calculating Qi for all species
h1 = histogram(hmp.D1, num_bins, 'EdgeColor', 'none', 'FaceColor', hist_color);
hold on;
ylimits = get(gca, 'Ylim'); 
ymin = ylimits(1); ymax = ylimits(2);
xlimits = get(gca, 'Xlim');
xmin = xlimits(1); xmax = xlimits(2);
% plot([D_key, D_key], [0, ymax/2]);
Annotate(gca, 'arrow', [xmin + (xmax-xmin)*(2/3), hmp.D1(ind_max_D1)], [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', bcolor);
text(xmin + (xmax-xmin)*(2/3), ymin + (ymax-ymin)/2 + 1, 'Candidate keystone',...
    'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.8], 'FontWeight', 'bold')
Annotate(gca, 'arrow', [xmin + (xmax-xmin)/3, hmp.D1(ind_rand_D1)],  [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', rcolor);
text(xmin + (xmax-xmin)/3, ymin + (ymax-ymin)/2 + 5, 'Random species',...
    'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.8], 'FontWeight', 'bold')
% annotation('arrow', [D_key, D_key], [0, ymax/2]);
xlabel('D_1');
ylabel('Number of species');
% ylabel('Number of species');
set(gca, 'FontSize', fontsize);
text(letterposx, letterposy, 'a', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
area([nanmean(hmp.D1) + 2*nanstd(hmp.D1) ,xmax+1], [40, 40], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlim(xlimits);
box off;
axis square;

num_D1 = sum(hmp.D1(not(isnan(hmp.D1))) > (nanmean(hmp.D1) + 2*nanstd(hmp.D1)))./sum(not(isnan(hmp.D1)))

% [~, ind] = sort(hmp.D1, 'descend');
% bar(hmp.D1(ind), 'FaceColor', [0 0 0], 'FaceAlpha', 0.5);
% ylabel('$D_1^i$', 'Interpreter', 'latex', 'Rotation', 0, 'HorizontalAlignment', 'right');
% xlabel('Species ranking', 'Interpreter', 'latex');
% set(gca, 'FontSize', fontsize);
% text(0, 1.1, '(a)', 'Units', 'normalized', 'FontSize', 18, 'Interpreter', 'latex');
% y_line_value = 2*std(hmp.D1) + mean(hmp.D1);
% hold on; plot([0, N], [y_line_value, y_line_value], 'Color', 'k', 'LineWidth', 2);
% xmin = 0; xmax = N; 
% ylimits = get(gca, 'Ylim'); 
% ymin = ylimits(1); ymax = ylimits(2);
% ind_max_Q_annot = find(ind==ind_max_D1);
% ind_rand_Q_annot = find(ind==ind_rand_D1);
% annot_const = 300;
% % Annotate(gca, 'arrow', [ind_max_Q_annot + annot_const, ind_max_Q_annot], [hmp.D1(ind_max_Q), hmp.D1(ind_max_Q)]);
% % Annotate(gca, 'arrow', [ind_rand_Q_annot + annot_const, ind_rand_Q_annot],  [hmp.D1(ind_rand_Q), hmp.D1(ind_rand_Q)], 'Linestyle', '--');
% % ylim([0.5 max(hmp.D1)]);
% axis square;

% Histogram D2 ------------------------------------------------------------
subplot(3, 3, 4);
% Calculating Qi for all species
h1 = histogram(hmp.D2, num_bins, 'EdgeColor', 'none', 'FaceColor', hist_color);
hold on;
ylimits = get(gca, 'Ylim'); 
ymin = ylimits(1); ymax = ylimits(2);
xlimits = get(gca, 'Xlim');
xmin = xlimits(1); xmax = xlimits(2);
% plot([D_key, D_key], [0, ymax/2]);
Annotate(gca, 'arrow', [xmin + (xmax-xmin)*(2/3), hmp.D2(ind_max_D2)], [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', bcolor);
Annotate(gca, 'arrow', [xmin + (xmax-xmin)/3, hmp.D2(ind_rand_D2)],  [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', rcolor);
% annotation('arrow', [D_key, D_key], [0, ymax/2]);
xlabel('D_2');
ylabel('Number of species');
% ylabel('Number of species');
set(gca, 'FontSize', fontsize);
text(letterposx, letterposy, 'd', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
area([nanmean(hmp.D2) + 2*nanstd(hmp.D2) ,xmax+1], [ymax, ymax], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlim(xlimits);
box off;
axis square;

num_D2 = sum(hmp.D2(not(isnan(hmp.D2))) > (nanmean(hmp.D2) + 2*nanstd(hmp.D2)))./sum(not(isnan(hmp.D2)))

% % axis square;
% [~, ind] = sort(hmp.D2, 'descend');
% bar(hmp.D2(ind), 'FaceColor', [0 0 0], 'FaceAlpha', 0.5);
% ylabel('$D_2^i$', 'Interpreter', 'latex', 'Rotation', 0, 'HorizontalAlignment', 'right');
% xlabel('Species ranking', 'Interpreter', 'latex');
% set(gca, 'FontSize', fontsize);
% text(0, 1.1, '(a)', 'Units', 'normalized', 'FontSize', 18, 'Interpreter', 'latex');
% y_line_value = 2*std(hmp.D2) + mean(hmp.D2);
% hold on; plot([0, N], [y_line_value, y_line_value], 'Color', 'k', 'LineWidth', 2);
% xmin = 0; xmax = N; 
% ylimits = get(gca, 'Ylim'); 
% ymin = ylimits(1); ymax = ylimits(2);
% ind_max_Q_annot = find(ind==ind_max_D2);
% ind_rand_Q_annot = find(ind==ind_rand_D2);
% annot_const = 300;
% % Annotate(gca, 'arrow', [ind_max_Q_annot + annot_const, ind_max_Q_annot], [hmp.D2(ind_max_Q), hmp.D2(ind_max_Q)]);
% % Annotate(gca, 'arrow', [ind_rand_Q_annot + annot_const, ind_rand_Q_annot],  [hmp.D2(ind_rand_Q), hmp.D2(ind_rand_Q)], 'Linestyle', '--');
% % ylim([0.5 max(hmp.D1)]);
% axis square;