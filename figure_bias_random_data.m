%% Calculating modularity/distance of real data with bootstrap/consistency

%% Parameters
min_perc_samples = 0.1; % Min num of samples per OTU
N = 500; % Total number of top OTUs
threshold_min_pres = 0.2; % Minimum num of samples per OTU
threshold_max_pres = 0.8; % Maximum num of samples per OTU
threshold_net = 0.5; % Percentile threshold for network construction of Q
fontsize = 15;
num_samples = 100;

%% Creating the data
% Random data with missing values
% p = 0.5; % Average number of species in each sample
% hmp.data = rand(N, num_samples).*double(rand(N, num_samples)<=0.8);
% hmp.datad = hmp.data>0;

% Random data with missing values/number of missing unique for each sample
% p = 0.2 + (0.8-0.2).*rand(num_samples, 1);
% data = rand(N, num_samples);
% for i = 1:num_samples
%     data(:, i) = data(:, i).*double(rand(N, 1) <= p(i));
% end
% hmp.data = data;
% hmp.datad = hmp.datad > 0;

% Random correlated data with missing values unique to sample
% https://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor
% k = 1;
% rho = factor(N, k);
% rho(1:N+1:end) = 1;
% rho_without_diag = rho(:);
% rho_without_diag(1:N+1:end) = [];
% data = copularnd('Gaussian', rho, num_samples)';
% p = 0.2 + (0.8-0.2).*rand(num_samples, 1);
% for i = 1:num_samples
%     data(:, i) = data(:, i).*double(rand(N, 1) <= p(i));
% end
% hmp.data = data;
% hmp.datad = hmp.datad > 0;

% Random correlated data with missing values unique to species
% https://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor
% k = 1;
% rho = factor(N, k);
% rho(1:N+1:end) = 1;
% rho_without_diag = rho(:);
% rho_without_diag(1:N+1:end) = [];
% data = copularnd('Gaussian', rho, num_samples)';
% p = 0.2 + (0.8-0.2).*rand(N, 1);
% for i = 1:N
%     data(i, :) = data(i, :) .* double(rand(1, num_samples) >= p(i));
% end
% hmp.data = data;
% hmp.datad = hmp.data > 0;

% Random uncorrelated data with missing values unique to species
data = rand(N, num_samples);
p = 0.2 + (0.8-0.2).*rand(N, 1);
for i = 1:N
    data(i, :) = data(i, :) .* double(rand(1, num_samples) >= p(i));
end
hmp.data = data;
hmp.datad = hmp.data > 0;


%% Shuffled data
hmp.data_shuffled = shuffle_data(data);
hmp.datad_shuffled = shuffle_data(hmp.datad);

%% Calculating modularity of each species
% Using Bray-Curtis

% Real data:
Q_i = nan(N, 1);
data = hmp.data;
datad = hmp.datad;
M = num_samples;
for i = 1:N
    i
    % Creating network of samples B_i without species i
    data_i = data; data_i(i, :) = [];
%     distances_i = pdist(data_i', @rJSD_pdist2);
    distances_i = pdist(data_i', @BC_pdist2);
    [cdf_dist, dist] = ecdf(distances_i);
    dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last')); % Finding distance threshold
    B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
    B_i(1:M+1:end) = 0;
    B_i_graph = graph(B_i);
    
    % Calculating modularity
    s_i = double(datad(i, :)); s_i(s_i==0) = -1; s_i = s_i';
    Q_i(i) = modularity_guy(B_i, s_i);
%     s_i_temp = s_i; s_i_temp(s_i_temp==-1) = 2;
%     Q_i(i) = modularity(B_i, s_i_temp);
end

% Shuffled data:
Q_i_shuffled = nan(N, 1);
data_shuffled = hmp.data_shuffled;
datad_shuffled = hmp.datad_shuffled;
M = num_samples;
for i = 1:N
    i
    % Creating network of samples B_i without species i
    data_i = data_shuffled; data_i(i, :) = [];
%     distances_i = pdist(data_i', @rJSD_pdist2);
    distances_i = pdist(data_i', @BC_pdist2);
    [cdf_dist, dist] = ecdf(distances_i);
    dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last')); % Finding distance threshold
    B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
    B_i(1:M+1:end) = 0;
    B_i_graph = graph(B_i);
    
    % Calculating modularity
    s_i = double(datad_shuffled(i, :)); s_i(s_i==0) = -1; s_i = s_i';
    Q_i_shuffled(i) = modularity_guy(B_i, s_i);
%     s_i_temp = s_i; s_i_temp(s_i_temp==-1) = 2;
%     Q_i_shuffled(i) = modularity(B_i, s_i_temp);
end

%% Plotting networks of top and rand modularity
% figure; 
% edgecolor = 0.8*[1 1 1];
% % Top modularity
% subplot(2, 3, 1);
% [~, ind] = max(Q_i);
% % Creating the networks
% data = hmp.data;
% datad = hmp.datad;
% data_i = data; data_i(ind, :) = [];
% distances_i = pdist(data_i', @BC_pdist2);
% [cdf_dist, dist] = ecdf(distances_i);
% dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last')); % Finding distance threshold
% B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
% B_i(1:M+1:end) = 0;
% B_i_graph = graph(B_i);
% % Plotting graph and highlighting
% h1 = plot(B_i_graph, 'EdgeColor', edgecolor);
% ind_pres = find(datad(ind, :));
% ind_abs = find(not(datad(ind, :)));
% highlight(h1, ind_pres, 'NodeColor', [1 0 0]);
% 
% 
% % Rand modularity
% subplot(2, 3, 2);

%% Plotting histogram of modularity
% subplot(2, 3, 3); hold on;
% real_color = 0.4*[1 1 1];
% shuffled_color = 0.1*[1 1 1];
% num_bins = 50;
% histogram(Q_i, num_bins, 'EdgeColor', 'none', 'FaceColor', real_color);
% histogram(Q_i_shuffled, num_bins, 'EdgeColor', 'none', 'FaceColor', shuffled_color);
% legend('Real data', 'Shuffled data');
% xlabel('$Q^i$', 'Interpreter', 'latex');
% ylabel('Number of species');
% set(gca, 'FontSize', fontsize);
% axis square;

%% Calculating distances
% Real data:
D_i_1 = nan(N, 1);
D_i_2 = nan(N, 1);
data = hmp.data;
datad = hmp.datad;
M = num_samples;
for i = 1:N
    i
    ind_pres = data(i, :) ~= 0;
    data_pres = data(:, ind_pres);
    data_abs = data(:, not(ind_pres));
    data_pres(i, :) = [];
    data_abs(i, :) = [];
    
%     D_i_1(i) = sum(sum(pdist2(data_pres', data_abs', @rJSD_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));
%     D_i_2(i) = pdist2(mean(data_pres'), mean(data_abs'), @rJSD_pdist2);
    
    D_i_1(i) = sum(sum(pdist2(data_pres', data_abs', @BC_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));
    D_i_2(i) = pdist2(mean(data_pres'), mean(data_abs'), @BC_pdist2);
    
end

% Shuffled data:
D_i_shuffled_1 = nan(N, 1);
D_i_shuffled_2 = nan(N, 1);
data_shuffled = hmp.data_shuffled;
datad_shuffled = hmp.datad_shuffled;
M = num_samples;
for i = 1:N
    i
    ind_pres = data_shuffled(i, :) ~= 0;
    data_pres = data_shuffled(:, ind_pres);
    data_abs = data_shuffled(:, not(ind_pres));
    data_pres(i, :) = [];
    data_abs(i, :) = [];
    
%     D_i_shuffled_1(i) = sum(sum(pdist2(data_pres', data_abs', @rJSD_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));
%     D_i_shuffled_2(i) = pdist2(mean(data_pres'), mean(data_abs'), @rJSD_pdist2);
    
    D_i_shuffled_1(i) = sum(sum(pdist2(data_pres', data_abs', @BC_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));
    D_i_shuffled_2(i) = pdist2(mean(data_pres'), mean(data_abs'), @BC_pdist2);
end

%% Calculating degree
% num_reshuffling = 1000;
% 
% % Real data
% data = hmp.data;
% datad = hmp.datad;
% % Creating co-occurence network -------------------------------------------
% D = 1 - pdist(data, 'spearman');
% D(isnan(D)) = 0;
% 
% % Calculating p-values with reshuffling -----------------------------------
% D_shuffle = zeros(num_reshuffling, N*(N-1)/2);
% for m = 1:num_reshuffling
%     data_shuffle = shuffle_data(hmp.data);
%     D_shuffle(m, :) = 1 - pdist(data_shuffle, 'spearman');
% end
% z = (D - mean(D_shuffle))./std(D_shuffle);
% pvals = 2*normcdf(-abs(z));
% 
% % Multiple comparison using Benjamini & Hochberg --------------------------
% q = 0.05;
% method = 'pdep';
% report = 'yes';
% [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report);
% D = D.*h;
% D = squareform(D);
% D(1:N+1:end) = 1;
% 
% % Converting nan to 0 SHOULD I DO THIS????? -------------------------------
% D(isnan(D)) = 0;
% 
% % Calculating topology ----------------------------------------------------
% G_D = graph(D);
% degree_i = G_D.degree;
% betweenness_centrality = centrality(G_D, 'betweenness');
% closenness_centrality = centrality(G_D, 'closeness');
% 
% 
% 
% % Shuffled data
% data = hmp.data_shuffled;
% datad = hmp.datad_shuffled;
% % Creating co-occurence network -------------------------------------------
% D = 1 - pdist(data, 'spearman');
% D(isnan(D)) = 0;
% 
% % Calculating p-values with reshuffling -----------------------------------
% D_shuffle = zeros(num_reshuffling, N*(N-1)/2);
% for m = 1:num_reshuffling
%     data_shuffle = shuffle_data(hmp.data_shuffled);
%     D_shuffle(m, :) = 1 - pdist(data_shuffle, 'spearman');
% end
% z = (D - mean(D_shuffle))./std(D_shuffle);
% pvals = 2*normcdf(-abs(z));
% 
% % Multiple comparison using Benjamini & Hochberg --------------------------
% q = 0.05;
% method = 'pdep';
% report = 'yes';
% [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report);
% D = D.*h;
% D = squareform(D);
% D(1:N+1:end) = 1;
% 
% % Converting nan to 0 SHOULD I DO THIS????? -------------------------------
% D(isnan(D)) = 0;
% 
% % Calculating topology ----------------------------------------------------
% G_D = graph(D);
% degree_i_shuffled = G_D.degree;
% betweenness_centrality_shuffled = centrality(G_D, 'betweenness');
% closenness_centrality_shuffled = centrality(G_D, 'closeness');


%% PCoA Real
% subplot(2, 3, 4);

%% PCoA shuffled
% subplot(2, 3, 5);

%% Histogram D_i
% subplot(2, 3, 6); hold on;
% real_color = 0.4*[1 1 1];
% shuffled_color = 0.1*[1 1 1];
% num_bins = 50;
% histogram(D_i_1, num_bins, 'EdgeColor', 'none', 'FaceColor', real_color);
% histogram(D_i_shuffled_1, num_bins, 'EdgeColor', 'none', 'FaceColor', shuffled_color);
% legend('Real data', 'Shuffled data');
% xlabel('$D^i_{\mathrm{cross}}$', 'Interpreter', 'latex');
% ylabel('Number of species');
% set(gca, 'FontSize', fontsize);
% axis square;

%% Other plots
% Bias of average distance between samples as a function of freq. why?
% D_i_1 - Mean of distances
% D_i_2 - Distances of mean
figure; 
set(gcf, 'Position', [20 100 1500 500], 'Units', 'centimeters');
markersize = 20;
fontsize = 20;
freq = sum(hmp.data>0, 2)/num_samples;
freq_shuffled = sum(data_shuffled>0, 2)/num_samples;

subplot(1, 3, 1);
plot(freq, D_i_1, '.k', 'MarkerSize', markersize);
axis square; box on;
xlabel('Relative frequency'); 
ylabel('D_1');
text(0, 1.1, 'a', 'Units', 'normalized', 'FontSize', 30, 'FontWeight', 'Bold');
set(gca, 'FontSize', fontsize);

subplot(1, 3, 2);
plot(freq, D_i_2, '.k', 'MarkerSize', markersize);
axis square; box on;
xlabel('Relative frequency'); 
ylabel('D_2');
text(0, 1.1, 'b', 'Units', 'normalized', 'FontSize', 30, 'FontWeight', 'Bold');
set(gca, 'FontSize', fontsize);

subplot(1, 3, 3); 
plot(freq, Q_i, '.k', 'MarkerSize', markersize);
axis square; box on;
xlabel('Relative frequency'); 
ylabel('Q');
text(0, 1.1, 'c', 'Units', 'normalized', 'FontSize', 30, 'FontWeight', 'Bold');
set(gca, 'FontSize', fontsize);


% figure; 
% subplot(2, 4, 1); %
% plot(freq, D_i_1, '.', 'MarkerSize', markersize);
% xlabel('Freq'); ylabel('D^i_1'); 
% set(gca, 'FontSize', fontsize);
% 
% subplot(2, 4, 5);
% plot(freq_shuffled, D_i_shuffled_1, '.', 'MarkerSize', markersize);
% xlabel('Freq'); ylabel('D^i_1 Shuffled'); 
% set(gca, 'FontSize', fontsize);
% 
% subplot(2, 4, 2);
% plot(freq, D_i_2, '.', 'MarkerSize', markersize);
% xlabel('Freq'); ylabel('D^i_2'); 
% set(gca, 'FontSize', fontsize);
% 
% subplot(2, 4, 6);
% plot(freq_shuffled, D_i_shuffled_2, '.', 'MarkerSize', markersize);
% xlabel('Freq'); ylabel('D^i_2 Shuffled'); 
% set(gca, 'FontSize', fontsize);
% 
% subplot(2, 4, 3);
% plot(freq, Q_i, '.', 'MarkerSize', markersize);
% xlabel('Freq'); ylabel('Q^i'); 
% set(gca, 'FontSize', fontsize);
% 
% subplot(2, 4, 7);
% plot(freq_shuffled, Q_i_shuffled, '.', 'MarkerSize', markersize);
% xlabel('Freq'); ylabel('Q^i Shuffled'); 
% set(gca, 'FontSize', fontsize);
% 
% subplot(2, 4, 4);
% plot(freq, degree_i, '.', 'MarkerSize', markersize);
% xlabel('Freq'); ylabel('Degree^i');
% set(gca, 'FontSize', fontsize);
% 
% subplot(2, 4, 8);
% plot(freq_shuffled, degree_i_shuffled, '.', 'MarkerSize', markersize);
% xlabel('Freq shuffled'); ylabel('Degree^i shuffled');
% set(gca, 'FontSize', fontsize);

%% Saving
saving_png_pdf('D_Q_bias');

%% Function
function S = factor(d,k)
W = randn(d,k);
S = W*W' + diag(rand(1,d));
S = diag(1./sqrt(diag(S))) * S * diag(1./sqrt(diag(S)));
end
