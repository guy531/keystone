%% Correlation between keystonenness and modularity on simulations
% Using BC distance and ER network
% v2 - Distribution of keystone factors
% v3 - Distribution and structural

%% INTERACTION BASED @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%% Parameters
clear;
N = 100; % Number of species
num_samples = 100; % Number of samples
M = num_samples;
t0 = 0; % Start time
tf = 100; % End time
threshold = 1e-5; % Threshold for dead species\
threshold_net = 0.5; % Percentile threshold for network construction of Q
num_reshuffling = 20; % Number of reshuffling for pvalue calculations D/Q
hist_color = [0.4, 0.4, 0.4];
num_bins = 100;
fontsize = 15;
bcolor = [0  113.9850  188.9550]/255;
rcolor = [216.7500   82.8750   24.9900]/255;
kcolor = [0 0 0];
markersize = 10;
num_reshuffle_degree = 1000;
textfont = 12;
textxvalue = 0.4;
textyvalue = 0.15;

%% Creating networks
% Erdos-Renyi
k = ones(N, 1); % Setting carrying capacity to 1. Differecnes in A. 
r = ones(N, 1);
d = 0.5;
prob_edge = (d*N)/(N-1);
factor = 0.1; % Interaction reduction factor to ensure stabilization. 
A = factor.*(2*(rand(N, N) - 0.5)) .*...
    (rand(N, N) < prob_edge); % "Deciphering" paper
A(1:N+1:end) = 1;

% "Universality" interactions. Not using it now. 
% sigma = 1;
% sigma_max = 0.1;
% var = (sigma*sigma_max)^2;
% A = randn(N, 1);
% A(1:N+1:end) = 1;

% Inverting the graph
A = A';

% Creating ER keystone - with distribution
mu = 0.5;
sigma = 1;
keystone_factor = lognrnd(mu, sigma, N, 1);
for i = 1:N
    A(i, 1:end) = A(i, 1:end).*keystone_factor(i);
    % A(1:end, ind_key) = 0; % THIS RUINES CORRELATION BUT NOT KEYSTONE
    r(i) = 1;
    A(i, i) = 1;
end

%% Creating samples
overlap_factor = 0.5;
lotka_volterra = @(t, x, Ai, ri, ki) x.*(ri - Ai'*x./ki);

% ER
X0 = (rand(N, num_samples).*(rand(N, num_samples)<=overlap_factor));
data = zeros(N, num_samples);
i = 1;
while i <= num_samples
    i
    [T, X] = ode45(@(t, x) lotka_volterra(t, x, A, r, k), [t0, tf], X0(:, i));
    if sum(X(end, :)>1000)==0
       data(:, i) = X(end, :);
       i = i + 1;
    end
end
data(data < threshold) = 0; % Setting threshold

%% Normalizing samples
data_norm = data./sum(data, 2);
data = data_norm;
datad = double(data>0);

%% Calculating keystonennes
% For each species in each sample, reverse the presence/absance of the
% species and calculate the distance before and after. Then, calculate
% the average difference and call that keystonennes.

% ER:
keystonennes = nan(N, num_samples);
for i = 1:num_samples
    i
    for j = 1:N
        % Simulating after removal/addition
        original_sample = data(:, i);
        X0 = original_sample;
        X0(j) = not(X0(j));
        [T, X] = ode45(@(t, x) lotka_volterra(t, x, A, r, k), [t0, tf], X0);
        X = X(end, :); X(X<threshold) = 0; 
        
        % Normalizing before/after
        original_sample = original_sample./sum(original_sample);
        X = X./sum(X);
        
        % Deleting species
        original_sample(j) = [];
        X(j) = [];
        
        % Renormalizing samples
        original_sample = original_sample./sum(original_sample);
        X = X./sum(X);
        
        % Calculating distance
        % rJSD:
%         keystonennes(j, i) = pdist2(original_sample', X, @rJSD_pdist2);

        % Spearman:
%         keystonennes(j, i) = pdist2(original_sample', X, 'spearman');

        % Manhattan (Equivilant to BC of normalized data)
%           keystonennes(j, i) = sum(abs(original_sample'-X));
          
        % BC:
        keystonennes(j, i) = pdist2(original_sample', X, @BC_pdist2);
    end
end
keystonennes = mean(keystonennes, 2);

%% Calculating Q/D1/D2 for all speices
D1 = nan(N, 1);
D2 = nan(N, 1);
for i = 1:N
    i
    % Q
    data_i = data; data_i(i, :) = [];
    
    % Renormalizing
    data_i = data_i./sum(data_i);
    
    distances_i = pdist(data_i', @BC_pdist2);
    [cdf_dist, dist] = ecdf(distances_i);
    dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last')); % Finding distance threshold
    B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
    B_i(1:M+1:end) = 0;
    B_i_graph = graph(B_i);
    s_i = double(datad(i, :)); s_i(s_i==0) = -1; s_i = s_i';
    Q(i) = modularity_guy(B_i, s_i);
    
    % D1/D2
    ind_pres = data(i, :) ~= 0;
    data_pres = data(:, ind_pres);
    data_abs = data(:, not(ind_pres));
    data_pres(i, :) = [];
    data_abs(i, :) = [];
    
    % Renormalizing
    data_pres = data_pres./sum(data_pres);
    data_abs = data_abs./sum(data_abs);
    
    D1(i) = sum(sum(pdist2(data_pres', data_abs', @BC_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));
    D2(i) = pdist2(mean(data_pres'), mean(data_abs'), @BC_pdist2);    
end

%% Calculating degree for all species
% Creating co-occurence network
D = 1 - pdist(data, 'spearman');
D(isnan(D)) = 0; % Should I do that?

% Calculating p-values with reshuffling
D_shuffle = zeros(num_reshuffle_degree, N*(N-1)/2);
for m = 1:num_reshuffle_degree
    m
    data_shuffle = shuffle_data(data);
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

%% Plot
figure_7 = figure;
set(gcf, 'Position', [20 100 800*2 300*2], 'Units', 'centimeters');
markersize = 20;
% textxvalue = 0.05;
% textyvalue = 0.85;
% textfont = 15;

subplot(2, 3, 3); box off;
plot(keystonennes, Q, '.k', 'MarkerSize', markersize);
xlabel('$I$', 'Interpreter', 'latex');
ylabel('Q');
text(0, 1.1, 'c', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
[r, p] = corr(keystonennes, Q');
text(textxvalue, textyvalue, ['r = ', ' ', num2str(r), newline,...
    'p = ', ' ', num2str(p)],...
    'Units', 'normalized', 'EdgeColor', 'none',...
    'FontSize', textfont, 'BackgroundColor', [1 1 1]);
set(gca, 'FontSize', fontsize);
axis square;
box off;

subplot(2, 3, 1); box off;
plot(keystonennes, D1, '.k', 'MarkerSize', markersize);
xlabel('$I$', 'Interpreter', 'latex');
ylabel('D_1');
text(0, 1.1, 'a', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
[r, p] = corr(keystonennes, D1);
text(textxvalue, textyvalue, ['r = ', ' ', num2str(r), newline,...
    'p = ', ' ', num2str(p)],...
    'Units', 'normalized', 'EdgeColor', 'none',...
    'FontSize', textfont, 'BackgroundColor',[1 1 1]);
set(gca, 'FontSize', fontsize);
set(gca, 'FontSize', fontsize);
axis square;
box off;
 
subplot(2, 3, 2);
plot(keystonennes, D2, '.k', 'MarkerSize', markersize);
xlabel('$I$', 'Interpreter', 'latex');
ylabel('D_2');
text(0, 1.1, 'b', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
[r, p] = corr(keystonennes, D2);
text(textxvalue, textyvalue, ['r = ', ' ', num2str(r), newline,...
    'p = ', ' ', num2str(p)],...
    'Units', 'normalized', 'EdgeColor', 'none',...
    'FontSize', textfont, 'BackgroundColor', [1 1 1]);
set(gca, 'FontSize', fontsize);
set(gca, 'FontSize', fontsize);
axis square;
box off;

% subplot(2, 4, 4);
% plot(keystonennes, degree, '.k', 'MarkerSize', markersize);
% xlabel('$I$', 'Interpreter', 'latex');
% ylabel('degree', 'Interpreter', 'latex');
% text(0, 1.1, '(d)', 'Units', 'normalized', 'FontSize', 18, 'Interpreter', 'latex');
% [r, p] = corr(keystonennes, degree);
% text(textxvalue, textyvalue, ['$r = $', ' ', num2str(r), newline,...
%     '$p_\mathrm{value} = $', ' ', num2str(p)],...
%     'Units', 'normalized', 'EdgeColor', 'k', 'Interpreter', 'latex',...
%     'FontSize', textfont, 'BackgroundColor', [1 1 1]);
% set(gca, 'FontSize', fontsize);
% axis square;









%% STRUCTURAL BASED @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%% Parameters
% clear;
N = 100; % Number of species
num_samples = 100; % Number of samples
M = num_samples;
t0 = 0; % Start time
tf = 100; % End time
threshold = 1e-5; % Threshold for dead species\
threshold_net = 0.5; % Percentile threshold for network construction of Q
num_reshuffling = 20; % Number of reshuffling for pvalue calculations D/Q
hist_color = [0.4, 0.4, 0.4];
num_bins = 100;
fontsize = 15;
bcolor = [0  113.9850  188.9550]/255;
rcolor = [216.7500   82.8750   24.9900]/255;
kcolor = [0 0 0];
markersize = 10;
num_reshuffle_degree = 1000;

%% Creating networks
% Creating BA keystone
% --- Barabasi-Albert ---
k = ones(N, 1);
r = ones(N, 1); % Change to ones for consistency 
% My new function:
m_0 = 3;
m = 2; % 2 edges per one (Mean degree =~= 2*m);
d = 0.1;
factor_BA = 0.5;
A = Barabasi_Albert_2_DIRECTED(N, m_0, m, d);
A = full(A.adjacency);
A = factor_BA.*A .* (2*(rand(N, N) - 0.5));
% A = A';
A(1:N+1:end) = 1;

%% Creating samples
overlap_factor = 0.5;
lotka_volterra = @(t, x, Ai, ri, ki) x.*(ri - Ai'*x./ki);

% ER
X0 = (rand(N, num_samples).*(rand(N, num_samples)<=overlap_factor));
data = zeros(N, num_samples);
i = 1;
while i <= num_samples
    i
    [T, X] = ode45(@(t, x) lotka_volterra(t, x, A, r, k), [t0, tf], X0(:, i));
    if sum(X(end, :)>1000)==0
       data(:, i) = X(end, :);
       i = i + 1;
    end
end
data(data < threshold) = 0; % Setting threshold

%% Normalizing samples
data_norm = data./sum(data, 2);
data = data_norm;
datad = double(data>0);

%% Calculating keystonennes
% For each species in each sample, reverse the presence/absance of the
% species and calculate the distance before and after. Then, calculate
% the average difference and call that keystonennes.

% ER:
keystonennes = nan(N, num_samples);
for i = 1:num_samples
    i
    for j = 1:N
        % Simulating after removal/addition
        original_sample = data(:, i);
        X0 = original_sample;
        X0(j) = not(X0(j));
        [T, X] = ode45(@(t, x) lotka_volterra(t, x, A, r, k), [t0, tf], X0);
        X = X(end, :); X(X<threshold) = 0; 
        
        % Normalizing before/after
        original_sample = original_sample./sum(original_sample);
        X = X./sum(X);
        
        % Deleting species
        original_sample(j) = [];
        X(j) = [];
        
        % Renormalizing
        original_sample = original_sample./sum(original_sample);
        X = X./sum(X);
        
        % Calculating distance
        % rJSD:
%         keystonennes(j, i) = pdist2(original_sample', X, @rJSD_pdist2);

        % Spearman:
%         keystonennes(j, i) = pdist2(original_sample', X, 'spearman');

        % Manhattan (Equivilant to BC of normalized data)
%           keystonennes(j, i) = sum(abs(original_sample'-X));
          
        % BC:
        keystonennes(j, i) = pdist2(original_sample', X, @BC_pdist2);
    end
end
keystonennes = mean(keystonennes, 2);

%% Calculating Q/D1/D2 for all speices
D1 = nan(N, 1);
D2 = nan(N, 1);
for i = 1:N
    i
    % Q
    data_i = data; data_i(i, :) = [];
    
    % Renormalizing
    data_i = data_i./sum(data_i);
    
    distances_i = pdist(data_i', @BC_pdist2);
    [cdf_dist, dist] = ecdf(distances_i);
    dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last')); % Finding distance threshold
    B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
    B_i(1:M+1:end) = 0;
    B_i_graph = graph(B_i);
    s_i = double(datad(i, :)); s_i(s_i==0) = -1; s_i = s_i';
    Q(i) = modularity_guy(B_i, s_i);
    
    % D1/D2
    ind_pres = data(i, :) ~= 0;
    data_pres = data(:, ind_pres);
    data_abs = data(:, not(ind_pres));
    data_pres(i, :) = [];
    data_abs(i, :) = [];
    
    % Renormalizing
    data_pres = data_pres./sum(data_pres);
    data_abs = data_abs./sum(data_abs);
    
    D1(i) = sum(sum(pdist2(data_pres', data_abs', @BC_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));
    D2(i) = pdist2(mean(data_pres'), mean(data_abs'), @BC_pdist2);    
end

%% Calculating degree for all species
% Creating co-occurence network
D = 1 - pdist(data, 'spearman');
D(isnan(D)) = 0; % Should I do that?

% Calculating p-values with reshuffling
D_shuffle = zeros(num_reshuffle_degree, N*(N-1)/2);
for m = 1:num_reshuffle_degree
    m
    data_shuffle = shuffle_data(data);
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

%% Plot

% figure;
% set(gcf, 'Position', [20 100 900*2 300*2], 'Units', 'centimeters');
markersize = 20;
% textxvalue = 0.5;
% textyvalue = 0.1;

subplot(2, 3, 6); 
plot(keystonennes, Q, '.k', 'MarkerSize', markersize);
xlabel('$I$', 'Interpreter', 'latex');
ylabel('Q');
text(0, 1.1, 'f', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
[r, p] = corr(keystonennes, Q');
text(textxvalue, textyvalue, ['r = ', ' ', num2str(r), newline,...
    'p = ', ' ', num2str(p)],...
    'Units', 'normalized', 'EdgeColor', 'none',...
    'FontSize', textfont, 'BackgroundColor', [1 1 1]);
set(gca, 'FontSize', fontsize);
axis square;
box off;

subplot(2, 3, 4); box off;
plot(keystonennes, D1, '.k', 'MarkerSize', markersize);
xlabel('$I$', 'Interpreter', 'latex');
ylabel('D_1');
text(0, 1.1, 'd', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
[r, p] = corr(keystonennes, D1);
text(textxvalue, textyvalue, ['r = ', ' ', num2str(r), newline,...
    'p = ', ' ', num2str(p)],...
    'Units', 'normalized', 'EdgeColor', 'none',...
    'FontSize', textfont, 'BackgroundColor', [1 1 1]);
set(gca, 'FontSize', fontsize);
set(gca, 'FontSize', fontsize);
axis square;
box off;

subplot(2, 3, 5);
plot(keystonennes, D2, '.k', 'MarkerSize', markersize);
xlabel('$I$', 'Interpreter', 'latex');
ylabel('D_2');
text(0, 1.1, 'e', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
[r, p] = corr(keystonennes, D2);
text(textxvalue, textyvalue, ['r = ', ' ', num2str(r), newline,...
    'p = ', ' ', num2str(p)],...
    'Units', 'normalized', 'EdgeColor', 'none',...
    'FontSize', textfont, 'BackgroundColor', [1 1 1]);
set(gca, 'FontSize', fontsize);
set(gca, 'FontSize', fontsize);
axis square;
box off;

% subplot(2, 4, 8);
% plot(keystonennes, degree, '.k', 'MarkerSize', markersize);
% xlabel('$I$', 'Interpreter', 'latex');
% ylabel('degree', 'Interpreter', 'latex');
% text(0, 1.1, '(h)', 'Units', 'normalized', 'FontSize', 18, 'Interpreter', 'latex');
% [r, p] = corr(keystonennes, degree);
% text(textxvalue, textyvalue, ['$r = $', ' ', num2str(r), newline,...
%     '$p_\mathrm{value} = $', ' ', num2str(p)],...
%     'Units', 'normalized', 'EdgeColor', 'k', 'Interpreter', 'latex',...
%     'FontSize', textfont, 'BackgroundColor', [1 1 1]);
% set(gca, 'FontSize', fontsize);
% axis square;





