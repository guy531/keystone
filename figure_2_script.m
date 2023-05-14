%% Figure 2
% Demonstration of a single strength-based and structre-based keystone in
% GLV dynamics

%% Parameters
N = 100; % Number of species
num_samples = 100; % Number of samples
t0 = 0; % Start time
tf = 100; % End time
threshold = 1e-5; % Threshold for dead species
fontsize = 15;
hist_color = [0.4, 0.4, 0.4];

% Distance function
% Kullback-Leibler Divergence
KLD=@(x,y) sum(x.*log(x./y));

% Jason-Shanon Divergence
rJSD=@(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2));


%% Creating networks
% Erdos-Renyi
k_ER = ones(N, 1); % Setting carrying capacity to 1. Differecnes in A. 
r_ER = ones(N, 1); % Change to ones for consistency 
d = 0.1;
prob_edge = (d*N)/(N-1);
factor_ER = 0.1; % Interaction reduction factor to ensure stabilization. 
A_ER = factor_ER.*(2*(rand(N, N) - 0.5)) .*...
    (rand(N, N) < prob_edge); % "Deciphering" paper
A_ER(1:N+1:end) = 1;

% "Universality" interactions. Not using it now. 
% sigma = 1;
% sigma_max = 0.1;
% var_ER = (sigma*sigma_max)^2;
% A_ER = randn(N, 1);
% A_ER(1:N+1:end) = 1;

% Inverting the graph
A_ER = A_ER';

% Creating ER keystone
num_keystone = 1;
keystone_factor = 10;
ind_key = randperm(N); ind_key = ind_key(1:num_keystone);
A_ER(ind_key, 1:end) = A_ER(ind_key, 1:end).*keystone_factor;
r_ER(ind_key) = 1;
% k(ind_key) = 200;
% A_ER(1:end, ind_key) = 0; % THIS RUINES CORRELATION BUT NOT KEYSTONE
A_ER(ind_key, ind_key) = 1;
% A(ind_key, 1:end) = (rand(num_keystone, N) - 0.5).*keystone_factor;




% Creating BA keystone
% --- Barabasi-Albert ---
k_BA = ones(N, 1);
r_BA = ones(N, 1); % Change to ones for consistency 
% My new function:
m_0 = 3;
m = 2; % 2 edges per one (Mean degree =~= 2*m);
d = 0.1;
factor_BA = 0.1;
A_BA = Barabasi_Albert_2_DIRECTED(N, m_0, m, d);
A_BA = full(A_BA.adjacency);
A_BA = factor_BA.*A_BA .* (2*(rand(N, N) - 0.5));
A_BA(1:N+1:end) = 1;


%% Creating samples
overlap_factor = 0.8;
lotka_volterra = @(t, x, Ai, ri, ki) x.*(ri - Ai'*x./ki);

% ER
X0 = (rand(N, num_samples).*(rand(N, num_samples)<=overlap_factor));
data_ER = zeros(N, num_samples);
i = 1;
while i <= num_samples
    i
    [T, X] = ode45(@(t, x) lotka_volterra(t, x, A_ER, r_ER, k_ER), [t0, tf], X0(:, i));
    if sum(X(end, :)>1000)==0
       data_ER(:, i) = X(end, :);
       i = i + 1;
    end
end
data_ER(data_ER < threshold) = 0; % Setting threshold

% BA
X0 = (rand(N, num_samples).*(rand(N, num_samples)<=overlap_factor));
data_BA = zeros(N, num_samples);
i = 1;
while i <= num_samples
    i
    [T, X] = ode45(@(t, x) lotka_volterra(t, x, A_BA, r_BA, k_BA), [t0, tf], X0(:, i));
    if sum(X(end, :)>1000)==0
       data_BA(:, i) = X(end, :);
       i = i + 1;
    end
end
data_BA(data_BA < threshold) = 0; % Setting threshold


%% Normalizing samples
data_ER_norm = data_ER./sum(data_ER, 2);
data_BA_norm = data_BA./sum(data_BA, 2);

%% Calculating keystonennes
% For each species in each sample, reverse the presence/absance of the
% species and calculate the distance before and after. Then, calculate
% the average difference and call that keystonennes.

% ER:
keystonennes_ER = nan(N, num_samples);
for i = 1:num_samples
    i
    for j = 1:N
        % Simulating after removal/addition
        original_sample = data_ER(:, i);
        X0 = original_sample;
        X0(j) = not(X0(j));
        [T, X] = ode45(@(t, x) lotka_volterra(t, x, A_ER, r_ER, k_ER), [t0, tf], X0);
        X = X(end, :); X(X<threshold) = 0; 
        
        % Normalizing before/after
        original_sample = original_sample./sum(original_sample);
        X = X./sum(X);
        
        % Deleting species
        original_sample(j) = [];
        X(j) = [];
        
        % Renormalizing after removal
        original_sample = original_sample./sum(original_sample);
        X = X./sum(X);
        
        % Calculating distance
        % rJSD:
%         keystonennes_ER(j, i) = pdist2(original_sample', X, @rJSD_pdist2);

        % Spearman:
%         keystonennes_ER(j, i) = pdist2(original_sample', X, 'spearman');

        % Manhattan (Equivilant to BC of normalized data)
%           keystonennes_ER(j, i) = sum(abs(original_sample'-X));
          
        % BC:
        keystonennes_ER(j, i) = pdist2(original_sample', X, @BC_pdist2);
    end
end
keystonennes_ER = mean(keystonennes_ER, 2);

% BA:
keystonennes_BA = nan(N, num_samples);
for i = 1:num_samples
    i
    for j = 1:N
        % Simulating after removal/addition
        original_sample = data_BA(:, i);
        X0 = original_sample;
        X0(j) = not(X0(j));
        [T, X] = ode45(@(t, x) lotka_volterra(t, x, A_BA, r_BA, k_BA), [t0, tf], X0);
        X = X(end, :); X(X<threshold) = 0;
        
        % Normalizing before/after
        original_sample = original_sample./sum(original_sample);
        X = X./sum(X);
        
        % Removing species
        X(j) = [];
        original_sample(j) = [];
        
        % Renormalizing after removal
        original_sample = original_sample./sum(original_sample);
        X = X./sum(X);
        
        % Calculating distance
        % rJSD:
%         keystonennes_BA(j, i) = pdist2(original_sample', X, @rJSD_pdist2); 
        
        % Spearmn:
%         keystonennes_BA(j, i) = pdist2(original_sample', X, 'spearman');

        % Manhattan (Equivilant to BC of normalized data)
        % keystonennes_BA(j, i) = sum(abs(original_sample'-X));   
        
        % BC:
        keystonennes_BA(j, i) = pdist2(original_sample', X, @BC_pdist2);
    end
end
keystonennes_BA = mean(keystonennes_BA, 2);

%% Plotting
% https://www.mathworks.com/matlabcentral/answers/607141-highlighting-edges-of-a-graph
rcolor = [216.7500   82.8750   24.9900]/255;
bcolor = 0.8.*[0  113.9850  188.9550]/255;
ncolor = [1 1 1].*0.5; 

figure_1 = figure; 
hold on;
set(gcf, 'Position', [20 100 900 900], 'Units', 'centimeters');
subplot(2, 2, 1); hold on; axis square;
temp = A_ER; temp(1:N+1:end) = 0;
temp = digraph(temp);
G1 = plot(temp, 'NodeLabel', [], 'EdgeColor', ncolor, 'NodeColor', ncolor);
G2 = plot(temp, 'EdgeAlpha', 1, 'NodeLabel', [], 'Marker', 'none', ...
    'LineWidth', 1.5, 'EdgeColor', bcolor, 'ArrowSize', 10);
% G2 = plot(temp, 'NodeLabel', [], 'EdgeCData', rand(1, 4), 'NodeColor', 'none');
n_key_1 = predecessors(temp, ind_key);
n_key_2 = successors(temp, ind_key);
% Predecessors - Don't need it

% Find all edges not connected to the one I want and delete them from the
% second graph
t = temp.Edges(:, 1:2);
t = table2array(t);
t(t(:, 1) == ind_key, :) = [];
t(:, 3) = [];

% Successors
for i = 1:size(t, 1)
    i
    highlight(G2, [t(i, 1), t(i, 2)], 'LineWidth', 1.5, 'LineStyle', 'none');
end
highlight(G1, ind_key, 'MarkerSize', 10, 'NodeColor', bcolor);

set(gca, 'XTick', []);
set(gca, 'YTick', []);
box off;
axis off;
set(gca, 'FontSize', fontsize);
title('Interaction network', 'FontSize', 20, 'FontWeight', 'normal');
text(-0.1, 1.1, 'a', 'Units', 'normalized', 'FontSize', 25, 'FontWeight', 'bold');

% Legened:
l = line(nan,nan, 'Color', bcolor, 'LineWidth', 2);
legend(l, 'Boosted interactions');

%
subplot(2, 2, 2); axis square;
num_bins = 100;
histogram(keystonennes_ER,num_bins, 'EdgeColor', 'none', 'FaceColor', hist_color);
xlabel('Presence-impact, $I^i$', 'Interpreter', 'latex');
ylabel('Number of species');

set(gca, 'FontSize', fontsize);
text(-0.25, 1.1, 'b', 'Units', 'normalized', 'FontSize', 25, 'FontWeight', 'bold');
box off;

%
subplot(2, 2, 3); hold on; axis square;
temp = A_BA; temp(1:N+1:end) = 0;
temp = digraph(temp);
G3 = plot(temp, 'NodeLabel', [], 'EdgeColor', ncolor, 'NodeColor', ncolor);
G4 = plot(temp, 'EdgeAlpha', 1, 'NodeLabel', [], 'Marker', 'none', ...
    'LineWidth', 1.5, 'EdgeColor', [0 0 0], 'ArrowSize', 10);
% Find highest keystonenness and highlight
[~, ind_key_BA] = max(keystonennes_BA);
n_key_1 = predecessors(temp, ind_key_BA);
n_key_2 = successors(temp, ind_key_BA);

t = temp.Edges(:, 1:2);
t = table2array(t);
t(t(:, 1) == ind_key_BA, :) = [];
t(:, 3) = [];

for i = 1:size(t, 1)
    i
    highlight(G4, [t(i, 1), t(i, 2)], 'LineWidth', 1.5, 'LineStyle', 'none');
end

highlight(G3, ind_key_BA, 'MarkerSize', 10, 'NodeColor', ncolor);

set(gca, 'XTick', []);
set(gca, 'YTick', []);


set(gca, 'FontSize', fontsize);
title('Interaction network', 'FontSize', 20, 'FontWeight', 'normal');
text(-0.1, 1.1, 'c', 'Units', 'normalized', 'FontSize', 25, 'FontWeight', 'bold');
box off;
axis off;
% Legened:
l = line(nan,nan, 'Color', [0 0 0], 'LineWidth', 2);
legend(l, 'Hub interactions');

subplot(2, 2, 4); axis square;
num_bins = 100;
histogram(keystonennes_BA, num_bins, 'EdgeColor', 'none', 'FaceColor', hist_color);
xlabel('Presence-impact, $I^i$', 'Interpreter', 'latex');
ylabel('Number of species');

set(gca, 'FontSize', fontsize);
text(-0.25, 1.1, 'd', 'Units', 'normalized', 'FontSize', 25, 'FontWeight', 'bold');
box off;

%% Calculating pvalues
% Interaction
z_vector = zscore(keystonennes_ER);
pvalue_1 = pvaluefromz(z_vector(ind_key));

% Structure
z_vector = zscore(keystonennes_BA);
pvalue_2 = pvaluefromz(z_vector(ind_key_BA));
