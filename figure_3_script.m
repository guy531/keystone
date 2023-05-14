%% Calculating distance/modularity for keystone/random with PCoA/network
% Only ER with Bray-Curtis?

%% Parameters
N = 100; % Number of species
num_samples = 100; % Number of samples
M = num_samples;
t0 = 0; % Start time
tf = 100; % End time
threshold = 1e-5; % Threshold for dead species\
hist_color = [0.4, 0.4, 0.4];
num_bins = 100;
fontsize = 15;

% Distance functions
% Kullback-Leibler Divergence
KLD = @(x,y) sum(x.*log(x./y));
% Jason-Shanon Divergence
rJSD = @(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2));

bcolor = [0  113.9850  188.9550]/255;
rcolor = [216.7500   82.8750   24.9900]/255;
kcolor = [0 0 0];
node_color_parm = 0.8;


%% Creating networks
% Erdos-Renyi
k_ER = ones(N, 1); % Setting carrying capacity to 1. Differecnes in A. 
r_ER = ones(N, 1);
d = 0.5;
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

%% Creating samples
overlap_factor = 0.5;
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


%% Normalizing samples
data_ER_norm = data_ER./sum(data_ER, 2);
data_ER = data_ER_norm;

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
        X(j) = [];
        original_sample(j) = [];
        
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

%% Choosing indecies
temp = 1:N; temp(ind_key) = [];
ind_rand = temp(randi([1, N-1], 1));

%% Calculating distance/PCoA
% Using two types of distances
figure_3 = figure; 
set(gcf, 'Position', [20 100 1000 1000], 'Units', 'centimeters');
letterposx = -0.25;
letterposy = 1.2;
% Using Bray-Curtis

% AVERAGE OF DISTANCES ----------------------------------------------------
% Key ---------------------------------------------------------------------
subplot(3, 3, 2); box off;
ind_pres = data_ER(ind_key, :) ~= 0;
data_pres = data_ER(:, ind_pres); 
data_abs = data_ER(:, not(ind_pres)); 
% Deleting species
data_pres(ind_key, :) = [];
data_abs(ind_key, :) = [];
% Calculating distances between all samples
distances_key = pdist2(data_ER_norm', data_ER_norm', @BC_pdist2);
distances_key(1:num_samples+1:end) = 0;
y = mdscale(distances_key, 2);
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

% legend('Present', 'Absent', 'FontSize', 13);

text(0.5, 1.1, 'Simulated keystone', 'Units', 'normalized', 'FontSize', 18, 'HorizontalAlignment', 'center');
axis square;
D_key = sum(sum(pdist2(data_pres', data_abs', @BC_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));



% Rand --------------------------------------------------------------------
subplot(3, 3, 3); box off;
ind_pres = data_ER(ind_rand, :) ~= 0;
data_pres = data_ER(:, ind_pres); 
data_abs = data_ER(:, not(ind_pres)); 
% Deleting species
data_pres(ind_rand, :) = [];
data_abs(ind_rand, :) = [];
% Calculating distances between all samples
distances_key = pdist2(data_ER_norm', data_ER_norm', @BC_pdist2);
distances_key(1:num_samples+1:end) = 0;
y = mdscale(distances_key, 2);
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
% legend('Present', 'Absent', 'FontSize', 10);
text(0.5, 1.1, 'Random species', 'Units', 'normalized', 'FontSize', 18, 'HorizontalAlignment', 'center');
axis square;
D_rand = sum(sum(pdist2(data_pres', data_abs', @BC_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));



% Histogram ---------------------------------------------------------------
subplot(3, 3, 1);
% Calculating Di for all species
Di = nan(N, 1);
for i = 1:N
    i;
    ind_pres = data_ER(i, :) ~= 0;
    data_pres = data_ER(:, ind_pres); 
    data_abs = data_ER(:, not(ind_pres));
    data_pres(i, :) = [];
    data_abs(i, :) = [];
    
    % Renormalizing the samples
    data_pres = data_pres./sum(data_pres);
    data_abs = data_abs./sum(data_abs);
    
    Di(i) = sum(sum(pdist2(data_pres', data_abs', @BC_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));
end
h1 = histogram(Di, num_bins, 'EdgeColor', 'none', 'FaceColor', hist_color);
hold on;
ylimits = get(gca, 'Ylim'); 
ymin = ylimits(1); ymax = ylimits(2);
xlimits = get(gca, 'Xlim');
xmin = xlimits(1); xmax = xlimits(2);
% plot([D_key, D_key], [0, ymax/2]);
Annotate(gca, 'arrow', [xmin + (xmax-xmin)*(2/3), D_key], [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', bcolor);
Annotate(gca, 'arrow', [xmin + (xmax-xmin)/3, D_rand],  [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', rcolor);
% annotation('arrow', [D_key, D_key], [0, ymax/2]);
% xlabel('$D^i_{1}$', 'Interpreter', 'latex');
xlabel('D_1');
ylabel('Number of species');
set(gca, 'FontSize', fontsize);
text(letterposx, letterposy, 'a', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
% text(-0.5, 0.5, 'D_1', 'Units', 'normalized', 'FontSize', 18, 'HorizontalAlignment', 'center');
box off;
axis square;







% DISTANCES BETWEEN AVERAGES ----------------------------------------------
% Key ---------------------------------------------------------------------
subplot(3, 3, 5); box off;
ind_pres = data_ER(ind_key, :) ~= 0;
data_pres = data_ER(:, ind_pres); 
data_abs = data_ER(:, not(ind_pres)); 
% Deleting species
data_pres(ind_key, :) = [];
data_abs(ind_key, :) = [];
data_pres_avg = mean(data_pres'); % The average of the pres samples
data_abs_avg = mean(data_abs'); % The average of the abs samples
% Calculating distances between all samples
distances_key = pdist2(data_ER_norm', data_ER_norm', @BC_pdist2);
distances_key(1:num_samples+1:end) = 0;
distances_avg_key = pdist2(data_pres_avg, data_abs_avg, @BC_pdist2);
y = mdscale(distances_key, 2);
hold on;
% plot(y(ind_pres, 1), y(ind_pres, 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*bcolor);
% plot(y(not(ind_pres), 1), y(not(ind_pres), 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*rcolor);
sz = 30;
scatter(y(ind_pres, 1), y(ind_pres, 2), sz, bcolor, 'filled');
scatter(y(not(ind_pres), 1), y(not(ind_pres), 2), sz, bcolor);
% Plotting average
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
D_key_2 = pdist2(data_pres_avg, data_abs_avg, @BC_pdist2);



% Rand --------------------------------------------------------------------
subplot(3, 3, 6); box off;
ind_pres = data_ER(ind_rand, :) ~= 0;
data_pres = data_ER(:, ind_pres); 
data_abs = data_ER(:, not(ind_pres)); 
% Deleting species
data_pres(ind_rand, :) = [];
data_abs(ind_rand, :) = [];
data_pres_avg = mean(data_pres'); % The average of the pres samples
data_abs_avg = mean(data_abs'); % The average of the abs samples
% Calculating distances between all samples
distances_key = pdist2(data_ER_norm', data_ER_norm', @BC_pdist2);
distances_key(1:num_samples+1:end) = 0;
y = mdscale(distances_key, 2);
hold on;
% plot(y(ind_pres, 1), y(ind_pres, 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*bcolor);
% plot(y(not(ind_pres), 1), y(not(ind_pres), 2), '.', 'MarkerSize', 20, 'Color', node_color_parm.*rcolor);
sz = 30;
scatter(y(ind_pres, 1), y(ind_pres, 2), sz, rcolor, 'filled');
scatter(y(not(ind_pres), 1), y(not(ind_pres), 2), sz, rcolor);
% Plotting average
y_pres_avg = [mean(y(ind_pres, 1)), mean(y(ind_pres, 2))];
y_abs_avg = [mean(y(not(ind_pres), 1)), mean(y(not(ind_pres), 2))];
plot(y_pres_avg(1), y_pres_avg(2), 'X', 'MarkerSize', 15, 'Color', 'k', 'LineWidth', 2);
plot(y_abs_avg(1), y_abs_avg(2), 'X', 'MarkerSize', 15, 'Color', 'k', 'LineWidth', 2)
plot([y_pres_avg(1), y_abs_avg(1)], [y_pres_avg(2), y_abs_avg(2)], 'k', 'LineWidth', 1.5)
xlabel('PC1'); ylabel('PC2');
set(gca, 'FontSize', fontsize);
set(gca, 'XTick', []); set(gca, 'YTick', []);
text(letterposx, letterposy, 'f', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
% legend('Present', 'Absent', 'FontSize', 10);
axis square;
D_rand_2 = pdist2(data_pres_avg, data_abs_avg, @BC_pdist2);


% Histogram ---------------------------------------------------------------
subplot(3, 3, 4);
% Calculating Di for all species
Di_2 = nan(N, 1);
for i = 1:N
    i;
    ind_pres = data_ER(i, :) ~= 0;
    data_pres = data_ER(:, ind_pres); 
    data_abs = data_ER(:, not(ind_pres));
    data_pres(i, :) = [];
    data_abs(i, :) = [];
    
    % Renormalizing the samples
    data_pres = data_pres./sum(data_pres);
    data_abs = data_abs./sum(data_abs);
    
    Di_2(i) = pdist2(mean(data_pres'), mean(data_abs'), @BC_pdist2);
end
h1 = histogram(Di_2, num_bins, 'EdgeColor', 'none', 'FaceColor', hist_color);
hold on;
ylimits = get(gca, 'Ylim'); 
ymin = ylimits(1); ymax = ylimits(2);
xlimits = get(gca, 'Xlim');
xmin = xlimits(1); xmax = xlimits(2);
% plot([D_key, D_key], [0, ymax/2]);
Annotate(gca, 'arrow', [xmin + (xmax-xmin)*(2/3), D_key_2], [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', bcolor);
Annotate(gca, 'arrow', [xmin + (xmax-xmin)/3, D_rand_2],  [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', rcolor);
% annotation('arrow', [D_key, D_key], [0, ymax/2]);
% xlabel('$D^i_{2}$', 'Interpreter', 'latex');
xlabel('D_2');
ylabel('Number of species');
set(gca, 'FontSize', fontsize);
text(letterposx, letterposy, 'd', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
% text(-0.1, 0.5, 'D_2^i', 'Units', 'normalized', 'FontSize', 18);
% text(-0.5, 0.5, 'D_2', 'Units', 'normalized', 'FontSize', 18, 'HorizontalAlignment', 'center');
box off;
axis square;








% Calculating modularity/network
linewidth = 1.5;
markersize = 3;
threshold_net = 0.1; % Precentile threshold for network construction

% Using Bray-Curtis

% Key ---------------------------------------------------------------------
subplot(3, 3, 8); hold on;
% Creating network of samples without species i
data_i = data_ER; data_i(ind_key, :) = [];

% Renormalzing the samples
data_i = data_i./sum(data_i);

distances_i = pdist(data_i', @BC_pdist2);
[cdf_dist,dist] = ecdf(distances_i);
dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last')); % Finding distance threshold
% Creating co-occurence network
B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
B_i(1:M+1:end) = 0;
B_i_graph = graph(B_i);

% Calculating modularity
s_i = double(data_ER(ind_key, :)==0); s_i(s_i==0) = -1; s_i = s_i';
Q_key = modularity_guy(B_i, s_i);

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

% % Plotting network
% h1 = plot(B_i_graph, 'NodeLabel', []);
% axis square;
% 
% % Highlighting nodes
% temp = find(s_i==1);
% highlight(h1, temp, 'NodeColor', node_color_parm.*bcolor, 'MarkerSize', markersize);
% temp = find(s_i==-1);
% highlight(h1, temp, 'NodeColor', node_color_parm.*rcolor, 'MarkerSize', markersize);
% 
% % Highlighting edges
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
% title('Sample-to-sample similarity network, B^i');
% text(0.5, 1.3, 'Sample-to-sample similarity network, B^i', 'Units', 'normalized', 'FontSize', 18, 'HorizontalAlignment', 'center');
text(0.5, -0.1 ,'Sample-to-sample network, B^i', 'Units', 'normalized', 'FontSize', 13, 'HorizontalAlignment', 'center');
box off;
axis off;
set(gca, 'FontSize', fontsize);

% Rand --------------------------------------------------------------------
subplot(3, 3, 9); hold on;
% Creating network of samples without species i
data_i = data_ER; data_i(ind_rand, :) = [];

% Renormalizing the samples
data_i = data_i./sum(data_i);

distances_i = pdist(data_i', @BC_pdist2);
[cdf_dist,dist] = ecdf(distances_i);
dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last')); % Finding distance threshold
% Creating co-occurence network
B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
B_i(1:M+1:end) = 0;
B_i_graph = graph(B_i);

% Calculating modularity
s_i = double(data_ER(ind_rand, :)==0); s_i(s_i==0) = -1;  s_i = s_i';
Q_rand = modularity_guy(B_i, s_i);

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
scatter(xvalues(temp), yvalues(temp), sz, rcolor, 'MarkerFaceColor', [1 1 1]);
temp = find(s_i==1);
scatter(xvalues(temp), yvalues(temp), sz, rcolor, 'filled');

% % Plotting network
% h1 = plot(B_i_graph, 'NodeLabel', []);
% axis square;
% 
% % Highlighting nodes
% temp = find(s_i==1);
% highlight(h1, temp, 'NodeColor', node_color_parm.*bcolor, 'MarkerSize', markersize);
% temp = find(s_i==-1);
% highlight(h1, temp, 'NodeColor', node_color_parm.*rcolor, 'MarkerSize', markersize);
% 
% % Highlighting edges
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
text(letterposx, letterposy, 'i', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');

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

% Histogram ---------------------------------------------------------------
subplot(3, 3, 7);
% Calculating Qi for all species
Qi = nan(N, 1);
for i = 1:N
    i
    
    % Creating network of samples without species i
    data_i = data_ER; data_i(i, :) = [];
    
    % Renormalizing the samples
    data_i = data_i./sum(data_i);
    
    distances_i = pdist(data_i', @BC_pdist2);
    [cdf_dist,dist] = ecdf(distances_i);
    dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last')); % Finding distance threshold
    % Creating co-occurence network
    B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
    B_i(1:M+1:end) = 0;
    B_i_graph = graph(B_i);
    
    % Calculating modularity
    s_i = double(data_ER(i, :)==0); s_i(s_i==0) = -1; s_i = s_i';
    Qi(i) = modularity_guy(B_i, s_i);
end
h1 = histogram(Qi, num_bins, 'EdgeColor', 'none', 'FaceColor', hist_color);
hold on;
ylimits = get(gca, 'Ylim'); 
ymin = ylimits(1); ymax = ylimits(2);
xlimits = get(gca, 'Xlim');
xmin = xlimits(1); xmax = xlimits(2);
% plot([D_key, D_key], [0, ymax/2]);
Annotate(gca, 'arrow', [xmin + (xmax-xmin)*(2/3), Q_key], [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', bcolor);
Annotate(gca, 'arrow', [xmin + (xmax-xmin)/3, Q_rand],  [ymin + (ymax-ymin)/2, ymin + (ymax-ymin)/10], 'Color', rcolor);
% annotation('arrow', [D_key, D_key], [0, ymax/2]);
% xlabel('$Q^i$', 'Interpreter', 'latex');
xlabel('Q');
ylabel('Number of species');
set(gca, 'FontSize', fontsize);
text(letterposx, letterposy, 'g', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');
box off;
% text(-0.5, 0.5, 'Q', 'Units', 'normalized', 'FontSize', 18, 'HorizontalAlignment', 'center');
axis square;











