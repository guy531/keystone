% Loading the data new
load HMP_stool_1.mat;
load HMP_stool_2.mat;

hmp1 = HMP_stool_1;
hmp2 = HMP_stool_2;

%% Parameters
min_perc_samples = 0.1; % Min num of samples per OTU
N = 1000; % Total number of top OTU
threshold_min_pres = 0.2; % Minimum num of samples per OTU - FOR ANALYSIS
threshold_max_pres = 0.8; % Maximum num of samples per OTU - FOR ANALYSIS
threshold_net = 0.5; % Percentile threshold for network construction of Q
fontsize = 15;
num_reshuffling = 1000; % Number of reshuffling for pvalue calculations of top species
threshold_stat_test = 0.05; % Percentile of top OTU to compare

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

%% Shuffleing the data
% [temp1, temp2] = shuffle_data_for_two(hmp1.data, hmp2.data);
hmp1.data = shuffle_data_ravid_2(hmp1.data);
hmp2.data = shuffle_data_ravid_2(hmp2.data);

%% Removing constant rows
[data,PS] = removeconstantrows(hmp1.data);
hmp1.data = data;
hmp1.OTU(PS.remove,:) = [];

[data,PS] = removeconstantrows(hmp2.data);
hmp2.data = data;
hmp2.OTU(PS.remove,:) = [];

%% Intersecting the samples by the subject ID - SHOULD I DO THAT?
% [C, ia, ib] = intersect(hmp1.metadata.RSID, hmp2.metadata.RSID);

%% Normalizing the data
hmp1.data = hmp1.data./sum(hmp1.data);
hmp2.data = hmp2.data./sum(hmp2.data);

%% Ordering by mean abundance
[~, temp] = sort(mean(hmp1.data, 2), 'descend');
hmp1.data = hmp1.data(temp, :);
hmp1.OTU = hmp1.OTU(temp, :);

[~, temp] = sort(mean(hmp2.data, 2), 'descend');
hmp2.data = hmp2.data(temp, :);
hmp2.OTU = hmp2.OTU(temp, :);

%% Choosing top N OTUs
hmp1.data = hmp1.data(1:N, :);
hmp1.OTU = hmp1.OTU(1:N, :);

hmp2.data = hmp2.data(1:N, :);
hmp2.OTU = hmp2.OTU(1:N, :);

%% Intersecting and reordering the OTUs
[C, ia, ib] = intersect(hmp1.OTU(:, 1), hmp2.OTU(:, 1));

hmp1.data = hmp1.data(ia, :);
hmp1.OTU = hmp1.OTU(ia, :);

hmp2.data = hmp2.data(ib, :);
hmp2.OTU = hmp2.OTU(ib, :);

N = length(C);
hmp1.num_samples = size(hmp1.data, 2);
hmp2.num_samples = size(hmp2.data, 2);

%% Renormalizing
hmp1.data = hmp1.data./sum(hmp1.data);
hmp2.data = hmp2.data./sum(hmp2.data);

%% Coverting to presence-absence data
hmp1.datad = double(hmp1.data>0);
hmp2.datad = double(hmp2.data>0);

%% Calculating relative freq of all species
hmp1.rel_freq = sum(hmp1.datad, 2)./hmp1.num_samples;
hmp2.rel_freq = sum(hmp2.datad, 2)./hmp2.num_samples;

%% Calculating impact before/after
I_12 = nan(N, 1);
I_12_cell = cell(N, 1);
k_i = nan(N, 1);
for i = 1:N
    i
    I_12_i = [];
    k = 0;
    for j = 1:hmp2.num_samples
        % Testing if the sample had a "before" sample
        [temp1, temp2] = ismember(hmp2.metadata.RSID(j), hmp1.metadata.RSID);
        if temp1
            data1 = hmp1.data(:, temp2);
            datad1 = hmp1.datad(:, temp2);
            data2 = hmp2.data(:, j);
            datad2 = hmp2.datad(:, j);
            
            if datad1(i) ~= datad2(i)
                % Removing the species
                data1_i = data1; data1_i(i) = [];
                data2_i = data2; data2_i(i) = [];
                
                % Renormalizing
                data1_i = data1_i./sum(data1_i);
                data2_i = data2_i./sum(data2_i);
                
                % Calculating the distance before/after
                k = k + 1;
                I_12_i(k) = pdist2(data1_i', data2_i', @BC_pdist2);
                
            end
        end
    end
    k_i(i) = k;
    % Averaging over all the samples
    I_12(i) = mean(I_12_i);
    I_12_cell{i} = I_12_i;
end

%% Removing I_12 with k_i smaller then number
num_k = 10;
temp = k_i < num_k;
I_12(temp) = nan;

%% Calculating impact only before
% Do this only for species in a specific region?

% Indecies of rel_freq
ind_rel = find(hmp1.rel_freq > threshold_min_pres & hmp1.rel_freq < threshold_max_pres);

D1_1 = nan(N, 1);
D2_1 = nan(N, 1);
Q_1 = nan(N, 1);

data = hmp1.data;
datad = hmp1.datad;
M = hmp1.num_samples;

hmp1.B_i_graph = cell(N, 1);
hmp1.s_i = cell(N, 1);
hmp1.dist_i = cell(N, 1);

for ii = 1:length(ind_rel)
    i = ind_rel(ii);
    i
    % Modularity
    data_i = data; data_i(i, :) = [];
    
    % Renormalization
    data_i = data_i./sum(data_i);
    
    % Calculating the network
    distances_i = pdist(data_i', @BC_pdist2);
    [cdf_dist, dist] = ecdf(distances_i);
    dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last')); % Finding distance threshold
    B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
    B_i(1:M+1:end) = 0;
    B_i_graph = graph(B_i);
    
    % Calculating the modularity
    s_i = double(datad(i, :)); s_i(s_i==0) = -1; s_i = s_i';
    Q_1(i) = modularity_guy(B_i, s_i);
    
    % Saving the result
    hmp1.B_i_graph{i} = B_i_graph;
    hmp1.s_i{i} = s_i;
    
    
    % D1 & D2
    ind_pres = data(i, :) ~= 0;
    data_pres = data(:, ind_pres);
    data_abs = data(:, not(ind_pres));
    data_pres(i, :) = [];
    data_abs(i, :) = [];
    
    % Renormalization
    data_pres = data_pres./sum(data_pres);
    data_abs = data_abs./sum(data_abs);
    
    hmp1.dist_i{i} = pdist2(data_i', data_i', @BC_pdist2);
    
    D1_1(i) = sum(sum(pdist2(data_pres', data_abs', @BC_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));
    D2_1(i) = pdist2(mean(data_pres'), mean(data_abs'), @BC_pdist2);    
    
end

hmp1.Q = Q_1;
hmp1.D1 = D1_1;
hmp1.D2 = D2_1;

%% Calculating impact only after
% Do this only for species in a specific region?

% Indecies of rel_freq
% ind_rel = find(hmp2.rel_freq > threshold_min_pres & hmp2.rel_freq < threshold_max_pres);
% 
% D1_2 = nan(N, 1);
% D2_2 = nan(N, 1);
% Q_2 = nan(N, 1);
% 
% data = hmp2.data;
% datad = hmp2.datad;
% M = hmp2.num_samples;
% 
% hmp2.B_i_graph = cell(N, 1);
% hmp2.s_i = cell(N, 1);
% hmp2.dist_i = cell(N, 1);
% 
% for ii = 1:length(ind_rel)
%     i = ind_rel(ii);
%     i
%     % Modularity
%     data_i = data; data_i(i, :) = [];
%     
%     % Renormalization
%     data_i = data_i./sum(data_i);
%     
%     % Calculating the network
%     distances_i = pdist(data_i', @BC_pdist2);
%     [cdf_dist, dist] = ecdf(distances_i);
%     dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last')); % Finding distance threshold
%     B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
%     B_i(1:M+1:end) = 0;
%     B_i_graph = graph(B_i);
%     
%     % Calculating the modularity
%     s_i = double(datad(i, :)); s_i(s_i==0) = -1; s_i = s_i';
%     Q_2(i) = modularity_guy(B_i, s_i);
%     
%     % Saving the result
%     hmp2.B_i_graph{i} = B_i_graph;
%     hmp2.s_i{i} = s_i;
%     
%     
%     % D1 & D2
%     ind_pres = data(i, :) ~= 0;
%     data_pres = data(:, ind_pres);
%     data_abs = data(:, not(ind_pres));
%     data_pres(i, :) = [];
%     data_abs(i, :) = [];
%     
%     % Renormalization
%     data_pres = data_pres./sum(data_pres);
%     data_abs = data_abs./sum(data_abs);
%     
%     hmp2.dist_i{i} = pdist2(data_i', data_i', @BC_pdist2);
%     
%     D1_2(i) = sum(sum(pdist2(data_pres', data_abs', @BC_pdist2)))./(sum(ind_pres).*sum(not(ind_pres)));
%     D2_2(i) = pdist2(mean(data_pres'), mean(data_abs'), @BC_pdist2);    
%     
% end
% 
% hmp2.Q = Q_2;
% hmp2.D1 = D1_2;
% hmp2.D2 = D2_2;

%% Statistical testing
% Removing nans from I_12
I_12_new = I_12;

OTU_12 = hmp1.OTU;
ind = isnan(I_12_new);
OTU_12(ind) = [];
I_12_new(ind) = [];

% Finding OTU names of top X% I_12 and shuffled
I_12_num_OTU = round((length(I_12_new).*threshold_stat_test));

[~, temp] = sort(I_12_new, 'descend');
temp = temp(1:I_12_num_OTU);

I_12_new = I_12_new(temp);
OTU_12 = OTU_12(temp);

% Finding OTU names of top Q/D1/D2
% Q:
Q = hmp1.Q;
OTU_Q = hmp1.OTU;
ind = isnan(Q);

OTU_Q(ind) = [];
Q(ind) = [];

% Finding OTU names of top X% Q 
Q_num_OTU = round((length(Q).*threshold_stat_test));

[~, temp] = sort(Q, 'descend');
temp = temp(1:Q_num_OTU);

Q = Q(temp);
OTU_Q = OTU_Q(temp);


% D1:
D1 = hmp1.D1;
OTU_D1 = hmp1.OTU;
ind = isnan(D1);

OTU_D1(ind) = [];
D1(ind) = [];

% Finding OTU names of top X% Q 
D1_num_OTU = round((length(D1).*threshold_stat_test));

[~, temp] = sort(D1, 'descend');
temp = temp(1:D1_num_OTU);

D1 = D1(temp);
OTU_D1 = OTU_D1(temp);


% D2:
D2 = hmp1.D2;
OTU_D2 = hmp1.OTU;
ind = isnan(D2);

OTU_D2(ind) = [];
D2(ind) = [];

% Finding OTU names of top X% Q 
D2_num_OTU = round((length(D2).*threshold_stat_test));

[~, temp] = sort(D2, 'descend');
temp = temp(1:D2_num_OTU);

D2 = D2(temp);
OTU_D2 = OTU_D2(temp);


% Intersecting tables and comparing
% Q:
OTU_Q_intersect = intersect(OTU_Q, OTU_12);
if isempty(OTU_Q_intersect)
    OTU_Q_intersect_num = 0;
else
    OTU_Q_intersect_num = length(OTU_Q_intersect);
end

% D1:
OTU_D1_intersect = intersect(OTU_D1, OTU_12);
if isempty(OTU_D1_intersect)
    OTU_D1_intersect_num = 0;
else
    OTU_D1_intersect_num = length(OTU_D1_intersect);
end

% D2:
OTU_D2_intersect = intersect(OTU_D2, OTU_12);
if isempty(OTU_D2_intersect)
    OTU_D2_intersect_num = 0;
else
    OTU_D2_intersect_num = length(OTU_D2_intersect);
end


%% Shuffeling
OTU_Q_intersect_num_shuffle = nan(num_reshuffling, 1);
OTU_D1_intersect_num_shuffle = nan(num_reshuffling, 1);
OTU_D2_intersect_num_shuffle = nan(num_reshuffling, 1);

for i = 1:num_reshuffling
    i
    % Removing nans from I_12
    I_12_shuffle = I_12;
    OTU_12_shuffle = hmp1.OTU;
    ind = isnan(I_12_shuffle);
    OTU_12_shuffle(ind) = [];
    I_12_shuffle(ind) = [];
    
    % Shuffeling
    temp = randperm(length(I_12_shuffle));
    I_12_shuffle = I_12_shuffle(temp);
    
    
    % Finding OTU names of top X% I_12 and shuffled
    I_12_num_OTU = round((length(I_12_shuffle).*threshold_stat_test));
    
    [~, temp] = sort(I_12_shuffle, 'descend');
    temp = temp(1:I_12_num_OTU);
    
    I_12_shuffle = I_12_shuffle(temp);
    OTU_12_shuffle = OTU_12_shuffle(temp);
    
    
%     OTU_12_shuffle
    % Intersecting tables and comparing
    % Q:
    % Also shuffling Q/D1/D2:
%     OTU_Q = OTU_Q(randperm(length(OTU_Q)));
    OTU_Q_intersect_shuffle = intersect(OTU_Q, OTU_12_shuffle);
    if isempty(OTU_Q_intersect_shuffle)
        OTU_Q_intersect_num_shuffle(i) = 0;
    else
        OTU_Q_intersect_num_shuffle(i) = length(OTU_Q_intersect_shuffle);
    end
    
    % D1:
%     OTU_D1 = OTU_D1(randperm(length(OTU_D1)));
    OTU_D1_intersect_shuffle = intersect(OTU_D1, OTU_12_shuffle);
    if isempty(OTU_D1_intersect_shuffle)
        OTU_D1_intersect_num_shuffle(i) = 0;
    else
        OTU_D1_intersect_num_shuffle(i) = length(OTU_D1_intersect_shuffle);
    end
    
    % D2:
%     OTU_D2 = OTU_D2(randperm(length(OTU_D2)));
    OTU_D2_intersect_shuffle = intersect(OTU_D2, OTU_12_shuffle);
    if isempty(OTU_D2_intersect_shuffle)
        OTU_D2_intersect_num_shuffle(i) = 0;
    else
        OTU_D2_intersect_num_shuffle(i) = length(OTU_D2_intersect_shuffle);
    end
end

% Calculating pvalues
Q_pvalue = sum(double(OTU_Q_intersect_num_shuffle >= OTU_Q_intersect_num))./num_reshuffling
D1_pvalue = sum(double(OTU_D1_intersect_num_shuffle >= OTU_D1_intersect_num))./num_reshuffling
D2_pvalue = sum(double(OTU_D2_intersect_num_shuffle >= OTU_D2_intersect_num))./num_reshuffling


%% Scatter plots
figure; 
subplot(1, 3, 1);
plot(I_12, D1_1, '.', 'MarkerSize', 20);
axis square;
[r, p] = corr(I_12, D1_1, 'rows', 'complete');
title("r = " + string(r) + ", p = " + string(p));
xlabel('I_{before-after}'); ylabel('D_1');
set(gca, 'FontSize', 20);

subplot(1, 3, 2);
plot(I_12, D2_1, '.', 'MarkerSize', 20);
axis square;
[r, p] = corr(I_12, D2_1, 'rows', 'complete');
title("r = " + string(r) + ", p = " + string(p));
xlabel('I_{before-after}'); ylabel('D_2');
set(gca, 'FontSize', 20);

subplot(1, 3, 3);
plot(I_12, Q_1, '.', 'MarkerSize', 20);
axis square;
[r, p] = corr(I_12, Q_1, 'rows', 'complete');
title("r = " + string(r) + ", p = " + string(p));
xlabel('I_{before-after}'); ylabel('Q');
set(gca, 'FontSize', 20);

%%
figure; 
plot(hmp1.rel_freq, I_12, '.', 'MarkerSize', 20);
axis square;
[r, p] = corr(I_12, hmp1.rel_freq, 'rows', 'complete');
title("r = " + string(r) + ", p = " + string(p));
xlabel('Relative freq'); ylabel('I_{before-after}');
set(gca, 'FontSize', 20);


%%
figure; 
plot(mean(hmp1.data, 2), I_12, '.', 'MarkerSize', 20);
axis square;
[r, p] = corr(I_12, mean(hmp1.data, 2), 'rows', 'complete');
title("r = " + string(r) + ", p = " + string(p));
xlabel('Mean abundance'); ylabel('I_{before-after}');
set(gca, 'FontSize', 20);

%% 
% figure;
% subplot(1, 3, 1);
% plot(hmp1.D1, hmp2.D1, '.', 'MarkerSize', 20);
% axis square;
% [r, p] = corr(hmp1.D1, hmp2.D1, 'rows', 'complete');
% title("r = " + string(r) + ", p = " + string(p));
% xlabel('D1_{before}'); ylabel('D1_{after}');
% set(gca, 'FontSize', 20);
% 
% subplot(1, 3, 2);
% plot(hmp1.D2, hmp2.D2, '.', 'MarkerSize', 20);
% axis square;
% [r, p] = corr(hmp1.D2, hmp2.D2, 'rows', 'complete');
% title("r = " + string(r) + ", p = " + string(p));
% xlabel('D2_{before}'); ylabel('D2_{after}');
% set(gca, 'FontSize', 20);
% 
% subplot(1, 3, 3);
% plot(hmp1.Q, hmp2.Q, '.', 'MarkerSize', 20);
% axis square;
% [r, p] = corr(hmp1.Q, hmp2.Q, 'rows', 'complete');
% title("r = " + string(r) + ", p = " + string(p));
% xlabel('Q_{before}'); ylabel('Q_{after}');
% set(gca, 'FontSize', 20);




%% Figure for paper 1 - I_12 vs D1 with shaded area
figure; 
set(gcf, 'Position', [20 100 800 800], 'Units', 'centimeters');
hold on;
facealpha = 0.3;
markersize = 20;
fontsize = 30;

% Keeping only cases where D1_1 and I_12 are not nan
ind = isnan(D1_1) | isnan(I_12);
D1_1_new = D1_1(not(ind));
I_12_new = I_12(not(ind));

plot(D1_1_new, I_12_new, '.k', 'MarkerSize', markersize);
xlabel('D_1');
ylabel('L', 'Rotation', 0);
set(gca, 'FontSize', fontsize);
axis square;

% Finding value of top percent I and shading:
I_12_num_OTU = round((length(I_12_new).*threshold_stat_test));
[temp, ~] = sort(I_12_new, 'descend');
temp = temp(I_12_num_OTU);

ymin = temp;
ymax = max(I_12);
xmin = min(D1_1);
xmax = max(D1_1);


patch('XData', [xmin, xmax, xmax, xmin], 'YData', [ymin, ymin, ymax, ymax],...
    'FaceColor', bcolor, 'EdgeColor', 'none', 'FaceAlpha', facealpha);
% area([xmin, xmax], [ymin, ymax], 'EdgeColor', 'none', 'FaceColor', bcolor, 'FaceAlpha', facealpha);

% Finding value of top percent D1 and shading:
% D1_1_new = D1_1;
% D1_1_new(isnan(D1_1_new)) = [];
D1_1_num_OTU = round((length(D1_1_new).*threshold_stat_test));
[temp, ~] = sort(D1_1_new, 'descend');
temp = temp(D1_1_num_OTU);

ymin = min(I_12);
ymax = max(I_12);
xmin = temp;
xmax = max(D1_1);

patch('XData', [xmin, xmax, xmax, xmin], 'YData', [ymin, ymin, ymax, ymax],...
    'FaceColor', rcolor, 'EdgeColor', 'none', 'FaceAlpha', facealpha);
% area([xmin, xmax], [ymin, ymax], 'EdgeColor', 'none', 'FaceColor', bcolor, 'FaceAlpha', facealpha);

plot(D1_1, I_12, '.k', 'MarkerSize', markersize);

xlim([min(D1_1), max(D1_1)]);
ylim([min(I_12), max(I_12)]);

%% Fisher test
% Removing nan
% ind = isnan(D1_1) | isnan(I_12);
% D1_1_new = D1_1(not(ind));
% I_12_new = I_12(not(ind));

% Testing if in top x%
% threshold_stat_test
% D1:
D1_test = zeros(length(D1_1_new), 1);

D1_num_OTU = round((length(D1_1_new).*threshold_stat_test));
[~, temp] = sort(D1_1_new, 'descend');
temp = temp(1:D1_num_OTU);

D1_test(temp) = 1;

% I12:
I_12_test = zeros(length(I_12_new), 1);

I12_num_OTU = round((length(I_12_new).*threshold_stat_test));
[~, temp] = sort(I_12_new, 'descend');
temp = temp(1:I12_num_OTU);

I_12_test(temp) = 1;

clear x;
% Creating contigency table
x(1, 1) = sum((D1_test==0) & (I_12_test==0));
x(1, 2) = sum((D1_test==0) & (I_12_test==1));
x(2, 1) = sum((D1_test==1) & (I_12_test==0));
x(2, 2) = sum((D1_test==1) & (I_12_test==1));

[h, p, stats] = fishertest(x);

%% Histogram of I_12_cell
% Choosting i's with at least 20 k_i. Top and bottom
ind = k_i>=20;
k_i_temp = k_i(ind);
I_12_temp = I_12(ind);
I_12_cell_temp = I_12_cell(ind);

[~, ind1] = min(I_12_temp);
[~, ind2] = max(I_12_temp);
[~, ind3] = sort(I_12_temp);
ind3 = ind3(round(end/2));

I_12_cell_temp_1 = I_12_cell_temp(ind1);
I_12_cell_temp_2 = I_12_cell_temp(ind2);
I_12_cell_temp_3 = I_12_cell_temp(ind3);

figure;
hold on;
num_bins = 10;
% histogram(I_12_cell_temp_1{:}, num_bins, 'Normalization', 'pdf');
% histogram(I_12_cell_temp_2{:}, num_bins, 'Normalization', 'pdf');
% histogram(I_12_cell_temp_3{:}, num_bins, 'Normalization', 'pdf');

% boxplot(I_12_cell_temp_1{:}, I_12_cell_temp_2{:}, I_12_cell_temp_3{:});

g = [ones(1, length(I_12_cell_temp_1{:})), 2.*ones(1, length(I_12_cell_temp_2{:})), 3.*ones(1, length(I_12_cell_temp_3{:}))];
x = [I_12_cell_temp_1{:}, I_12_cell_temp_3{:}, I_12_cell_temp_2{:}];
boxplot(x, g);


xlabel('L^i_k');
set(gca, 'FontSize', 20);
% legend('Low L^i', 'High L^i', 'Midd L^i');
