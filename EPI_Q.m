function Q = EPI_Q(S, threshold_net)
%EPI_Q Emperical Presence-Impact of a cohort by using the modularity meausre, Q
%   Q = EPI(data) returns the Q EPI values of each of the N species in a
%   cohort of M samples, using the Bray-Curtis dissimilarity measure. 
%   S is a cohort of samples, represented by a N-by-M matrix of abundances.
%   threshold_net is the distance threshold for the constructing of the
%   sample-to-sample network, by percentile
%   Q is an N-by-1 vector of Q values.

% Initialization
[N, M] = size(S);
S_01 = double(S>0);
Q = nan(N, 1);

for i = 1:N
    % If the species is always present/absent, Q is undefined
    if sum(S_01(i, :), 2) ~= 0 || sum(S_01(i, :), 2) ~= M
        
        % Removing the i species
        S_i = S;
        S_i(i, :) = [];
        
        % Normalizing
        S_i = S_i./sum(S_i);
        
        % Building the network
        distances_i = pdist(S_i', @BC_pdist2);
        [cdf_dist, dist] = ecdf(distances_i);
        dist_threshold = dist(find(cdf_dist<=threshold_net, 1, 'last'));
        B_i = squareform(distances_i, 'tomatrix')<=dist_threshold;
        B_i(1:M+1:end) = 0;
        s_i = double(S_01(i, :)); 
        s_i(s_i==0) = -1; 
        s_i = s_i';
        
        % Calculating
        Q(i) = modularity(B_i, s_i);
    end
end

function Dpdist = BC_pdist2(ZI, ZJ)
    BC = @(x, y) sum((abs(x - y)))./(sum(x) + sum(y));
    m2 = size(ZJ, 1);
    Dpdist = nan(m2, 1);
    x = ZI;
    for k = 1:m2
        y = ZJ(k, :);
        Dpdist(k) = BC(x, y);
    end
end

function Qmod = modularity(B, s)
    B_graph = graph(B);
    d = degree(B_graph); % Degree of each sample
    q = B_graph.numedges;
    Qmod = (s' * (B - (d*d')/(2*q)) * s) / (4*q);
end

end