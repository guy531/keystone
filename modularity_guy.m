function Q = modularity_guy(B, s)
% B - MXM Interaction network.
% s - Mx1 Assignment vector (1 present, -1 absent). 
% Make sure the diagonal of B is 0. 
B_graph = graph(B);
d = degree(B_graph); % Degree of each sample
q = B_graph.numedges;
Q = (s' * (B - (d*d')/(2*q)) * s) / (4*q);
end