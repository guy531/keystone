function D2 = BC_pdist2(ZI, ZJ)
% Calculating Bray-Curtis distance for pdist2 function

% Definition from here: http://www.econ.upf.edu/~michael/stanford/maeb5.pdf

% https://metarabbit.wordpress.com/2018/01/16/bray-curtis-dissimilarity-on-relative-abundance-data-is-the-manhattan-distance-aka-l%E2%82%81-distance/
% For normalized data Bray-Curtis = Manhattan distance
% Manhattan distance: D = sum(abs(x-y))
% Bray Curtis: 

% ZI is a 1-by-n vector containing a single observation.
% ZJ is an m2-by-n matrix containing multiple observations. distfun must accept a matrix ZJ with an arbitrary number of observations.
% D2 is an m2-by-1 vector of distances, and D2(k) is the distance between observations ZI and ZJ(k,:).
% If your data is not sparse, you can generally compute distance more quickly by using a built-in distance instead of a function handle.

BC = @(x, y) sum((abs(x - y)))./(sum(x) + sum(y));


m2 = size(ZJ, 1);
D2 = nan(m2, 1);
x = ZI;
for k = 1:m2
    y = ZJ(k, :);
    D2(k) = BC(x, y);
end

end