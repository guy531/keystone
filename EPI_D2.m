function D2 = EPI_D2(S)
%EPI_D2 Emperical Presence-Impact of a cohort by using the D2 measure
%   D2 = EPI(data) returns the D2 EPI values of each of the N species in a
%   cohort of M samples, using the Bray-Curtis dissimilarity measure. 
%   S is a cohort of samples, represented by a N-by-M matrix of abundances.
%   D2 is an N-by-1 vector of D2 values.

% Initialization
[N, M] = size(S);
S_01 = double(S>0);
D2 = nan(N, 1);

for i = 1:N
    % If the species is always present/absent, D2 is undefined
    if sum(S_01(i, :), 2) ~= 0 || sum(S_01(i, :), 2) ~= M
        
        % Dividing into the two groups
        ind_pres = S_01(i, :) ~= 0;
        S_pres = S(:, ind_pres);
        S_abs = S(:, not(ind_pres));
        
        % Removing the i species
        S_pres(i, :) = [];
        S_abs(i, :) = [];
        
        % Normalizing
        S_pres = S_pres./sum(S_pres);
        S_abs = S_abs./sum(S_abs);
        
        % Calculating D2
        D2(i) = pdist2(mean(S_pres'), mean(S_abs'), @BC_pdist2);
        
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

end