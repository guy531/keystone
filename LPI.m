function [L, k_i] = LPI(S1, S2)
%LPI Longtitudinal Presence-Impact of all the species using the samples in 
%   the first and second collection (S1 and S2)
%   L = LPI(data) returns the L values of each of the N species in the two
%   cohorts of M subjects using the Bray-Curtis dissimilarity measure. 
%   S1 is the cohort of samples in the first collection, and S2 is of the 
%   second collection, represented by a N-by-M matrices of abundances.
%   L is an N-by-1 vector of L values. If the species had no reversed
%   presence-state in all the subjects, its L value is NaN.
%   [L, k_i] = LPI(S1, S2) also returns the number of subjects for which
%   the L value of each species was calculated for (the number of cases
%   where the species had a reversed presence-state).

[N, M] = size(S1);
L = nan(N, 1);
k_i = nan(N, 1); % Index number of samples
for i = 1:N
    L_i = [];
    k = 0;
    for j = 1:M
        data1 = S1(:, j);
        data2 = S2(:, j);
        datad1 = double(data1>0);
        datad2 = double(data2>0);
        % Checking if the presence-state of species i is reversed
        if datad1(i) ~= datad2(i)
            % Removing the species
            data1_i = data1;
            data1_i(i) = [];
            data2_i = data2;
            data2_i(i) = [];
            
            % Normalizing
            data1_i = data1_i./sum(data1_i);
            data2_i = data2_i./sum(data2_i);
            
            % Calculating the distance before/after
            k = k + 1;
            L_i(k) = pdist2(data1_i', data2_i', @BC_pdist2);
        end
    end
    k_i(i) = k;
    L(i) = mean(L_i);
end

function Dpdist = BC_pdist2(ZI, ZJ)
    BC = @(x, y) sum((abs(x - y)))./(sum(x) + sum(y));
    m2 = size(ZJ, 1);
    Dpdist = nan(m2, 1);
    x = ZI;
    for kk = 1:m2
        y = ZJ(kk, :);
        Dpdist(kk) = BC(x, y);
    end
end

end