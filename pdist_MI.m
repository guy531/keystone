function MI = pdist_MI(data)
% Calculating distance based on MI
N = size(data,1);
MI = zeros(N,N);
for i = 1:N-1
    i
    X = data(i,:);
    for j = i:N
        Y = data(j,:);
        MI(i,j) = nmi(X,Y);
    end
end
% Converting to square matrix
MI = triu(MI) + triu(MI,1)';
% MI = 1 - MI; 
MI(1:N+1:end) = 0;
% MI = squareform(MI,'tovector');
end