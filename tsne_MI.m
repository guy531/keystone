function D2 = tsne_MI(ZI,ZJ)
[m,~] = size(ZJ);
D2 = zeros(m,1);
for i = 1:m
    D2(i) = 1 - nmi(ZI,ZJ(i,:));
end
end