function [shuffle_data1, shuffle_data2] = shuffle_data_for_two(data1, data2)
% Crearting a shuffled matrix
[num_genes,num_cells] = size(data1);
shuffle_data1 = nan(size(data1));
shuffle_data2 = nan(size(data2));
for j = 1:num_genes
    temp_ind = randperm(num_cells);
    temp1 = data1(j,:);
    shuffle_data1(j,:) = temp1(temp_ind);
    temp_ind = randperm(num_cells);
    temp2 = data2(j,:);
    shuffle_data2(j,:) = temp2(temp_ind);
end

end