function shuffle_data_output = shuffle_data_ravid_2(data)
% Crearting a shuffled with weight to samples by gene prevalence 
% For 1/0 data only

[num_genes,num_cells] = size(data);
shuffle_data_output = nan(size(data));

% Calculating weight of each sample
weights = sum(double(data>0));

% data_shuffle_temp = shuffle_data(


for j = 1:num_genes
    temp2 = data(j,:);
    % Random choosing of samples
    num_gene_present = sum(double(temp2>0));
    temp2 = shuffle_data(data(j, :));
    temp = zeros(size(temp2));
    ind = datasample(1:num_cells, num_gene_present, 'Weights', weights, 'Replace', false);
    temp(ind) = temp2(ind);
    shuffle_data_output(j,:) = temp;
end


end