function [data_matrix,data_norm] = to_norm1(data_matrix)

data_norm = zeros(size(data_matrix,1),1);
for ii = 1:size(data_matrix,1)
    data_norm(ii) = norm(data_matrix(ii,:));
    data_matrix(ii,:) = data_matrix(ii,:)/data_norm(ii);
end
end