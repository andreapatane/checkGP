function one_hot_labels = one_hot_encoding(labels)
one_hot_labels = -0.1*ones(size(labels,1),10);
for jj = 1:length(labels)
    one_hot_labels(jj,labels(jj)+1) = 0.9;
end
end