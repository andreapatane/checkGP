function data_matrix = prioritize_pixel2mod(data_matrix,pix2Mod)
%I put feature pixels at the beginning of the image vectors. This makes
%things easier in the spherical space.

m = size(data_matrix,2);
for ii = 1:size(data_matrix,1)
    currPoint = [data_matrix(ii,ismember(1:m,pix2Mod)) , data_matrix(ii,~ismember(1:m,pix2Mod))];
    data_matrix(ii,:) = currPoint;
end


end
