function data_matrix = get_cartesian_coordinates(data_matrix_angular,pixel2Modify)
%function for obtaining back cartesian coordinates of testPoint.
%The pixels that were selected for modification will be put back at the
%original position.

m = size(data_matrix_angular,2) + 1 ;
data_matrix = zeros(size(data_matrix_angular,1),m);

for ii = 1:size(data_matrix,1)
    
    %transorming spherical coordiantes into cartesian coordinates
    auxVal = 1.0;
    for jj = 1:(m-1)

        data_matrix(ii,jj) = auxVal * cos(data_matrix_angular(ii,jj));
        auxVal = auxVal*sin(data_matrix_angular(ii,jj)); 
    end
    data_matrix(ii,jj) = auxVal;


    data_matrix(ii,:) = data_matrix(ii,:)/norm(data_matrix(ii,:)); %accounting for numerical issues...
    %putting pixels back into their original position
    for jj = 1:length(pixel2Modify)
        currIdx = m - length(pixel2Modify) + jj;
        newIdx = pixel2Modify(jj);
        auxData = [data_matrix(ii,1:newIdx-1),data_matrix(ii,currIdx) , data_matrix(ii,newIdx:end)];
        auxData(currIdx + 1 ) = [];
        data_matrix(ii,:) = auxData;

        
    end
    
    
    
    
    
end




end