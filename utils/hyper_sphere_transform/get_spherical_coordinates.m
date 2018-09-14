function data_matrix_angular = get_spherical_coordinates(data_matrix,pixel2Modify)
%function for spherical coordinates of testPoint. The pixels that will be
%modified will be put at the end
m = size(data_matrix,2);
data_matrix_angular = zeros(size(data_matrix,1),m-1);

for ii = 1:size(data_matrix,1)
    currPoint = [data_matrix(ii,~ismember(1:m,pixel2Modify)) , data_matrix(ii,ismember(1:m,pixel2Modify)) ];


    
    for jj = 1:(m-2)
        if norm(currPoint(jj:end))
            data_matrix_angular(ii,jj) = acot(currPoint(jj)/norm(currPoint( (jj + 1):end)));
        else
            data_matrix_angular(ii,jj) = 0;
        end

    end
    jj = m-1;
    
    if currPoint(jj+1) >0
        data_matrix_angular(ii,jj) = acos(currPoint(jj)/norm(currPoint(jj:end)));
    elseif currPoint(jj+1) < 0
        data_matrix_angular(ii,jj) = - acos(currPoint(jj)/norm(currPoint(jj:end)));
    else
        data_matrix_angular(ii,jj) = 0;
    end
    
    
end




end