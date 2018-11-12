function data_matrix_angular = get_spherical_coordinates(data_matrix,pixel2Modify)
%function for spherical coordinates of testPoint. The pixels that will be
%modified will be put at the end


if nargin < 2
    pixel2Modify = [];
end

m = size(data_matrix,2);
data_matrix_angular = zeros(size(data_matrix,1),m-1);

for ii = 1:size(data_matrix,1)
    currPoint = [data_matrix(ii,~ismember(1:m,pixel2Modify)) , data_matrix(ii,ismember(1:m,pixel2Modify)) ];
    % the only thing I need to make the following computations safely, is
    % that the last pixel is differnt from zero. I check if it is,
    % otherwise I add a tiny jitter to it.
    %if currPoint(end) == 0
        %currPoint(end) = 0.00001;
    %end
    for jj = 1:(m-2)
        %if currPoint(jj) ~= 0
        if norm(currPoint( (jj + 1) :end)) > 0
            data_matrix_angular(ii,jj) = acos(currPoint(jj)/norm(currPoint(jj:end)));
            %data_matrix_angular(ii,jj) = acot(currPoint(jj)/norm(currPoint( (jj + 1):end)));
            %disp('here')
        elseif currPoint(jj) >= 0
            data_matrix_angular(ii,jj) = 0;
            %disp('here1')
        else
            data_matrix_angular(ii,jj) = pi;
            %disp('here2')
        end
        %else
        %    data_matrix_angular(ii,jj) = 0.5*pi;
        %end
    end
    jj = m-1;
    
    if currPoint(jj+1) >0
        data_matrix_angular(ii,jj) = acos(currPoint(jj)/norm(currPoint(jj:end)));
    elseif currPoint(jj+1) < 0
        data_matrix_angular(ii,jj) = 2*pi - acos(currPoint(jj)/norm(currPoint(jj:end)));
        %disp('here3')
    else
        data_matrix_angular(ii,jj) = 0;
        %disp('here4')
    end
    
    
end




end