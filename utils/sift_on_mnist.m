function [feat_pix_cell,magnitudes] = sift_on_mnist(img,maxFeatNum,maxFeatSize,vl_setup_path)

run([vl_setup_path,'vl_setup']);
row=14;
colum=14;
img = reshape(img,row,colum)';

close all



img_single = single(img);

[f,~] = vl_sift(img_single,'Levels',20,'EdgeThresh',100,'WindowSize',100) ;

%figure
%imshow(reshape(img,14,14)','InitialMagnification' ,'fit')
%hold on
%h1 = vl_plotframe(f(:,1:end)) ;
%h2 = vl_plotframe(f(:,1:end)) ;
%set(h1,'color','k','linewidth',3) ;
%set(h2,'color','y','linewidth',2) ;

[feat_pix_cell,magnitudes] = get_feat_pixels(f,row,maxFeatNum,maxFeatSize);
%for ii = 1:length(feat_pix_cell)
%    disp(feat_pix_cell{ii}') 
%end

end


function [feat_pix_cell,magnitudes] = get_feat_pixels(f,nPixs,maxFeatNum,maxFeatSize)
feat_pix_cell = cell(size(f,2),1);
magnitudes = f(3,:);
for ii = 1:length(feat_pix_cell)
    center = f([2,1],ii);
    center(2) = center(2);
    radius = f(3,ii);
    lb_x = max(round(center(1) - radius),1);
    ub_x = min(round(center(1) + radius),nPixs);
    lb_y = max(round(center(2) - radius),1);
    ub_y = min(round(center(2) + radius),nPixs);
    feat_pix_cell{ii} = [];
    for curr_x = lb_x:ub_x
        for curr_y = lb_y:ub_y
            if norm([curr_x; curr_y  ] - center ,1)  <= radius
                feat_pix_cell{ii} = [feat_pix_cell{ii}; sub2ind([nPixs,nPixs],curr_y,curr_x)];
            end
        end
    end
end
[feat_pix_cell,magnitudes] = keep_unique(feat_pix_cell,magnitudes);
if length(feat_pix_cell) > maxFeatNum
    feat_pix_cell = feat_pix_cell(1:magnitudes);
    magnitudes = magnitudes(1:magnitudes);
end
for ii = 1:length(feat_pix_cell)
    if length(feat_pix_cell{ii}) > maxFeatSize
        feat_pix_cell{ii} = feat_pix_cell{ii}(1:maxFeatSize);
    end
end
end

function [feat_pix_cell_unique,magnitudes_unique] = keep_unique(feat_pix_cell,magnitudes)
feat_pix_cell_unique = {};
magnitudes_unique = [];
for ii = 1:length(feat_pix_cell)
    flag = true;
    for jj = (ii+1):length(feat_pix_cell)
        if isequal(feat_pix_cell{ii},feat_pix_cell{jj})
            flag = false;
        end
    end
    if flag
        feat_pix_cell_unique{end+1} = feat_pix_cell{ii};
        magnitudes_unique(end+1) = magnitudes(ii);
    end
end
feat_pix_cell_unique = feat_pix_cell_unique';

end
