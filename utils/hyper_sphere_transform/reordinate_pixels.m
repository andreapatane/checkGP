function data = reordinate_pixels(data,pixel2Modify)


m = size(data,2) ;
for ii = 1:size(data,1)
    %putting pixels back into their original position
    for jj = 1:length(pixel2Modify)
        currIdx = m - length(pixel2Modify) + jj;
        newIdx = pixel2Modify(jj);
        auxData = [data(ii,1:newIdx-1),data(ii,currIdx) , data(ii,newIdx:end)];
        auxData(currIdx + 1 ) = [];
        data(ii,:) = auxData;  
    end 
end
end