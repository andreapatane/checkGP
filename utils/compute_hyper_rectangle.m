function [x_L, x_U] = compute_hyper_rectangle(epsilon,testPoint, pixel2Modify,dataSet)


x_L = testPoint;
x_U = testPoint;

for ii = 1:length(pixel2Modify)
    currPixel = pixel2Modify(ii);
    if strcmp(dataSet,'mnist')
        x_L(currPixel) = max(x_L(currPixel) - epsilon,0); % I take care that it is still a valid pixel
        x_U(currPixel) = min(x_U(currPixel) + epsilon,1.0); % I take care that it is still a valid pixel
    elseif strcmp(dataSet,'x_products')
        x_L(currPixel) = x_L(currPixel) - epsilon; % I take care that it is still a valid pixel
        x_U(currPixel) = x_U(currPixel) + epsilon; % I take care that it is still a valid pixel
    end
end




end