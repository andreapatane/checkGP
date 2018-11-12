function [Y,std_Y] = generate_Y_values(X,s,n,noise_std,std_Y)
if nargin < 5
    std_Y = [];
end
% generate nonlinear function of X as the sum of all cross-products
J =  1:s;    % select the first s variables
Y = zeros(n,1);
for i=1:s
	for j=1:i-1
		Y = Y + X(:,J(i)) .* X(:,J(j));
	end
end

% add some noise with known standard deviation
Y =  Y + randn(n,1) * noise_std;

% normalize to unit standard deviation
if isempty(std_Y)
    std_Y =  std(Y);
end
Y = Y / std_Y;
%rescaling
Y = Y /10;

end