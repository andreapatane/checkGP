function [m,M] = max_min_of_beta(c1,c2,x_L,x_U)
%get minimum and maximum of (c1 + c2 t)/(c1+c2), for t in [x_L,x_U]. Code
%assume c1>0, c2>0

m = (c1 + c2 * x_L)/(c1+c2);
M = (c1 + c2 * x_U)/(c1+c2);

end