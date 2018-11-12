function [a,b] = line_thru_points(x_l,f_l,x_u,f_u)

b = (f_l - f_u) / (x_l - x_u);

a = f_l -b*x_l;

end