function [a,b] = line_thru_point_with_ang_coeff(x,f_x,df_x)
b = df_x;
a = - b*x + f_x;
end