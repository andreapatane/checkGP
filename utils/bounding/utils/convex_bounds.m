function [a_lower,b_lower,a_upper,b_upper] = convex_bounds(x_l,x_u, fun, d_fun)
% in the convex part of the function, the lower bound is given be
% the by the tangent line thru the mid point

x_c = 0.5*(x_l + x_u);
df_c = d_fun(x_c);
f_c = fun(x_c);
[a_lower,b_lower] = line_thru_point_with_ang_coeff(x_c,f_c , df_c);

%while the upper bound is given by line that connects the two function extremest
f_l = fun(x_l);
f_u = fun(x_u);
[a_upper,b_upper] = line_thru_points(x_l,f_l,x_u,f_u);
end