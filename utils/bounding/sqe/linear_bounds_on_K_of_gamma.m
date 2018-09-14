function [a_lower, b_lower, a_upper, b_upper] = linear_bounds_on_K_of_gamma(gamma_l,gamma_u)
%find linear and upper bounds "a+ b*gamma" on k(gamma) = sin(gamma) + cos(gamma) (pi - gamma)

gamma_f = 1.11283481547935901489;


if gamma_u <= gamma_f
    [a_lower,b_lower,a_upper,b_upper] = concave_bounds(gamma_l,gamma_u,@kernelFun,@derivative_kernelFun);
elseif gamma_l >= gamma_f
    [a_lower,b_lower,a_upper,b_upper] = convex_bounds(gamma_l,gamma_u,@kernelFun,@derivative_kernelFun);
else
    [a_lower_cnc,b_lower_cnc,a_upper_cnc,b_upper_cnc] = concave_bounds(gamma_l,gamma_f,@kernelFun,@derivative_kernelFun);
    [a_lower_cvx,b_lower_cvx,a_upper_cvx,b_upper_cvx] = convex_bounds(gamma_f,gamma_u,@kernelFun,@derivative_kernelFun);
    [a_lower,b_lower ] = merge_cnc_cvx_lower(a_lower_cnc,b_lower_cnc,a_lower_cvx,b_lower_cvx,gamma_l,gamma_u,false);
    [a_upper,b_upper ] = merge_cnc_cvx_upper(a_upper_cnc,b_upper_cnc,a_upper_cvx,b_upper_cvx,gamma_l,gamma_u,false);
end

  
end

function out = kernelFun(x)
    out = sin(x) + cos(x)*(pi - x);
end

function out = derivative_kernelFun(x)
    out = sin(x) * (x-pi);
end


