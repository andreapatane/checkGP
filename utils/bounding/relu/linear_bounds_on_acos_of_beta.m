function [a_lower, b_lower, a_upper, b_upper] = linear_bounds_on_acos_of_beta(beta_l,beta_u)
beta_f = 0;



if beta_u <= beta_f
    [a_lower,b_lower,a_upper,b_upper] = convex_bounds(beta_l,beta_u,@acosFun,@acosFunDerivative);
elseif beta_l >= beta_f
    [a_lower,b_lower,a_upper,b_upper] = concave_bounds(beta_l,beta_u,@acosFun,@acosFunDerivative);
else
    [a_lower_cvx,b_lower_cvx,a_upper_cvx,b_upper_cvx] = convex_bounds(beta_l,beta_f,@acosFun,@acosFunDerivative);
    [a_lower_cnc,b_lower_cnc,a_upper_cnc,b_upper_cnc] = concave_bounds(beta_f,beta_u,@acosFun,@acosFunDerivative);
    [a_lower,b_lower ] = merge_cnc_cvx_lower(a_lower_cnc,b_lower_cnc,a_lower_cvx,b_lower_cvx,beta_l,beta_u,1);
    [a_upper,b_upper ] = merge_cnc_cvx_upper(a_upper_cnc,b_upper_cnc,a_upper_cvx,b_upper_cvx,beta_l,beta_u,1);
end

end

function out = acosFun(x)
out = acos(x);
end

function out = acosFunDerivative(x)
out =  - 1/sqrt( 1 - x^2);
end