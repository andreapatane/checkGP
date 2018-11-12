function [beta_l,beta_u] = min_max_of_kernel_function(gamma_l,gamma_u,k1,k2,normalisingFactor)

beta_u = k1 + k2*(sin(gamma_l) + (pi - gamma_l)*cos(gamma_l) );
beta_l = k1 + k2*(sin(gamma_u) + (pi - gamma_u)*cos(gamma_u) );

if nargin >=5
    beta_u = beta_u/normalisingFactor;
    beta_l = beta_l/normalisingFactor;
end

end