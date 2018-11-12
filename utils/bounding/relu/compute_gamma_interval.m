function [gamma_l,gamma_u] = compute_gamma_interval(beta_l,beta_u)
% compute minimum and maximum of the function:
% gamma(phi) = acos((c1+c2cos(cartesian_of(phi) \dot cartesian_of(psi)))/(c1+c2) )
% subject to the constraint that phi is angular coordiante in the surface of S_m

[gamma_l,gamma_u] = min_max_of_gamma(beta_l,beta_u);



end