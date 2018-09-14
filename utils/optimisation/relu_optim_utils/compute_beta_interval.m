function [beta_l,beta_u,dot_L,dot_U] = compute_beta_interval(psi,phi_L,phi_U,sigma_b_2,sigma_w_2)

m = length(psi) + 1;

%I invert phi_m and phi_M, because the acos will invert them at some point
[dot_L,~,dot_U,~] = min_max_of_dot_product_hypersphere(psi,phi_L,phi_U);

[beta_l,beta_u] = max_min_of_beta(sigma_b_2,sigma_w_2/m,dot_L,dot_U);




end