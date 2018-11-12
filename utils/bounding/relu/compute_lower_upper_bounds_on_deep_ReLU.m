function [a_L,b_L,a_U,b_U,beta_l,beta_u] =...
    compute_lower_upper_bounds_on_deep_ReLU(psi,phi_L,phi_U,sigma_b_2,sigma_w_2,m,depth)



[beta_l,beta_u] = compute_beta_interval(psi,phi_L,phi_U,sigma_b_2,sigma_w_2);

k1 = sigma_b_2;
a_L_vec = zeros(depth,1);
b_L_vec = zeros(depth,1);
a_U_vec = zeros(depth,1);
b_U_vec = zeros(depth,1);
K_diag = sigma_b_2 + (sigma_w_2/m);

for hh = 1:depth
    
    [gamma_l,gamma_u] = compute_gamma_interval(beta_l,beta_u);
    k2 = (sigma_w_2/(2*pi) ) * K_diag;
    
    [a_L_vec(hh),b_L_vec(hh)] = compute_linear_under_bound_ReLU(1,gamma_l,gamma_u,beta_l,beta_u,k1,k2);
    [a_U_vec(hh),b_U_vec(hh)] = compute_linear_under_bound_ReLU(-1,gamma_l,gamma_u,beta_l,beta_u,k1,k2);
    a_U_vec(hh) = - a_U_vec(hh);
    b_U_vec(hh) = - b_U_vec(hh);
    
    if hh < depth
        K_diag =  sigma_b_2 + (sigma_w_2/2)* K_diag;
        a_L_vec(hh) = a_L_vec(hh)/K_diag;
        a_U_vec(hh) = a_U_vec(hh)/K_diag;
        b_L_vec(hh) = b_L_vec(hh)/K_diag;
        b_U_vec(hh) = b_U_vec(hh)/K_diag;
    else
        K_diag = 1;
    end
    [beta_l,beta_u] = min_max_of_kernel_function(gamma_l,gamma_u,k1,k2,K_diag);
end


a_L = a_L_vec(depth);
b_L = b_L_vec(depth);

a_U = a_U_vec(depth);
b_U = b_U_vec(depth);


for hh = (depth-1):-1:1
    if b_L >=0
        a_L = a_L + b_L*a_L_vec(hh);
        b_L = b_L * b_L_vec(hh);
    else
        a_L = a_L + b_L*a_U_vec(hh);
        b_L = b_L * b_U_vec(hh);
    end
    
    if b_U >=0
        a_U = a_U + b_U*a_U_vec(hh);
        b_U = b_U * b_U_vec(hh);
    else
        a_U = a_U + b_U*a_L_vec(hh);
        b_U = b_U * b_L_vec(hh);
    end
end

a_L = a_L + b_L * sigma_b_2/(sigma_b_2 + sigma_w_2/m);
b_L = b_L * (sigma_w_2/m)/(sigma_b_2 + sigma_w_2/m);

a_U = a_U + b_U * sigma_b_2/(sigma_b_2 + sigma_w_2/m);
b_U = b_U * (sigma_w_2/m)/(sigma_b_2 + sigma_w_2/m);

end