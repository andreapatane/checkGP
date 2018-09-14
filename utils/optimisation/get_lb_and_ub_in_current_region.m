function [c_lb,c_ub,phi_m] = get_lb_and_ub_in_current_region(mode,kernel,mult_coefficients,c_x_l,...
    c_x_u,theta_vec,sigma_prior,training_data,outputIdx,testPoint,r_star,R,training_data_cartesian,testPoint_cartesian,phi_m_star,depth)
%getting lower and upper bounds in the current region. The fucntion is
%basically just a wrapper, getting into account the various modes and
%kernels used.

if strcmp(kernel,'N-ReLU')
    [c_x_l,c_x_u] = compute_angle_space_hyper_rectangle(c_x_l,c_x_u);
end

phi_m = [];
if strcmp(mode,'mu')
    if strcmp(kernel,'sqe')
        
        [c_lb,c_ub] = compute_lower_bound_mu_sqe(mult_coefficients,c_x_l,c_x_u,...
            theta_vec,sigma_prior,training_data,outputIdx,0,testPoint);
        
    elseif strcmp(kernel,'N-ReLU')
                 [c_lb,c_ub,phi_m] = compute_lower_bound_mu_ReLU(mult_coefficients,c_x_l,c_x_u,...
                     theta_vec,sigma_prior,training_data,outputIdx,testPoint,phi_m_star,training_data_cartesian,depth);
    end
elseif strcmp(mode,'xi')
    if strcmp(kernel,'sqe')
        [c_lb,c_ub] = compute_lower_bound_xi_sqe(mult_coefficients,c_x_l,c_x_u,...
            theta_vec,sigma_prior,training_data,outputIdx,0,testPoint,r_star,R);
    elseif strcmp(kernel,'N-ReLU')
        [c_lb,c_ub,phi_m] = compute_lower_bound_xi_ReLU(mult_coefficients,c_x_l,c_x_u,...
            theta_vec,sigma_prior,training_data,training_data_cartesian,...
            testPoint,testPoint_cartesian,r_star,phi_m_star,depth);
    end
elseif strcmp(mode,'sigma')
    if strcmp(kernel,'sqe')
        [c_lb,c_ub] = compute_lower_bound_sigma_sqe(mult_coefficients,c_x_l,c_x_u,...
            theta_vec,sigma_prior,training_data,outputIdx,0,testPoint,R);
    elseif strcmp(kernel,'N-ReLU')
        [c_lb,c_ub,phi_m] = compute_lower_bound_sigma_ReLU(mult_coefficients,c_x_l,c_x_u,...
            theta_vec,sigma_prior,training_data,training_data_cartesian,phi_m_star,depth);
    end
    
end