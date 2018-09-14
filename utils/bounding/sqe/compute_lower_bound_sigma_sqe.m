function [lb_on_sigma,ub_on_sigma] = compute_lower_bound_sigma_sqe(R_inv,x_L,x_U,theta_vec,sigma_prior,training_data,...
                                                outputIdx,debug,testPoint,R)
%
%Function for variance computation on SQE kernel
%
     
                                            

%scaling everything by a-priori variance
R_inv = sigma_prior * R_inv;
R  = R / sigma_prior;

n = size(training_data,1);
m = size(training_data,2);

z_i_L_vec = zeros(1,n);
z_i_U_vec = zeros(1,n);
                                            
B_il =  zeros(n,n);
a_il_sum = 0;

%precomputing all the z-intervals
for ii = 1:n
    [z_i_L_vec(ii),z_i_U_vec(ii)] = compute_z_interval(training_data(ii,:),x_L,x_U,theta_vec);
end                                           


%iterating over all the couples of training points (things are symmetric though).
for ii = 1:n
    for ll = 1:ii
        z_il_L = z_i_L_vec(ii) + z_i_L_vec(ll);
        z_il_U = z_i_U_vec(ii) + z_i_U_vec(ll);
        [a_il,B_il(ii,ll)] = compute_linear_under_approx(R_inv(ii,ll),z_il_L,z_il_U);        
        if ll < ii
            a_il = 2*a_il;
            B_il(ll,ii) = B_il(ii,ll);
        end
        a_il_sum = a_il_sum + a_il;
    end
end


%using  lower approximations to define hessian and f vector and bias for
%quadratic program
B_sum = sum(B_il,2);

C = 0;
for ii = 1:n
    C = C + B_sum(ii) * dot(theta_vec,training_data(ii,:).^2);
end
C = 2*C;
H = 4*sum(B_sum) * theta_vec;

f = zeros(m,1);
for  jj = 1:m
    f(jj) = -4*theta_vec(jj) * dot(B_sum,training_data(:,jj));
end


[x_star, f_val] = separate_quadprog(H',f,x_L,x_U);
f_val = f_val + a_il_sum + C;
lb = f_val;

%getting actual prediction
lb_on_sigma = sigma_prior*get_actual_prediction(training_data,x_star',theta_vec,R);
ub_on_sigma = sigma_prior*(1 - lb);
end



function sigma_test = get_actual_prediction(training_data,x_star,theta_vec,R)

[Kstar,Kstarstar] = compute_kernel_on_test(x_star,training_data,'sqe',0,1,theta_vec,{});
[~,sigma_test] = make_predictions_on_test(Kstar,[],Kstarstar,R);
end

