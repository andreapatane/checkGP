function [lb,ub] = compute_lower_bound_mu_sqe(trainedSystem,x_L,x_U,theta_vec,sigma_prior,training_data,outputIdx)
%computation for lower and upper bound on the optimisation problem for mu^0
%using sqe kernel



%getting intersting output index only
y_i_vec = trainedSystem(:,outputIdx);
%scaling output with a priori variance
y_i_vec = y_i_vec*sigma_prior;

n = size(training_data,1);


m = size(training_data,2);
a_i_sum = 0;
b_i_vec = zeros(n,1);
%looping over the training set
for ii = 1:n
    y_i = y_i_vec(ii);
    x_i = training_data(ii,:);
    [z_i_L, z_i_U] = compute_z_interval(x_i,x_L,x_U,theta_vec);
    [a_i,b_i] = compute_linear_under_approx(y_i,z_i_L,z_i_U);
    b_i_vec(ii) = b_i;

    
    a_i_sum = a_i_sum + a_i;
end

%defining matrix and vector for separable quadratic optimisation
H = 2*sum(b_i_vec)*theta_vec;

f = zeros(m,1);
for jj = 1:m
    f(jj) = - 2 * theta_vec(jj)*dot(training_data(:,jj),b_i_vec);
end
C = 0;
for ii = 1:n
    for jj = 1:m
        C = C + theta_vec(jj) * b_i_vec(ii) *training_data(ii,jj)^2;
    end
end

[x_star, f_val] = separate_quadprog(H,f,x_L,x_U);
f_val = f_val + a_i_sum + C;


%computing actual prediction on optimal point for the lower bound function
ub = get_actual_prediction(training_data,x_star',n,theta_vec,y_i_vec);
lb = f_val;


end




function mu = get_actual_prediction(training_data,testPoint,n,theta_vec,y_i_vec)
mu = 0;
for ii = 1:n
    mu = mu +  exp(- dot(theta_vec, (training_data(ii,:) - testPoint).^2))*y_i_vec(ii) ;
end

end


