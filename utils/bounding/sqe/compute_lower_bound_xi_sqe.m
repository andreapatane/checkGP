function [lb_on_xi,ub_on_xi] = compute_lower_bound_xi_sqe(R_inv,x_L,x_U,theta_vec,sigma_prior,training_data,...
                                                outputIdx,debug,testPoint,r_star,R,output_mode)



R_inv = sigma_prior * R_inv;
R  = R / sigma_prior;

n = size(training_data,1);
m = size(training_data,2);

a_il_sum = 0;
B_il =  zeros(n,n);
b = zeros(n,1);
a_i_sum = 0;

% computing a_i and b_i
z_i_L_vec = zeros(1,n);
z_i_U_vec = zeros(1,n);

%tic
%r_star = compute_kernel_on_test(training_data,testPoint,'sqe',0,1,theta_vec,{});
F = -2 * r_star'*R_inv;
for ii = 1:n
    [z_i_L_vec(ii),z_i_U_vec(ii)] = compute_z_interval(training_data(ii,:),x_L,x_U,theta_vec);
    %F_i = -2 * dot( R_inv(:,ii), r_star  );
    [a_i,b_i] = compute_linear_under_approx(F(ii),z_i_L_vec(ii),z_i_U_vec(ii));
    b(ii) = b_i;
    a_i_sum = a_i_sum + a_i;
end
%t1 = toc;

% computing a_il and b_il
for ii = 1:n
    for ll = 1:ii
        z_il_L = z_i_L_vec(ii) + z_i_L_vec(ll);
        z_il_U = z_i_U_vec(ii) + z_i_U_vec(ll);
        %tic
        % this is the function call that I should put here
        %[a_il,B_il(ii,ll)] = compute_linear_under_approx(R_inv(ii,ll),z_il_L,z_il_U);
        %however I inline it in the following:
        if R_inv(ii,ll) >= 0
            z_il_M = 0.5*(z_il_L + z_il_U);
            a_il = (1 + z_il_M)*R_inv(ii,ll)*exp(-z_il_M);
            B_il(ii,ll) = -R_inv(ii,ll)*exp(-z_il_M);
            
        else
            a_il = R_inv(ii,ll) * ( exp(-z_il_L)  - z_il_L * (  exp(-z_il_L) -  exp(-z_il_U) )/( z_il_L - z_il_U   )  );
            B_il(ii,ll) = R_inv(ii,ll) * (  exp(-z_il_L) -  exp(-z_il_U) )/( z_il_L - z_il_U   );
            %if z_i_L == z_i_U
            %    disp('');
            %end
        end
        %I know it's ugly but it gives a significant speed-up for each bnb iteration...
        %t1 = t1 + toc;
        %taking into account symmetry...
        if ll < ii
            a_il = 2*a_il;
            B_il(ll,ii) = B_il(ii,ll);
        end
        a_il_sum = a_il_sum + a_il;
    end
end
%t2 = toc;
% computing a_* and b_*

%tic
[z_star_L , z_star_U ] = compute_z_interval(testPoint,x_L,x_U,theta_vec);
[a_star, b_star] = compute_linear_under_approx(2,z_star_L,z_star_U);
%t3 = toc;
% computing H, f and C for quadratic program
%tic
B_sum = sum(B_il)';
H = zeros(m,1);
f = zeros(m,1);
C = 0;
for jj = 1:m
    H(jj) = theta_vec(jj) * (2*sum(B_sum) + sum(b) + b_star );
    f(jj) = -2*theta_vec(jj) * (  dot(2*B_sum + b  , training_data(:,jj)) + b_star * testPoint(jj));
    C = C + theta_vec(jj) * (  dot( 2*B_sum + b, training_data(:,jj).^2   ) + b_star * testPoint(jj)^2);
end

C = C + a_star + a_il_sum + a_i_sum;
%t4 = toc;
%solving the quadratic program
%tic
[x_star, f_val] = separate_quadprog(2*H',f,x_L,x_U);
f_val = f_val + C;
lb = f_val;

%getting upperbound as well
ub = get_actual_prediction(training_data,x_star',testPoint,theta_vec,R,output_mode);
if ub < 0
    disp('....')
%    stop
end
%scaling stuff with prior variance
%lb = sigma_prior*lb;
%ub = sigma_prior*ub;
aux = (2 - r_star' * R_inv * r_star);
ub_on_xi = sigma_prior * (aux - lb);
lb_on_xi = sigma_prior * ub;
%t5 = toc;
%disp('Time per phase:')
%disp(t1)
%disp(t2)
%disp(t3)
%disp(t4)
%disp(t5)

% % Old implementation, does not take into account covariance...
% % t1 = 0;
% % t2 = 0;
% % t3 = 0;
% % 
% % for ii = 1:n
% %     
% %     for ll = 1:ii
% %         %H = ;
% %         %f = - theta_vec'.*(training_data(ii,:) + training_data(ll,:));
% %         tic
% %         f = - (training_data(ii,:) + training_data(ll,:));
% %         %[~, z_il_L,~, z_il_U] = separate_quadprog(H,f,x_L,x_U);
% %         [z_il_L, z_il_U] = simult_separate_quadprog_Hg0(theta_vec,f,x_L,x_U);
% %         t1 = t1 + toc;
% %         tic
% %         xxx =  dot(theta_vec,training_data(ii,:).^2 + training_data(ll,:).^2);
% %         z_il_L = 2*z_il_L + xxx;
% %         z_il_U = 2*z_il_U + xxx;
% %         t2 = t2 + toc;
% %         %tic 
% %         %[~ , z_il_U] = separate_quadprog(-H,-f,x_L,x_U);
% %         %t2 = t2 + toc;
% %         %z_il_U = -2*z_il_U + dot(theta_vec,training_data(ii,:).^2 + training_data(ll,:).^2);
% %         
% %         tic
% %         [a_il,B_il(ii,ll)] = compute_linear_under_approx(R_inv(ii,ll),z_il_L,z_il_U);
% %         t3 = t3 + toc;
% %         %plot_debug(z_il_L,z_il_U,a_il,b_il,R_inv(ii,ll))
% %         a_il_sum = a_il_sum + a_il;
% %     end
% %     
% %     
% % end
% % 
% % B_il = B_il + B_il' - diag(diag(B_il));
% % B_sum = sum(B_il,2);
% % 
% % C = 0;
% % for ii = 1:n
% %     C = C + B_sum(ii) * dot(theta_vec,training_data(ii,:).^2);
% % end
% % C = 2*C;
% % H = 2*sum(B_sum) * theta_vec;
% % 
% % f = zeros(m,1);
% % for  jj = 1:m
% %     f(jj) = -4*theta_vec(jj) * dot(B_sum,training_data(:,jj));
% % end
% % 
% % 
% % [x_star, f_val] = separate_quadprog(H',f,x_L,x_U);
% % f_val = f_val + a_il_sum + C;
% % lb = f_val;
% % 
% % ub = get_actual_prediction(training_data,x_star',theta_vec,R_inv);
% % 
% % lb = sigma_prior*lb;
% % ub = sigma_prior*ub;
end




function sigma = get_actual_prediction(training_data,x_star,testPoint,theta_vec,R,output_mode)

gp_training_opts.kernel = 'sqe';
gp_training_opts.mle = false;

gp_training_opts.kernel_params.sigma = 1;
gp_training_opts.kernel_params.theta_vec = theta_vec;


[Kstar,Kstarstar] = compute_kernel_on_test([x_star;testPoint],training_data,gp_training_opts);
[~,sigma_test] = make_predictions_on_test(Kstar,[],Kstarstar,R,output_mode);
sigma = sigma_test(1,1) + sigma_test(2,2) - 2 * sigma_test(1,2);

end



