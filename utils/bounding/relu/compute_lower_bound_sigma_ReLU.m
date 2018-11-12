function [lb_on_sigma,ub_on_sigma,x_opt_T] = compute_lower_bound_sigma_ReLU(S,phi_L,phi_U,...
    sigma_b_2,sigma_w_2,training_data_angular,training_data,x_opt_T_prev,depth)

m = length(phi_L) + 1;
n = size(training_data_angular,1);

Q = [ zeros(m,m) , zeros(m,n); zeros(n,m), S ];
f = zeros(size(Q,1),1);

A = zeros(2*n,m+n);
b = zeros(2*n,1);


%k1 = sigma_b_2;
%k2 = sigma_w_2 * (sigma_b_2 + sigma_w_2/m) ;
%kernel_fun = @(t)(k1 + k2 * (sin(t) + (pi-t).*cos(t) )   );


kernel_low_training = zeros(n,1);
kernel_up_training = zeros(n,1);
for ii = 1:n
    %[gamma_l_ii,gamma_u_ii,~,~,beta_l_ii,beta_u_ii,~,~] = ...
    %    compute_gamma_interval(training_data_angular(ii,:),phi_L,phi_U,sigma_b_2,sigma_w_2);
    %kernel_low_training(ii) = kernel_fun(gamma_u_ii);
    %kernel_up_training(ii) = kernel_fun(gamma_l_ii);
    %[a_ii,b_ii] = compute_linear_under_bound_ReLU(1,gamma_l_ii,gamma_u_ii,beta_l_ii,beta_u_ii,sigma_b_2,sigma_w_2,m);
    [a_ii_low,b_ii_low,a_ii_up,b_ii_up,kernel_low_training(ii),kernel_up_training(ii)] =...
        compute_lower_upper_bounds_on_deep_ReLU(training_data_angular(ii,:),phi_L,phi_U,sigma_b_2,sigma_w_2,m,depth);

    f_ii_minus = b_ii_low*training_data(ii,:);
    A(2*ii - 1 ,:) = [f_ii_minus,zeros(1,n)];
    A(2*ii - 1 ,m+ii) = -1;
    b(2*ii - 1) = - a_ii_low;

    %[a_ii,b_ii] = compute_linear_under_bound_ReLU(-1,gamma_l_ii,gamma_u_ii,beta_l_ii,beta_u_ii,sigma_b_2,sigma_w_2,m);
    f_ii_plus = - b_ii_up * training_data(ii,:);
    
    A(2*ii ,:) = [f_ii_plus,zeros(1,n)];
    A(2*ii , m+ii) = 1;
    b(2*ii ) =  a_ii_up;
end

[x_L_cart,x_U_cart] = cartesian_encolising_hypercube(phi_L,phi_U);


x_L = [x_L_cart;kernel_low_training];
x_U = [x_U_cart;kernel_up_training];


opts = optimoptions('quadprog');
opts.Display = 'off';
%tic
[x_opt,lb,exitFlag] = quadprog(2*Q,f,A,b,[],[],x_L,x_U,[],opts);
%toc
switch exitFlag
    case 0
        disp('number of iterations for quadprog exceeded')
    case -2
        error('problem unfeasible?! This should not happen really...')
    case -3
        error('Unbounded problem?! This should not happen really...')
    case 2
        disp('small step, but I was still slowly improving...')
    case -6
        error('Non-convex problem... Check for numerical cancellation errors')
    case -8
        error('Some random stuff related to step direction. Please do not occur!')
end

x_opt = x_opt(1:m);
x_opt = x_opt/norm(x_opt);

phi_opt = get_spherical_coordinates(x_opt',[]);
phi_opt = min(max([phi_opt;phi_L]),phi_U);
x_opt_T = get_cartesian_coordinates(phi_opt,[]);
%x_opt_T = x_opt;

ub = get_actual_prediction(training_data,x_opt_T,sigma_b_2,sigma_w_2,S,depth);
if ~isempty(x_opt_T_prev)
    if all(get_spherical_coordinates(x_opt_T_prev,[]) >= phi_L) && all(get_spherical_coordinates(x_opt_T_prev,[]) <= phi_U)
    %if all(x_opt_T_prev >= x_L_cart) && all(x_opt_T_prev <= x_U_cart)
        ub1 = get_actual_prediction(training_data,x_opt_T_prev,sigma_b_2,sigma_w_2,S);
        [ub,idx] = min([ub,ub1]);
        if idx == 2
            x_opt_T = x_opt_T_prev;
        end
    end
end

%priori_var = sigma_b_2 + pi * sigma_w_2 * (sigma_b_2 + sigma_w_2/m );
priori_var = sigma_b_2 + (sigma_w_2/m);

for hh = 1:depth
    priori_var =  sigma_b_2 + (sigma_w_2/2)* priori_var;
end

ub_on_sigma = priori_var - lb;
lb_on_sigma = ub;




% % % % %%%DEBUGGING ON THREE ANGLES
% % % % %code below assume the first three angles are being changed
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % %preparing discrete grid
% % % phi_L1 = phi_L(1);
% % % phi_U1 = phi_U(1);
% % % 
% % % phi_L2 = phi_L(2);
% % % phi_U2 = phi_U(2);
% % % 
% % % phi_L3 = phi_L(3);
% % % phi_U3 = phi_U(3);
% % % 
% % % t1 = linspace(phi_L1,phi_U1,20);
% % % t2 = linspace(phi_L2,phi_U2,20);
% % % t3 = linspace(phi_L3,phi_U3,20);
% % % 
% % % x_debug = cell(length(t1),length(t2),length(t3));
% % % disc_sigma = zeros(size(x_debug));
% % % for ii = 1:length(t1)
% % %     for jj = 1:length(t2)
% % %         for kk = 1:length(t3)
% % %             phi_current = phi_L;
% % %             phi_current(1) = t1(ii);
% % %             phi_current(2) = t2(jj);
% % %             phi_current(3) = t3(jj);
% % %             x_debug{ii,jj,kk} = get_cartesian_coordinates(phi_current,[]);
% % %             disc_sigma(ii,jj,kk) = get_actual_prediction(training_data,x_debug{ii,jj,kk},sigma_b_2,sigma_w_2,S);
% % %         end
% % %     end
% % % end
% % % disp(lb_on_sigma)
% % % disp(max(max(max(disc_sigma))))
% % % disp(ub_on_sigma)


end



function sigma = get_actual_prediction(training_data,x_opt_T,sigma_b_2,sigma_w_2,S,depth)

gp_training_opts.kernel = 'ReLU';
gp_training_opts.kernel_params.sigma_b_2 = sigma_b_2;
gp_training_opts.kernel_params.sigma_w_2 = sigma_w_2;
gp_training_opts.kernel_params.depth = depth;

[Kstar,Kstarstar] = compute_kernel_on_test(x_opt_T,training_data,gp_training_opts);

sigma = Kstarstar - Kstar * S * Kstar';
end

