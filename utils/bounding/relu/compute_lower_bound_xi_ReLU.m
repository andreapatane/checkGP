function [lb_on_xi,ub_on_xi,x_opt_T] = compute_lower_bound_xi_ReLU(S,phi_L,phi_U,...
    sigma_b_2,sigma_w_2,training_data_angular,training_data, ...
    testPoint_angular,testPoint,r_star_star,x_opt_T_prev,depth)


m = length(phi_L) + 1;
n = size(training_data_angular,1);
f = [ zeros(m,1) ; 2 ; -2 * S * r_star_star ];
Q = [ zeros(m+1,m+1) , zeros(m+1,n); zeros(n,m+1), S ];

A = zeros(2+2*n,m+1+n);
b = zeros(2+2*n,1);

%k1 = sigma_b_2;
%k2 = sigma_w_2 * (sigma_b_2 + sigma_w_2/m) ;
%fix
%k2 = k2/(2*pi);
%
%kernel_fun = @(t)(k1 + k2 * (sin(t) + (pi-t).*cos(t) )   );


% [gamma_l_star,gamma_u_star,~,~,beta_l_star,beta_u_star] = compute_gamma_interval(testPoint_angular,phi_L,phi_U,sigma_b_2,sigma_w_2);
% [a_star_low,b_star_low] = compute_linear_under_bound_ReLU(1,gamma_l_star,gamma_u_star,beta_l_star,beta_u_star,sigma_b_2,sigma_w_2,m);
% [a_star_up,b_star_up] = compute_linear_under_bound_ReLU(-1,gamma_l_star,gamma_u_star,beta_l_star,beta_u_star,sigma_b_2,sigma_w_2,m);
% 
% 
% kernel_low_test = kernel_fun(gamma_u_star);
% kernel_up_test = kernel_fun(gamma_l_star);
[a_star_low,b_star_low,a_star_up,b_star_up,kernel_low_test,kernel_up_test] =...
    compute_lower_upper_bounds_on_deep_ReLU(testPoint_angular,phi_L,phi_U,sigma_b_2,sigma_w_2,m,depth);

f_star_minus = b_star_low*testPoint;
A(1,:) = [f_star_minus,-1,zeros(1,n)];
b(1) = - a_star_low;

f_star_plus = - b_star_up * testPoint;
A(2,:) = [f_star_plus,1,zeros(1,n)];
b(2) =  a_star_up;



kernel_low_training = zeros(n,1);
kernel_up_training = zeros(n,1);
for ii = 1:n
%     [gamma_l_ii,gamma_u_ii,~,~,beta_l_ii,beta_u_ii,~,~] = ...
%         compute_gamma_interval(training_data_angular(ii,:),phi_L,phi_U,sigma_b_2,sigma_w_2);
%     kernel_low_training(ii) = kernel_fun(gamma_u_ii);
%     kernel_up_training(ii) = kernel_fun(gamma_l_ii);
%     [a_ii,b_ii] = compute_linear_under_bound_ReLU(1,gamma_l_ii,gamma_u_ii,beta_l_ii,beta_u_ii,sigma_b_2,sigma_w_2,m);
    [a_ii_low,b_ii_low,a_ii_up,b_ii_up,kernel_low_training(ii),kernel_up_training(ii)] =...
        compute_lower_upper_bounds_on_deep_ReLU(training_data_angular(ii,:),phi_L,phi_U,sigma_b_2,sigma_w_2,m,depth);

    f_ii_minus = b_ii_low*training_data(ii,:);
    A(1 + 2*ii ,:) = [f_ii_minus,0,zeros(1,n)];
    A(1 + 2*ii ,m+1+ii) = -1;
    b(1 + 2*ii) = - a_ii_low;

%     [a_ii,b_ii] = compute_linear_under_bound_ReLU(-1,gamma_l_ii,gamma_u_ii,beta_l_ii,beta_u_ii,sigma_b_2,sigma_w_2,m);
    f_ii_plus = - b_ii_up*training_data(ii,:);
    A(1 + 2*ii + 1 ,:) = [f_ii_plus,0,zeros(1,n)];
    A(1 + 2*ii + 1 , m+1+ii) = 1;
    b(1 + 2*ii + 1) =  a_ii_up;
end

[x_L,x_U] = cartesian_encolising_hypercube(phi_L,phi_U);


x_L = [x_L;kernel_low_test;kernel_low_training];
x_U = [x_U;kernel_up_test;kernel_up_training];

%x_L = [x_L;-inf;-inf(n,1)];
%x_U = [x_U;inf;inf(n,1)];

opts = optimoptions('quadprog');
opts.Display = 'off';

[x_opt,lb,exitFlag] = quadprog(2*Q,f,A,b,[],[],x_L,x_U,[],opts);
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

ub = get_actual_prediction(training_data,x_opt_T,testPoint,sigma_b_2,sigma_w_2,S,depth);
if ~isempty(x_opt_T_prev)
    if all(get_spherical_coordinates(x_opt_T_prev,[]) >= phi_L) && all(get_spherical_coordinates(x_opt_T_prev,[]) <= phi_U)
        ub1 = get_actual_prediction(training_data,x_opt_T_prev,testPoint,sigma_b_2,sigma_w_2,S,depth);
        [ub,idx] = min([ub,ub1]);
        if idx == 2
            x_opt_T = x_opt_T_prev;
        end
    end
end
%tic
%opts = optimoptions('quadprog');
%[x_opt,c_lb,exitFlag] = quadprog(2*Q,f,[],[],[],[],x_L,x_U,x0,opts);
%toc


priori_var = sigma_b_2 + (sigma_w_2/m);

for hh = 1:depth
    priori_var =  sigma_b_2 + (sigma_w_2/2)* priori_var;
end
%priori_var = sigma_b_2 + 0.5*sigma_w_2 * (sigma_b_2 + sigma_w_2/m );
aux = (2*priori_var  - r_star_star' * S * r_star_star);
ub_on_xi = aux - lb;
lb_on_xi = ub;

% 
% %%%DEBUGGING ON THREE ANGLES
% %code below assume the first three angles are being changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% % % disc_xi = zeros(size(x_debug));
% % % for ii = 1:length(t1)
% % %     for jj = 1:length(t2)
% % %         for kk = 1:length(t3)
% % %             phi_current = phi_L;
% % %             phi_current(1) = t1(ii);
% % %             phi_current(2) = t2(jj);
% % %             phi_current(3) = t3(jj);
% % %             x_debug{ii,jj,kk} = get_cartesian_coordinates(phi_current,[]);
% % %             disc_xi(ii,jj,kk) = get_actual_prediction(training_data,x_debug{ii,jj,kk},testPoint,sigma_b_2,sigma_w_2,S);
% % %         end
% % %     end
% % % end
% % % disp(lb_on_xi)
% % % disp(max(max(max(disc_xi))))
% % % disp(ub_on_xi)


end



function sigma = get_actual_prediction(training_data,x_opt_T,testPoint,sigma_b_2,sigma_w_2,S,depth)

gp_training_opts.kernel = 'ReLU';
gp_training_opts.kernel_params.sigma_b_2 = sigma_b_2;
gp_training_opts.kernel_params.sigma_w_2 = sigma_w_2;
gp_training_opts.kernel_params.depth = depth;
[Kstar,Kstarstar] = compute_kernel_on_test([x_opt_T;testPoint],training_data,gp_training_opts);
sigma_test = Kstarstar - Kstar * S * Kstar';
sigma = sigma_test(1,1) + sigma_test(2,2) - 2 * sigma_test(1,2);
end


