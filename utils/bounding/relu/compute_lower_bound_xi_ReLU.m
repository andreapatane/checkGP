function [lb_on_xi,ub_on_xi,x_opt_T] = compute_lower_bound_xi_ReLU(S,phi_L,phi_U,...
    sigma_b_2,sigma_w_2,training_data_angular,training_data, ...
    testPoint_angular,testPoint,r_star_star,x_opt_T_prev,depth)
%
%Function for xi bounding for ReLU kernel
%

m = length(phi_L) + 1;
n = size(training_data_angular,1);
%defining f vector and hessian for quadratic program
f = [ zeros(m,1) ; 2 ; -2 * S * r_star_star ];
Q = [ zeros(m+1,m+1) , zeros(m+1,n); zeros(n,m+1), S ];

A = zeros(2+2*n,m+1+n);
b = zeros(2+2*n,1);

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
    [a_ii_low,b_ii_low,a_ii_up,b_ii_up,kernel_low_training(ii),kernel_up_training(ii)] =...
        compute_lower_upper_bounds_on_deep_ReLU(training_data_angular(ii,:),phi_L,phi_U,sigma_b_2,sigma_w_2,m,depth);

    f_ii_minus = b_ii_low*training_data(ii,:);
    A(1 + 2*ii ,:) = [f_ii_minus,0,zeros(1,n)];
    A(1 + 2*ii ,m+1+ii) = -1;
    b(1 + 2*ii) = - a_ii_low;

    f_ii_plus = - b_ii_up*training_data(ii,:);
    A(1 + 2*ii + 1 ,:) = [f_ii_plus,0,zeros(1,n)];
    A(1 + 2*ii + 1 , m+1+ii) = 1;
    b(1 + 2*ii + 1) =  a_ii_up;
end

[x_L,x_U] = cartesian_encolising_hypercube(phi_L,phi_U);


x_L = [x_L;kernel_low_test;kernel_low_training];
x_U = [x_U;kernel_up_test;kernel_up_training];



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



priori_var = sigma_b_2 + (sigma_w_2/m);

for hh = 1:depth
    priori_var =  sigma_b_2 + (sigma_w_2/2)* priori_var;
end
aux = (2*priori_var  - r_star_star' * S * r_star_star);
ub_on_xi = aux - lb;
lb_on_xi = ub;

end



function sigma = get_actual_prediction(training_data,x_opt_T,testPoint,sigma_b_2,sigma_w_2,S,depth)

[Kstar,Kstarstar] = compute_kernel_on_test([x_opt_T;testPoint],training_data,'N-ReLU',depth,sigma_w_2,sigma_b_2,{});

sigma_test = Kstarstar - Kstar * S * Kstar';
sigma = sigma_test(1,1) + sigma_test(2,2) - 2 * sigma_test(1,2);
end


