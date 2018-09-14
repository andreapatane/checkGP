function [lb_on_sigma,ub_on_sigma,x_opt_T] = compute_lower_bound_sigma_ReLU(S,phi_L,phi_U,...
    sigma_b_2,sigma_w_2,training_data_angular,training_data,x_opt_T_prev,depth)
%
%Lower bound for Variance analysis on the ReLU kernel
%
m = length(phi_L) + 1;
n = size(training_data_angular,1);

%defining matrix and vector for quadratic programming
%hessian
Q = [ zeros(m,m) , zeros(m,n); zeros(n,m), S ];
%first degree weigth vector
f = zeros(size(Q,1),1);

%Ax<=b matrix and vector defined below:
A = zeros(2*n,m+n);
b = zeros(2*n,1);


kernel_low_training = zeros(n,1);
kernel_up_training = zeros(n,1);
for ii = 1:n
    %allocating lower and upper bounds on ReLU kernel
    [a_ii_low,b_ii_low,a_ii_up,b_ii_up,kernel_low_training(ii),kernel_up_training(ii)] =...
        compute_lower_upper_bounds_on_deep_ReLU(training_data_angular(ii,:),phi_L,phi_U,sigma_b_2,sigma_w_2,m,depth);
    
    %updating constraints
    f_ii_minus = b_ii_low*training_data(ii,:);
    A(2*ii - 1 ,:) = [f_ii_minus,zeros(1,n)];
    A(2*ii - 1 ,m+ii) = -1;
    b(2*ii - 1) = - a_ii_low;
    
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
tic
[x_opt,lb,exitFlag] = quadprog(2*Q,f,A,b,[],[],x_L,x_U,[],opts);
toc
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


%computing upper bound using actual prediction below:

%first I project the solution into the actual region under analysis:
x_opt = x_opt(1:m);
x_opt = x_opt/norm(x_opt);
phi_opt = get_spherical_coordinates(x_opt',[]);
phi_opt = min(max([phi_opt;phi_L]),phi_U);

%I convert the projected solution into cartesian coordinates
x_opt_T = get_cartesian_coordinates(phi_opt,[]);
%and compute its true prediction
ub = get_actual_prediction(training_data,x_opt_T,sigma_b_2,sigma_w_2,S,depth);

%check if x_opt_T_prev is better than the new point
if ~isempty(x_opt_T_prev)
    if all(get_spherical_coordinates(x_opt_T_prev,[]) >= phi_L) && all(get_spherical_coordinates(x_opt_T_prev,[]) <= phi_U)
        ub1 = get_actual_prediction(training_data,x_opt_T_prev,sigma_b_2,sigma_w_2,S);
        [ub,idx] = min([ub,ub1]);
        if idx == 2
            x_opt_T = x_opt_T_prev;
        end
    end
end


%some arithmetic computation to take into account for the a-priori variance
priori_var = sigma_b_2 + (sigma_w_2/m);

for hh = 1:depth
    priori_var =  sigma_b_2 + (sigma_w_2/2)* priori_var;
end

ub_on_sigma = priori_var - lb;
lb_on_sigma = ub;





end



function sigma = get_actual_prediction(training_data,x_opt_T,sigma_b_2,sigma_w_2,S,depth)

[Kstar,Kstarstar] = compute_kernel_on_test(x_opt_T,training_data,'N-ReLU',depth,sigma_w_2,sigma_b_2,{});

sigma = Kstarstar - Kstar * S * Kstar';
end

