function [lb,ub] = compute_lower_bound_mu_sqe(trainedSystem,x_L,x_U,theta_vec,sigma_prior,training_data,outputIdx,debug,testPoint)

if nargin < 9
    testPoint = [];
    if nargin < 8
        debug = false;
    end
end

analytic_bound = false;

%getting only ouput of interest
%one_hot_labels = one_hot_encoding(training_labels);
%y_i_vec = one_hot_labels(:,outputIdx);
y_i_vec = trainedSystem(:,outputIdx);
%scaling output with a priori variance
y_i_vec = y_i_vec*sigma_prior;

n = size(training_data,1);

%debugging mu formula
%mu = get_actual_prediction(training_data,testPoint,n,theta_vec,y_i_vec);

m = size(training_data,2);
a_i_sum = 0;
theta_f_sum = 0;
b_i_vec = zeros(n,1);
for ii = 1:n
    y_i = y_i_vec(ii);
    x_i = training_data(ii,:);
    [z_i_L, z_i_U] = compute_z_interval(x_i,x_L,x_U,theta_vec);
    [a_i,b_i] = compute_linear_under_approx(y_i,z_i_L,z_i_U);
    b_i_vec(ii) = b_i;
    %b_i
    %plot_debug(z_i_L,z_i_U,a_i,b_i,y_i)
    %drawnow
    %pause(0.5)
    
    if analytic_bound
        inner_j_sum = 0;
        for jj = 1:m
            if b_i >= 0
                if (x_i(jj) <= x_U(jj)) && (x_i(jj) >= x_L(jj))
                    fstar_j_i = 0;
                else
                    fstar_j_i =  min( b_i *(x_i(jj) -  x_L(jj)  )^2 , b_i *(x_i(jj) -  x_U(jj)  )^2   );
                end
            else
                fstar_j_i =  min( b_i *(x_i(jj) -  x_L(jj)  )^2 , b_i *(x_i(jj) -  x_U(jj)  )^2   );
            end
            
            inner_j_sum = inner_j_sum + theta_vec(jj)*fstar_j_i;
            
        end
        theta_f_sum = theta_f_sum + inner_j_sum;
        
    end
    
    a_i_sum = a_i_sum + a_i;
end
if analytic_bound
    lb = a_i_sum + theta_f_sum;
end

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
%opts = optimoptions('quadprog');
%opts.MaxIterations = 200;
%opts.StepTolerance = 1.0e-03;
%opts.OptimalityTolerance = 1.0e-03;
%opts.ConstraintTolerance = 1.0e-03;
%if debug
%    opts.Display = 'iter';
%else
%    opts.Display = 'iter';
%end
%opts.Algorithm = 'trust-region-reflective';
%[x_star, f_val,exit_val] = quadprog(H,f,[],[],[],[],x_L',x_U',[],opts);
%f_val = f_val + a_i_sum + C;
%if exit_val ~= 1
%    disp('oh oh...')
%end
[x_star, f_val] = separate_quadprog(H,f,x_L,x_U);
f_val = f_val + a_i_sum + C;
%if sum(b_i_vec) >= 0
%    disp('oh oh...')
%end

ub = get_actual_prediction(training_data,x_star',n,theta_vec,y_i_vec);
lb = f_val;
% 
% if debug
%     %debugging
%     pix = linspace(x_L(103),x_U(103),10);
%     out_gp_debug = zeros(size(pix));
%     out_quadratic_debug = zeros(size(pix));
%     out_fu_under_debug = zeros(size(pix));
% 
%     for ii = 1:length(pix)
%         x_star(103) = pix(ii);
%         out_gp_debug(ii) = get_actual_prediction(training_data,x_star',n,theta_vec,y_i_vec);
%         out_quadratic_debug(ii) = 0.5*x_star'*H*x_star + f'* x_star + a_i_sum + C;
%         out_fu_under_debug(ii) = f_under(a_i_sum,b_i_vec,theta_vec,training_data,x_star');
%     end
%     figure
%     hold on
%     plot(pix,out_gp_debug)
%     plot(pix,lb*ones(size(pix)))
%     plot(pix,out_quadratic_debug)
%     plot(pix,f_val*ones(size(pix)))
%     plot(pix,out_fu_under_debug)
% 
%     legend({'gp','lb','quadratic','f_val','f_under'})
% end


end


function [f,g,H] = f_under(a_i_sum,b_i_vec,theta_vec,training_data,x)
inner_vec = zeros(size(training_data,1),1);
for ii = 1:size(training_data,1)
   inner_vec(ii) =  dot(theta_vec, (training_data(ii,:) - x).^2  );
end
f = a_i_sum + dot(b_i_vec,inner_vec);

grad = zeros(length(x));
for jj = 1:length(x)
    grad(jj) = -2*theta_vec(jj)* dot(b_i_vec, training_data(:,jj) - x(jj)); 
end

g = zeros(length(x));
for jj = 1:length(x)
    g(jj) = -2*theta_vec(jj)* dot(b_i_vec, training_data(:,jj) - x(jj)); 
end

H = speye(length(x));

for jj = 1:length(x)
    H(jj,jj) = 2*theta_vec(jj)*sum(b_i_vec);
end


end


function mu = get_actual_prediction(training_data,testPoint,n,theta_vec,y_i_vec)
mu = 0;
for ii = 1:n
    mu = mu +  exp(- dot(theta_vec, (training_data(ii,:) - testPoint).^2))*y_i_vec(ii) ;
end

end


