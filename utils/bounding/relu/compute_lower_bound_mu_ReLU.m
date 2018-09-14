function [lb,ub,phi_m] = compute_lower_bound_mu_ReLU(trainedSystem,phi_L,phi_U,sigma_b_2,sigma_w_2,training_data_angular,outputIdx,testPoint_angular,phi_m_star_prev...
    ,training_data,depth )
%
%compute_lower_bound_mu_ReLU computes the lower bound on the mu associated
%problem for the ReLU kernel
%

%getting the only ouput of interest
y_i_vec = trainedSystem(:,outputIdx);

n = size(training_data_angular,1);
m = size(training_data_angular,2) + 1;

a_i_sum = 0;
b_i_vec = zeros(n,1);


%looping over all the training points
for ii = 1:n
    y_i = y_i_vec(ii);
    psi_i = training_data_angular(ii,:);
    
    %getting initial betas
    [beta_l,beta_u,~,~] = compute_beta_interval(psi_i,phi_L,phi_U,sigma_b_2,sigma_w_2);

    
    
    k1 = sigma_b_2;
    
    a_L = zeros(depth,1);
    b_L = zeros(depth,1);
    
    a_U = zeros(depth,1);
    b_U = zeros(depth,1);
    K_diag = sigma_b_2 + (sigma_w_2/m);
    
    for hh = 1:depth
        %I compute the variation interval of gamma
        [gamma_l,gamma_u] = compute_gamma_interval(beta_l,beta_u);
        %Computing lower bound on current layer kernel
        
        k2 = (sigma_w_2/(2*pi) ) * K_diag;
        
        
        %looping over the number of layers. The approximation is done
        %iteratively on the layers
        if hh < depth
            %getting lower approximation
            [a_L(hh),b_L(hh),a_lower,b_lower,a_upper,b_upper] = compute_linear_under_bound_ReLU(1,gamma_l,gamma_u,beta_l,beta_u,k1,k2);
            %getting upper approximation
            [a_U(hh),b_U(hh)] = compute_linear_under_bound_ReLU(-1,gamma_l,gamma_u,beta_l,beta_u,k1,k2,a_lower,b_lower,a_upper,b_upper);
            a_U(hh) = - a_U(hh);
            b_U(hh) = - b_U(hh);

            %scaling by current layer a-priori variance
            K_diag =  sigma_b_2 + (sigma_w_2/2)* K_diag;
            a_L(hh) = a_L(hh)/K_diag;
            a_U(hh) = a_U(hh)/K_diag;
            b_L(hh) = b_L(hh)/K_diag;
            b_U(hh) = b_U(hh)/K_diag;
            %computing minimal and maximal value of current layer kernel
            %function
            [beta_l,beta_u] = min_max_of_kernel_function(gamma_l,gamma_u,k1,k2,K_diag);
        else
            %if this is the final layer, I just need the lower bound
            [a_L(hh),b_L(hh)] = compute_linear_under_bound_ReLU(y_i,gamma_l,gamma_u,beta_l,beta_u,k1,k2);
        end
        
        
    end
    %getting overall linear lower bound coefficients by iterative unrolling
    %of each layer coefficients
    a_i = a_L(depth);
    b_i = b_L(depth);
    for hh = (depth-1):-1:1
        if b_i >=0
            a_i = a_i + b_i*a_L(hh);
            b_i = b_i * b_L(hh);
        else
            a_i = a_i + b_i*a_U(hh);
            b_i = b_i * b_U(hh);
        end
    end
    a_i = a_i + b_i * sigma_b_2/(sigma_b_2 + sigma_w_2/m);
    b_i = b_i * (sigma_w_2/m)/(sigma_b_2 + sigma_w_2/m);
    
    
    
    
    
    b_i_vec(ii) = b_i;
    a_i_sum = a_i_sum + a_i;
end


bigSumVector = zeros(1,m);
for ii = 1:size(training_data_angular,1)
    bigSumVector = bigSumVector + b_i_vec(ii)*training_data(ii,:);
end
normOfBigSumVector = norm(bigSumVector);
bigSumVector = bigSumVector/normOfBigSumVector;
bigSumVector = get_spherical_coordinates(bigSumVector,[]);

[lb,phi_m,~,~] = min_max_of_dot_product_hypersphere(bigSumVector,phi_L,phi_U);
lb = a_i_sum + normOfBigSumVector*lb;
ub = get_actual_prediction(phi_m,training_data,sigma_b_2,sigma_w_2,y_i_vec,depth);

%checking if current result is actually better than the previous one
%obtained in the same region (if any). This can be tricky because of numerical cancellation erorr,
%when the region is very small. 
if ~isempty(phi_m_star_prev)
    if all(phi_m_star_prev>=phi_L) && (all(phi_m_star_prev<=phi_U))
        ub_prev = get_actual_prediction(phi_m_star_prev,training_data,sigma_b_2,sigma_w_2,y_i_vec,depth);
        [ub,ind] = min([ub,ub_prev]);
        phi_ms = [phi_m;phi_m_star_prev];
        phi_m = phi_ms(ind,:);
    end
end



end

function [mu,dotProdSum] = get_actual_prediction(testPoint,training_data,sigma_b_2,sigma_w_2,y_i_vec,depth)
%getting GP prediction on the input point testPoint
mu = 0;
dotProdSum = 0;
testPoint = get_cartesian_coordinates(testPoint,[]);
for ii = 1:size(training_data,1)
    [K,~] = getKernel(training_data(ii,:), testPoint   ,'N-ReLU',depth,sigma_w_2,sigma_b_2,[],[]);
    mu = mu +  K *y_i_vec(ii) ;
end

end