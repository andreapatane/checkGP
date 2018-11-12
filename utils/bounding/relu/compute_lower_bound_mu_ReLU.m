function [lb,ub,phi_m] = compute_lower_bound_mu_ReLU(trainedSystem,phi_L,phi_U,sigma_b_2,sigma_w_2,training_data_angular,outputIdx,testPoint_angular,phi_m_star_prev...
    ,training_data,depth )



%getting only ouput of interest
%one_hot_labels = one_hot_encoding(training_labels);
%y_i_vec = one_hot_labels(:,outputIdx);
y_i_vec = trainedSystem(:,outputIdx);
%scaling output with a priori variance

n = size(training_data_angular,1);
m = size(training_data_angular,2) + 1;

a_i_sum = 0;
b_i_vec = zeros(n,1);
for ii = 1:n
    y_i = y_i_vec(ii);
    psi_i = training_data_angular(ii,:);
    
    %getting initial betas
    [beta_l,beta_u,dot_L,dot_U] = compute_beta_interval(psi_i,phi_L,phi_U,sigma_b_2,sigma_w_2);
    beta_l_init = beta_l;
    beta_u_init = beta_u;
    %k2 = sigma_b_2 + sigma_w_2/m;
    k1 = sigma_b_2;
    
    a_L = zeros(depth,1);
    b_L = zeros(depth,1);
    
    a_U = zeros(depth,1);
    b_U = zeros(depth,1);
    K_diag = sigma_b_2 + (sigma_w_2/m);
    kernelFun1 = @(beta_x,k1_p,k2_p) (k1_p+k2_p*(sin( acos(beta_x)) + (pi-acos(beta_x)).*cos(acos(beta_x)) ));
    kernelFun2 = @(beta_x,k12_p,k22_p,k1_p,k2_p,K_diag) (k12_p + k22_p * ( sin(acos( kernelFun1(beta_x,k1_p,k2_p)/K_diag ) +...
                  (pi- acos(kernelFun1(beta_x,k1_p,k2_p)/K_diag )).*cos(acos(kernelFun1(beta_x,k1_p,k2_p)/K_diag )) )  )  );
    for hh = 1:depth
        %I compute the variation interval of gamma
        [gamma_l,gamma_u] = compute_gamma_interval(beta_l,beta_u);
        %Computing lower bound on current layer kernel
        
        k2 = (sigma_w_2/(2*pi) ) * K_diag;
        
        if hh == 1
            k2_init = k2;
        end
        
        
        if hh < depth
            [a_L(hh),b_L(hh),a_lower,b_lower,a_upper,b_upper] = compute_linear_under_bound_ReLU(1,gamma_l,gamma_u,beta_l,beta_u,k1,k2);
            [a_U(hh),b_U(hh)] = compute_linear_under_bound_ReLU(-1,gamma_l,gamma_u,beta_l,beta_u,k1,k2,a_lower,b_lower,a_upper,b_upper);
            a_U(hh) = - a_U(hh);
            b_U(hh) = - b_U(hh);
            if 0
                t_beta = linspace(beta_l,beta_u,10);
                figure
                hold on
                kerCore = @(x) (sin(x) + (pi-x).*cos(x));
                kernelFun = @(beta_x,k1,k2) (k1 + k2*(kerCore(acos(beta_x))));
                plot(t_beta,a_L(hh)+b_L(hh).*t_beta);
                plot(t_beta,a_U(hh)+b_U(hh).*t_beta);
                plot(t_beta,kernelFun(t_beta,k1,k2));
                legend({'lower','upper','true'})
                hold off
            end
            K_diag =  sigma_b_2 + (sigma_w_2/2)* K_diag;
            a_L(hh) = a_L(hh)/K_diag;
            a_U(hh) = a_U(hh)/K_diag;
            b_L(hh) = b_L(hh)/K_diag;
            b_U(hh) = b_U(hh)/K_diag;
            [beta_l,beta_u] = min_max_of_kernel_function(gamma_l,gamma_u,k1,k2,K_diag);
            %[beta_l,beta_u]
            %disp('.')
        else
            [a_L(hh),b_L(hh)] = compute_linear_under_bound_ReLU(y_i,gamma_l,gamma_u,beta_l,beta_u,k1,k2);
            if 0
                close
                kerCore = @(x) (sin(x) + (pi-x).*cos(x));
                kernelFun = @(beta_x,k1,k2) (k1 + k2*(kerCore(acos(beta_x))));
                t_beta = linspace(beta_l,beta_u,10);
                figure
                hold on
                plot(t_beta,a_L(hh)+b_L(hh).*t_beta);
                plot(t_beta,kernelFun(t_beta,k1,k2));
                legend({'lower','true'})
                drawnow
                pause(1.0)
                %close
            end
        end
        
        
        %[a_i,b_i] = compute_linear_under_bound_ReLU(y_i,gamma_l,gamma_u,beta_l,beta_u,sigma_b_2,sigma_w_2,m);
    end
    
    a_i = a_L(depth);
    b_i = b_L(depth);
    for hh = (depth-1):-1:1
        if b_i >=0
            %disp('positivo...')
            a_i = a_i + b_i*a_L(hh);
            b_i = b_i * b_L(hh);
        else
            %disp('negativo...')
            a_i = a_i + b_i*a_U(hh);
            b_i = b_i * b_U(hh);
        end
    end
    if 0
        figure 
        hold on
        beta_init_x = linspace(beta_l_init,beta_u_init,10);
        plot(beta_init_x,a_i + b_i*beta_init_x );
        plot(beta_init_x,kernelFun2(beta_init_x,k1,k2,k1,k2_init,K_diag))
        legend({'lower','actual'})
        hold off
    end
    
    a_i = a_i + b_i * sigma_b_2/(sigma_b_2 + sigma_w_2/m);
    b_i = b_i * (sigma_w_2/m)/(sigma_b_2 + sigma_w_2/m);
    
    if 0
        figure 
        hold on
        dot_x = linspace(dot_L,dot_U,10);
        plot(dot_x,a_i + b_i*dot_x );
        plot(dot_x,kernelFun2(sigma_b_2/(sigma_b_2 + sigma_w_2/m) + dot_x * (sigma_w_2/m)/(sigma_b_2 + sigma_w_2/m)  ,k1,k2,k1,k2_init,K_diag))
        legend({'lower','actual'})
        hold off
    end
    
    
    b_i_vec(ii) = b_i;
    if 0
        %disp(b_i)
        %pixels2Modify = find(phi_U > phi_L);
        %this assume the only angle being modified is the first one
        tt = linspace(phi_L(1),phi_U(1),10);
        mu = zeros(size(tt));
        mu_low = zeros(size(tt));
        for jj = 1:length(tt)
            phi_test = phi_L;
            phi_test(1) = tt(jj);
            mu(jj) = get_actual_prediction(phi_test,training_data(ii,:),sigma_b_2,sigma_w_2,y_i_vec(ii),depth);
            dotProd = dot(get_cartesian_coordinates(phi_test,[]),get_cartesian_coordinates(training_data_angular(ii,:),[]));
            
            mu_low(jj) = a_i + b_i * dotProd;
        end
        %assert(all(mu >= mu_low))
        %disp([min(mu - mu_low),max(mu - mu_low)])
        if any(mu < mu_low)
            disp('accipicchia!')
            %disp(y_i)
        else
            disp('questo ok...')
            %disp(y_i)
        end
        figure
        hold on
        plot(tt,mu);
        plot(tt,mu_low);
        legend({'actual','lower'})
        drawnow
        pause(1.0)
         close
    end
     a_i_sum = a_i_sum + a_i;
end

if 0
    tt = linspace(phi_L(1),phi_U(1),10);
    mu = zeros(size(tt));
    mu_low = zeros(size(tt));
    for jj = 1:length(tt)
        phi_test = phi_L;
        phi_test(end) = tt(jj);
        mu(jj) = get_actual_prediction(phi_test,training_data,sigma_b_2,sigma_w_2,y_i_vec,depth);
        dotProd = 0;
        sumVec_debug = b_i_vec(1) * training_data(1,:);
        for ii = 2:size(training_data_angular,1)
            sumVec_debug = sumVec_debug + b_i_vec(ii) * training_data(ii,:);
        end
        dotProd = dot(get_cartesian_coordinates(phi_test,[]), sumVec_debug);
        %for ii = 1:size(training_data_angular,1)
            %dotProd = dotProd + b_i_vec(ii) * dot(get_cartesian_coordinates(phi_test,[]),get_cartesian_coordinates(training_data_angular(ii,:),[]));
        %    dotProd = dotProd + dot(get_cartesian_coordinates(phi_test,[]), b_i_vec(ii) * training_data(ii,:));
        %end
        mu_low(jj) = a_i_sum +  dotProd;
    end
    figure
    hold on
    plot(tt,mu);
    plot(tt,mu_low);
    
    drawnow
    pause(2.0)
    %close
    
end

bigSumVector = zeros(1,m);
for ii = 1:size(training_data_angular,1)
    bigSumVector = bigSumVector + b_i_vec(ii)*training_data(ii,:);
end
normOfBigSumVector = norm(bigSumVector);
bigSumVector = bigSumVector/normOfBigSumVector;
bigSumVector = get_spherical_coordinates(bigSumVector);

[lb,phi_m,~,~] = min_max_of_dot_product_hypersphere(bigSumVector,phi_L,phi_U);
lb = a_i_sum + normOfBigSumVector*lb;
ub = get_actual_prediction(phi_m,training_data,sigma_b_2,sigma_w_2,y_i_vec,depth);

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

gp_training_opts.kernel = 'ReLU';
gp_training_opts.kernel_params.sigma_b_2 = sigma_b_2;
gp_training_opts.kernel_params.sigma_w_2 = sigma_w_2;
gp_training_opts.kernel_params.depth = depth;

mu = 0;
dotProdSum = 0;
testPoint = get_cartesian_coordinates(testPoint,[]);
for ii = 1:size(training_data,1)
    [K,~] = getKernel(training_data(ii,:),testPoint,gp_training_opts,false);
    mu = mu +  K *y_i_vec(ii) ;
end

end