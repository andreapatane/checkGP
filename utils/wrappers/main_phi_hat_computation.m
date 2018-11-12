function [p_ub_cell,sup_mu,sup_xi,pixel2Modify_cell] = main_phi_hat_computation(upper_bound_comp_opts,numOfThreads,startIdx,endIdx,training_dataInit,test_data,test_norm,...
    test_labels,mu_test,Kstar,trainedSystem,K_inv,R)
%main_phi_hat_computation: main function for the computation of
%\hat{\phi}_1 and \hat{\phi}_2.
%%%%%
%Inputs:
% - upper_bound_comp_opts: Matlab structure used to specify settings for the computations performed.
% - numOfThreads: Number of concurrent threads used by quadprog.
% - startIdx and endIdx: given the test set "test_data", the analysis are
%       performed on the test points from row startIdx to row endIdx.
% - training_dataInit: training data matrix
% - test_data: test data matrix
% - test_norm: [used only for ReLU kernel], specify the norm of the test
%       point (prior to projection to the unit sphere)
% - test_labels: labels of test set
% - mu_test: a-posteriori mean on test set
% - Kstar: a-posteriori covariance between test and train test
% - trainedSystem: linear coefficients on a-posteriori mean estimation
% - K_inv: inverse of the a-posteriori covariance function on the training
%       set
% - R: choleski decomposition of covariance matrix



%setting maximum number of concurrent threads to be used by quadprog
maxNumCompThreads(numOfThreads);

%keeping only test points that we want to analyse.
test_data = test_data(startIdx:endIdx,:);
if strcmp(upper_bound_comp_opts.kernel,'ReLU')
    test_norm = test_norm(startIdx:endIdx);
end
test_labels = test_labels(startIdx:endIdx);
mu_test = mu_test(startIdx:endIdx,:);
Kstar =  Kstar(startIdx:endIdx,:);

pixel2Modify_cell = cell((endIdx-startIdx+1),1);
p_ub_cell = cell((endIdx-startIdx+1),1);

%for loop on test points that need to be analysed
for testPointIdx = 1:(endIdx-startIdx+1)
    testPoint = test_data(testPointIdx,:);
    testPointInit = testPoint;
    
    
    
    
    if strcmp(upper_bound_comp_opts.pix_2_mod,'sift')
        %computing sift features
        pixel2Modify = sift_on_mnist(testPoint,upper_bound_comp_opts.max_feats,upper_bound_comp_opts.max_size,upper_bound_comp_opts.vl_setup_script_path);
    else
        pixel2Modify = upper_bound_comp_opts.pix_2_mod;
    end
    pixel2Modify_cell{testPointIdx} = pixel2Modify;
    
    
    p_ub = cell(length(pixel2Modify),1);
    
    if strcmp(upper_bound_comp_opts.kernel,'ReLU')
        normTestPoint = test_norm(testPointIdx);
    else
        normTestPoint = 1;
    end
    testLabel = test_labels(testPointIdx);
    test_prediction = mu_test(testPointIdx,:);
    
    if strcmp(upper_bound_comp_opts.kernel,'sqe') %If sqe kernel, than I normalise everything by the variance hyperparameter.
        r_star = Kstar(testPointIdx,:)'/upper_bound_comp_opts.kernel_params.sigma;
    elseif strcmp(upper_bound_comp_opts.kernel,'ReLU') %else I do nothing...
        r_star = Kstar(testPointIdx,:)';
    end
    
    %vector to save up partial results
    Xs = cell(length(pixel2Modify),2); %hyper-rectangle T
    sup_mu = zeros(size(pixel2Modify)); %sup on a-posteriori mean
    sup_xi = zeros(size(pixel2Modify)); %sup on a-posteriori \xi
    
    
    
    for ll = 1:length(pixel2Modify) %for loops ranging across the various features
        
        disp(['ll: ', int2str(ll)])
        %initilising result matrix:
        p_ub{ll} = zeros(length(upper_bound_comp_opts.delta),length(upper_bound_comp_opts.epsilons));
        
        
        for jj = 1:length(upper_bound_comp_opts.epsilons)
            
            if strcmp(upper_bound_comp_opts.kernel,'ReLU') %in case of ReLU kernel I do some extra work to facilate working on hyper-spherical coordinates
                %I re-arrange pixels so that the current feature is at the beginning of the point vector
                testPoint = prioritize_pixel2mod(testPointInit,pixel2Modify{ll} );
                training_data = prioritize_pixel2mod(training_dataInit,pixel2Modify{ll} );
                
                % getting sperical representation of training data
                training_data = get_spherical_coordinates(training_data,[]);
                
                %now the feature is just at the beginning of the vector so I update
                %its position vector:
                pix2Mod =  1:length(pixel2Modify{ll});
            elseif strcmp(upper_bound_comp_opts.kernel,'sqe')
                %in the case of sqe, I do not do much here...
                testPoint =  testPointInit;
                training_data = training_dataInit;
                pix2Mod = pixel2Modify{ll};
            end
            disp(['jj: ',int2str(jj)])
            
            %I compute the hyper-rectangle on the cartesian coordinate around
            %testPoint
            [x_L, x_U] = compute_hyper_rectangle(upper_bound_comp_opts.epsilons(jj),testPoint*normTestPoint, pix2Mod,upper_bound_comp_opts.constrain_2_one);
            
            if strcmp(upper_bound_comp_opts.kernel,'ReLU') %if I am using relu kernel, then I move to spherical space.
                testPoint = get_spherical_coordinates(testPoint,[]);
            end
            
            
            %depending on the task, I define the relevant output index
            if strcmp(upper_bound_comp_opts.output_mode,'class_via_regress')
                outIdx = testLabel;
            elseif strcmp(upper_bound_comp_opts.output_mode,'regression')
                outIdx = 1;
            end
            
            Xs{ll,1} = x_L;
            Xs{ll,2} = x_U;
            %%%%%
            %%%branch and bound for \mu computations
            disp('mu:')
            [~,ub_mu] = b_n_b_for_lb_computation(trainedSystem,x_L,x_U,training_data,...
                outIdx,testPoint,'mu',test_prediction,[],[],upper_bound_comp_opts);
            
            if upper_bound_comp_opts.invariance %in case of invariance property, I also look at the maximum variance difference
                [~,ub_mu_minus] = b_n_b_for_lb_computation(-trainedSystem,x_L,x_U,training_data,...
                    outIdx,testPoint,'mu',-test_prediction,[],[],upper_bound_comp_opts);
                ub_mu = max(ub_mu,ub_mu_minus); % and then I keep the worst case.
            end
            sup_mu(ll) = ub_mu;
            
            %%%%%
            %%%branch and bound for \xi computations
            disp('xi:')
            [~,ub_xi] = b_n_b_for_lb_computation(K_inv,x_L,x_U,training_data,...
                outIdx,testPoint,'xi',test_prediction,r_star,R,upper_bound_comp_opts);
            sup_xi(ll) = ub_xi;
            
            
            %computing the K constant
            if strcmp(upper_bound_comp_opts.kernel,'sqe')
                K = 1.414*upper_bound_comp_opts.kernel_params.sigma;
            elseif strcmp(upper_bound_comp_opts.kernel,'ReLU')
                [K,~,alpha_U] =  compute_K_constant_for_ReLU(x_L,x_U,upper_bound_comp_opts.kernel_params.sigma_b_2,...
                    upper_bound_comp_opts.kernel_params.sigma_w_2,length(x_L)+1,upper_bound_comp_opts.kernel_params.depth);
            end
            
            %overapproximating sup_d
            sup_d = 2.0*ub_xi;
            
            %define the integrand function
            if strcmp(upper_bound_comp_opts.kernel,'sqe')
                fun = @(z)( sqrt( length(pix2Mod)* log( sqrt(length(pix2Mod))*K*2*upper_bound_comp_opts.epsilons(jj)./z + 1 )  ) );
            elseif strcmp(upper_bound_comp_opts.kernel,'ReLU')
                fun = @(z)( sqrt( log( 2*K*alpha_U./z + 1 )  ) );
            end
            %numerical computation of integral
            eta_int = 12*integral(fun,0,0.5*sup_d);
            
            
            %final computations for the over-approximations
            for ii = 1:length(upper_bound_comp_opts.delta)
                eta = upper_bound_comp_opts.delta(ii) - (ub_mu + eta_int);
                if eta <=0
                    p_ub{ll}(ii,jj) = 1;
                else
                    p_ub{ll}(ii,jj) = exp( - eta^2 / (2*ub_xi) );
                    if upper_bound_comp_opts.invariance
                        p_ub{ll}(ii,jj) = min(2*p_ub{ll}(ii,jj),1);
                    end
                end
            end
            
            
        end
        
    end
    p_ub_cell{testPointIdx} = p_ub;
    
end
end