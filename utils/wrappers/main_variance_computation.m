function var_cell = main_variance_computation(var_comp_opts,numOfThreads,startIdx,endIdx,training_dataInit,test_data,test_norm,...
    test_labels,mu_test,Kstar,trainedSystem,K_inv,R)
%main_variance_computation is the main wrapper function for the variance
%computation. It's very similar to main_phi_hat_computation


%setting maximum number of concurrent threads to be used by quadprog
maxNumCompThreads(numOfThreads);

%keeping only test points that we want to analyse.
test_data = test_data(startIdx:endIdx,:);
if strcmp(var_comp_opts.kernel,'ReLU')
    test_norm = test_norm(startIdx:endIdx);
end
test_labels = test_labels(startIdx:endIdx);
mu_test = mu_test(startIdx:endIdx,:);
Kstar =  Kstar(startIdx:endIdx,:);

pixel2Modify_cell = cell((endIdx-startIdx+1),1);
var_cell = cell((endIdx-startIdx+1),1);
for testPointIdx = 1:(endIdx-startIdx+1)
    testPoint = test_data(testPointIdx,:);
    testPointInit = testPoint;

    
    
    
    if strcmp(var_comp_opts.pix_2_mod,'sift')
        pixel2Modify = sift_on_mnist(testPoint,var_comp_opts.max_feats,var_comp_opts.max_size,var_comp_opts.vl_setup_script_path);
    else
        pixel2Modify = var_comp_opts.pix_2_mod;
    end
    pixel2Modify_cell{testPointIdx} = pixel2Modify;
    
       
    var_sub_cell = cell(length(pixel2Modify),1);
 
    if strcmp(var_comp_opts.kernel,'ReLU')
        normTestPoint = test_norm(testPointIdx);
    else
        normTestPoint = 1;
    end
    testLabel = test_labels(testPointIdx);
    test_prediction = mu_test(testPointIdx,:);
    
    
    %If sqe kernel, than I normalise by the variance hyperparameter.
    if strcmp(var_comp_opts.kernel,'sqe')
        r_star = Kstar(testPointIdx,:)'/var_comp_opts.kernel_params.sigma;
    elseif strcmp(var_comp_opts.kernel,'ReLU')
        r_star = Kstar(testPointIdx,:)';
    end
    
  
    
    
    for ll = 1:length(pixel2Modify) %for loops ranging across the various features
        
        disp(['ll: ', int2str(ll)])
        %initilising result matrix:
        var_sub_cell{ll} = zeros(length(var_comp_opts.epsilons),2);

        
        for jj = 1:length(var_comp_opts.epsilons)
            if strcmp(var_comp_opts.kernel,'ReLU')
                
                %I re-arrange pixels so that the current feature is at the beginning of the point vector
                testPoint = prioritize_pixel2mod(testPointInit,pixel2Modify{ll} );
                training_data = prioritize_pixel2mod(training_dataInit,pixel2Modify{ll} );
                
                % getting sperical representation of training data
                training_data = get_spherical_coordinates(training_data,[]);
                
                %now the feature is just at the beginning of the vector so I update
                %its position vector:
                pix2Mod =  1:length(pixel2Modify{ll});
            elseif strcmp(var_comp_opts.kernel,'sqe')
                %in the case of sqe, I do not do much here...
                testPoint =  testPointInit;
                training_data = training_dataInit;
                pix2Mod = pixel2Modify{ll};
            end
            disp(['jj: ',int2str(jj)])
            
            %I compute the hyper-rectangle on the cartesian coordinate around
            %testPoint
            [x_L, x_U] = compute_hyper_rectangle(var_comp_opts.epsilons(jj),testPoint*normTestPoint, pix2Mod,var_comp_opts.constrain_2_one);
            
            %if I am using relu kernel, then I move to spherical space.
            if strcmp(var_comp_opts.kernel,'ReLU')
                testPoint = get_spherical_coordinates(testPoint,[]);
            end
            
            
            %branch and bound for \mu^0 computations
            if strcmp(var_comp_opts.output_mode,'class_via_regress')
                outIdx = testLabel;
            elseif strcmp(var_comp_opts.output_mode,'regression')
                outIdx = 1;
            end
            
            
            disp('sigma:')
            [lb_sigma,ub_sigma] = b_n_b_for_lb_computation(K_inv,x_L,x_U,training_data,...
                outIdx,testPoint,'sigma',test_prediction,r_star,R,var_comp_opts);
            var_sub_cell{ll}(jj,1) = lb_sigma;
            var_sub_cell{ll}(jj,2) = ub_sigma;
            
            
            
        end
        
    end
    var_cell{testPointIdx} = var_sub_cell;
end
end