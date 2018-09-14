function maxVarianceAnalysis(testPointIdx)
%
%Main function for performing analysis of maximal variance
%

maxNumCompThreads(get_num_of_threads());

doAnalyses = true;
dataSet = 'mnist';
depth_vec = [1:1:10];
sigma_w_2_vec = [3.19];
actualTrainingSamples_vec = [100:100:2000];
lls = zeros(length(sigma_w_2_vec),length(actualTrainingSamples_vec),length(depth_vec));
uus = zeros(length(sigma_w_2_vec),length(actualTrainingSamples_vec),length(depth_vec));
refs = zeros(length(sigma_w_2_vec),length(actualTrainingSamples_vec),length(depth_vec));
for ii_sigma = 1:length(sigma_w_2_vec)
    for ii_samples = 1:length(actualTrainingSamples_vec)
        for ii_depth = 1:length(depth_vec)
        disp([ii_sigma,ii_samples,ii_depth])
        sigma_w_2 = sigma_w_2_vec(ii_sigma);
        actualTrainingSamples = actualTrainingSamples_vec(ii_samples);
        depth = depth_vec(ii_depth);
        [sigma_prior,theta_vec,test_data,test_norm,test_labels,mu_test,kernel,...
            Kstar,training_data,trainedSystem,K_inv,actualTrainingSamples,R,std_Y,sigma_test] = ...
            main_trainGP(dataSet,depth,sigma_w_2,actualTrainingSamples);
        if doAnalyses
            training_dataInit = training_data;
            %and clearing up memory space.
            clear training_data
            epsilon = 0.15;
            testPoint = test_data(testPointIdx,:);
            testPointInit = testPoint;
            testLabel = test_labels(testPointIdx);
            test_prediction = mu_test(testPointIdx,:);
            [pixel2Modify,magnitudes] = sift_on_mnist(testPoint,10,11);
            [~, max_idx] = max(magnitudes);
            pixel2Modify = pixel2Modify(max_idx);
            if strcmp(kernel,'N-ReLU')
                normTestPoint = test_norm(testPointIdx);
            else
                normTestPoint = 1;
            end
            if strcmp(kernel,'sqe')
                r_star = Kstar(testPointIdx,:)'/sigma_prior;
            elseif strcmp(kernel,'N-ReLU')
                r_star = Kstar(testPointIdx,:)';
            end
            
            for hh = 1:length(pixel2Modify)
                
                if strcmp(kernel,'N-ReLU')
                    
                    %I re-arrange pixels so that the current feature is at the beginning of the point vector
                    testPoint = prioritize_pixel2mod(testPointInit,pixel2Modify{hh} );
                    training_data = prioritize_pixel2mod(training_dataInit,pixel2Modify{hh} );
                    
                    % getting sperical representation of training data
                    training_data = get_spherical_coordinates(training_data,[]);
                    
                    %now the feature is just at the beginning of the vector so I update
                    %its position vector:
                    pix2Mod =  1:length(pixel2Modify{hh});
                elseif strcmp(kernel,'sqe')
                    %in the case of sqe, I do not do much here...
                    testPoint =  testPointInit;
                    training_data = training_dataInit;
                    pix2Mod = pixel2Modify{hh};
                end
                
                
                %I compute the hyper-rectangle on the cartesian coordinate around
                %testPoint
                [x_L, x_U] = compute_hyper_rectangle(epsilon,testPoint*normTestPoint, pix2Mod,dataSet);
                
                %if I am using relu kernel, then I move to spherical space, and
                %over-approximate the cartesian hyper-cube with an hypercube in
                %spherical coordinates.
                if strcmp(kernel,'N-ReLU')
                    testPoint = get_spherical_coordinates(testPoint,[]);
                    %[x_L,x_U] = compute_angle_space_hyper_rectangle(x_L,x_U);
                end
                ref_value = sigma_test(testPointIdx,testPointIdx);
                [lb_sigma,ub_sigma] = b_n_b_for_lb_computation(K_inv,x_L,x_U,theta_vec,sigma_prior,training_data,...
                    testLabel + 1,testPoint,'sigma',test_prediction,r_star,[],kernel,dataSet,depth,ref_value);
                lls(ii_sigma,ii_samples,ii_depth) = lb_sigma;
                uus(ii_sigma,ii_samples,ii_depth) = ub_sigma;
                refs(ii_sigma,ii_samples,ii_depth) = ref_value;
            end
        end
        end
    end
end
if doAnalyses
    T.lls = lls;
    T.uus = uus;
    T.refs = refs;
    T.sigma_w_2_vec = sigma_w_2_vec;
    T.actualTrainingSamples_vec = actualTrainingSamples_vec;
    T.testPointIdx = testPointIdx;
    T.epsilon = epsilon;
    T.depth = depth_vec;
    result_folder = get_result_folder();
    save([result_folder,'variance_and_depth_ReLU_mnist_',int2str(testPointIdx),'.mat'],'T');
end
end
