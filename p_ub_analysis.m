function p_ub_analysis(depth,startIdx)
%p_ub_analysis is the main function for the experiments to computed
%over-approximation for adversarial.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input depth gives the depth of the ReLU kernel (in case this is used)
%Input startIdx specify the first test image to start analysing.
%No output is defined. Results are actually saved to disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%providing default inputs if none were provided
if nargin < 2
    startIdx = 1;
    if nargin < 1
        depth = 1;
    end
end

addpath(genpath('utils/'))

%Setting the maximum num of threads, accordingly to what provided in the
%get_num_of_threads util function. Parallel toolbox may be internatlly used
%when solving quadratic programs
maxNumCompThreads(get_num_of_threads());

%Flag for invariance analysis. When false, one-sided adversarial analysis
%is performed. Otherwise invariance property is over-approximated.
invariance = false;

%dataSet string name. This defines the settings for the experiments used in
%the paper. That is 'x_products' corresponds to the running example on GP
%regression, while 'mnist' corresponds to the analysis on the MNIST
%dataset.
dataSet = 'mnist';


%performing training of the GP
[sigma_prior,theta_vec,test_data,test_norm,test_labels,mu_test,kernel,...
    Kstar,training_data,trainedSystem,K_inv,actualTrainingSamples,R] = main_trainGP(dataSet,depth);


%defining epsilons and delta values for analysis
epsilons = [0.05,0.1,0.15]; 
delta = 0:0.01:1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Saving the initial testPoint vector so that I can modify testPoint in the
%for loop without risk of headache...
%doing the same with the training data
training_dataInit = training_data;
%and clearing up memory space.
clear training_data


%looping over the first 100 test point, and performing all analyses on them iteratively 
for testPointIdx = startIdx:100
    testPoint = test_data(testPointIdx,:);
    testPointInit = testPoint;
    
    %%%%%%%%%%%
    %below we define pixel2Modify. This is a Cell array of image patches,
    %of the form {featVec1,featVec2,....,featVecK}, where each featVec is a
    %vector of pixel indexes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(dataSet,'mnist')
        try
            %try to use sift, if this is not installed, I require manual
            %entries for features
            pixel2Modify = sift_on_mnist(testPoint,5,11);
        catch ME
            warning(ME.message)
            error('Sift not working, because of message above. Please manually insert features here.');
        end
    elseif strcmp(dataSet,'x_products')
        pixel2Modify = {[1,2]};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %p_ub is the cell array that will store the analyses results
    p_ub = cell(length(pixel2Modify),1);
    
    %Taking into account initial norm is using ReLu kernel...
    if strcmp(kernel,'N-ReLU')
        normTestPoint = test_norm(testPointIdx);
    else
        normTestPoint = 1;
    end
    testLabel = test_labels(testPointIdx);
    test_prediction = mu_test(testPointIdx,:);
    
    %If sqe kernel, than I normalise everything by the variance
    %hyperparameter, as the latter is kind of shared by everything and
    %cancels out easily.
    if strcmp(kernel,'sqe')
        r_star = Kstar(testPointIdx,:)'/sigma_prior;
    elseif strcmp(kernel,'N-ReLU')
        r_star = Kstar(testPointIdx,:)';
    end
    
    %aux vectors to store partial results
    Xs = [];
    sup_mu = [];
    sup_xi = [];
    %
    for ll = 1:length(pixel2Modify) %for loops ranging across the various "features"
        
        disp(['ll: ', int2str(ll)])
        %initilising result matrix:
        p_ub{ll} = zeros(length(delta),length(epsilons));
        %
        for jj = 1:length(epsilons)
            %
            if strcmp(kernel,'N-ReLU')
                %%%%%%%%%
                %I re-arrange pixels so that the current feature is at the beginning of the point vector
                %This helps in working with hyper-spherical coordinates
                testPoint = prioritize_pixel2mod(testPointInit,pixel2Modify{ll} );
                training_data = prioritize_pixel2mod(training_dataInit,pixel2Modify{ll} );
                %%%%%%%%%
                % getting sperical representation of training data
                training_data = get_spherical_coordinates(training_data,[]);
                %now the feature is just at the beginning of the vector so I update
                %its position vector:
                pix2Mod =  1:length(pixel2Modify{ll});
                %
            elseif strcmp(kernel,'sqe')
                %in the case of sqe, I do not do much here... just making
                %some copies really
                testPoint =  testPointInit;
                training_data = training_dataInit;
                pix2Mod = pixel2Modify{ll};
            end
            disp(['jj: ',int2str(jj)])
            
            %I compute the hyper-rectangle on the cartesian coordinate
            %around the testPoint under analysis.
            [x_L, x_U] = compute_hyper_rectangle(epsilons(jj),testPoint*normTestPoint, pix2Mod,dataSet);
            
            %if I am using relu kernel, I now convert the testPoint into
            %spherical coordinates.
            if strcmp(kernel,'N-ReLU')
                testPoint = get_spherical_coordinates(testPoint,[]);
            end
            
            
            %defining the output index that I want to analyse
            if strcmp(dataSet,'mnist')
                outIdx = testLabel + 1;
            elseif strcmp(dataSet,'x_products')
                outIdx = 1;
            end
            
            Xs = [Xs;x_L;x_U];
            %branch and bound on the mean
            disp('Performing mu computation:')
            [~,ub_mu] = b_n_b_for_lb_computation(trainedSystem,x_L,x_U,theta_vec,sigma_prior,training_data,...
                outIdx,testPoint,'mu',test_prediction,[],[],kernel,dataSet,depth);
            if invariance
                %in invariance case I need the sup of the absolute value of
                %mean. So I compute also the min of -mean and get the
                %maximum of both.
                [lb_mu_minus,ub_mu_minus] = b_n_b_for_lb_computation(-trainedSystem,x_L,x_U,theta_vec,sigma_prior,training_data,...
                    outIdx,testPoint,'mu',-test_prediction,[],[],kernel,dataSet,depth);
                ub_mu = max(ub_mu,ub_mu_minus);
            end
            sup_mu = [sup_mu;ub_mu];
            
            disp('Performing xi computation:')
            %branch and bound on \xi
            [~,ub_xi] = b_n_b_for_lb_computation(K_inv,x_L,x_U,theta_vec,sigma_prior,training_data,...
                outIdx,testPoint,'xi',test_prediction,r_star,R,kernel,dataSet,depth);
            sup_xi = [sup_xi;ub_xi];
            
            %computing K here:
            if strcmp(kernel,'sqe')
                K = 1.414*sigma_prior;
            elseif strcmp(kernel,'N-ReLU')
                [K,~,alpha_U] =  compute_K_constant_for_ReLU(x_L,x_U,theta_vec,sigma_prior,length(x_L)+1,depth);
            end
            
            %overapproximating sup_d here:
            sup_d = 2.0*ub_xi;
            
            %defining the integrand function here:
            if strcmp(kernel,'sqe')
                fun = @(z)( sqrt( length(pix2Mod)* log( 4*K*(jj)./z + 1 )  ) );
            elseif strcmp(kernel,'N-ReLU')
                fun = @(z)( sqrt( log( 2*K*alpha_U./z + 1 )  ) );
            end
            %numerical 1-d integration:
            eta_int = 12*integral(fun,0,0.5*sup_d);
            
            
            %simple algebra to get final value of p_ub:
            for ii = 1:length(delta)
                eta = delta(ii) - (ub_mu + eta_int);
                if eta <=0
                    p_ub{ll}(ii,jj) = 1;
                else
                    p_ub{ll}(ii,jj) = exp( - eta^2 / (2*ub_xi) );
                    if invariance
                        p_ub{ll}(ii,jj) = min(2*p_ub{ll}(ii,jj),1);
                    end
                end
            end
            
        end
    end
    
    %surf(epsilons,delta,p_ub{1})
    %plot(delta,p_ub{1})
    %ylabel('delta')
    %xlabel('epsilon')
    
    %defining structure to store on disk
    T.trainSamples = actualTrainingSamples;
    T.testPointIdx = testPointIdx;
    T.sigma_prior = sigma_prior;
    T.theta = theta_vec(1);
    T.pixel2Modify = pixel2Modify;
    T.epsilons = epsilons;
    T.delta = delta;
    T.p_ub = p_ub;
    T.sup_mu = sup_mu;
    T.sup_xi = sup_xi;
    result_folder = get_result_folder();
    save([result_folder,'p_ub_ReLU_mnist_',int2str(depth),'_',int2str(testPointIdx),'.mat'],'T');
    
end
end