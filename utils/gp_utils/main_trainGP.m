function [test_data,test_norm,test_labels,mu_test,...
    Kstar,training_data,trainedSystem,K_inv,R,gp_training_opts,sigma_test]=main_trainGP(gp_training_opts)
%main_trainGP is the main wrapper function for GP training. All the inputs
%are structured in the gp_training_opts MATLAB structure. The fields of the
%latter include:
% - dataSet: Specifies the dataset to be used. (Current implementation is either 'mnist' or 'x_products' [the latter is the running example on the paper] )
% - data_folder: [Used for the mnist dataset] Specifies the directory in
% which training and testing data are stored.
% - x_train_file: [Used for the mnist dataset] name of the .csv file in
% which training set mnist image are stored (similarly for y_train_file, x_test_file, y_test_file)
% - num_ts: number of training samples 
% - num_test_samples: number of test samples
% - kernel: string for the kernel to use inside the GP model (current implementation account only for 'ReLU' and 'sqe').
% - kernel_params: sub-structure for kernel specific hyper-parameters.
% - output_mode: If 'regression' then the GP is used to solve a regression
% problem. If 'class_via_regress' then a regression GP model is still
% trained, but for classification purposes.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loading training and test set%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(gp_training_opts.dataSet,'mnist')
    training_data = csvread([gp_training_opts.data_folder,gp_training_opts.x_train_file]);
    training_labels = csvread([gp_training_opts.data_folder,gp_training_opts.y_train_file]);
    training_data = training_data(1:gp_training_opts.num_ts,:);
    training_labels = training_labels(1:gp_training_opts.num_ts,:);
    [training_data,~]  = averagePooling(training_data); %subsampling pixels by a factor of 2.
    test_data = csvread([gp_training_opts.data_folder,gp_training_opts.x_test_file]);
    test_data = test_data(1:gp_training_opts.num_test_samples,:);
    test_data = averagePooling(test_data); %subsampling pixels by a factor of 2.
    test_labels = csvread([gp_training_opts.data_folder,gp_training_opts.y_test_file]);
    test_labels = test_labels(1:gp_training_opts.num_test_samples);
elseif strcmp(gp_training_opts.dataSet,'x_products')
    [training_data,training_labels,test_data,test_labels] = x_product_data(gp_training_opts); %generating data for the x_products example
end






if strcmp(gp_training_opts.kernel,'ReLU')
    training_data = to_norm1(training_data); %projecting every point to the unit sphere for the ReLU kernel
    [test_data,test_norm] = to_norm1(test_data);
else
    test_norm = [];
end

disp('training the GP...')
[trainedSystem,R,K_inv,gp_training_opts] = trainGP(training_data,training_labels,gp_training_opts);

disp('Computing kernels on test set')
[Kstar,Kstarstar] = compute_kernel_on_test(test_data,training_data,gp_training_opts);

disp('finally making predictions')
[mu_test, sigma_test,y_test_hat] = make_predictions_on_test(Kstar,trainedSystem,Kstarstar,R,gp_training_opts.output_mode);



if strcmp(gp_training_opts.output_mode,'class_via_regress')
    
    test_accuracy = mean(y_test_hat ==  test_labels);
    disp('Accuracy obtained on test set:')
    disp(test_accuracy)
    
    disp('I will keep only correctly classified images for further analysis...') %Is there any point in looking for adversarial in images that are already mis-classified???
    goodIdxs = find(y_test_hat ==  test_labels);
    test_labels = test_labels(goodIdxs);
    test_data = test_data(goodIdxs,:);
    test_norm = test_norm(goodIdxs);
    mu_test = mu_test(goodIdxs,:);
    Kstar = Kstar(goodIdxs,:);
    
elseif strcmp(gp_training_opts.output_mode,'regression')
    
    rmse = sqrt(  mean( (mu_test - test_labels).^2   )  );
    disp('root mean squared error obtained on test set:')
    disp(rmse)
    
end




end



