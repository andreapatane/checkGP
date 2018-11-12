%%%
%Script for analysis and plotting of Figure 4.

clear all
close all
clc
addpath(genpath('utils/'))

depths = 1:10; %parametric analysis on kernel depth from 1 to 10
training_set_samples = 100:100:2000; %similiarly for training set size.
testPointIdx = 1; %selecting only first test sample

gp_training_opts.kernel = 'ReLU'; 
gp_training_opts.dataSet = 'mnist';
gp_training_opts.num_test_samples = 5;
gp_training_opts.output_mode = 'class_via_regress';

gp_training_opts.data_folder = '../data/flattened_inputs/mnist/';
gp_training_opts.x_train_file = 'x_train2000.csv';
gp_training_opts.y_train_file = 'y_train2000.csv';
gp_training_opts.x_test_file = 'x_test.csv';
gp_training_opts.y_test_file = 'y_test.csv';

gp_training_opts.kernel_params.sigma_w_2 = 3.19;
gp_training_opts.kernel_params.sigma_b_2 = 0.00;
gp_training_opts.kernel_params.sig_eps = 1e-10;


var_comp_opts.epsilons = 0.10;
var_comp_opts.invariance = false; %this is completely unrelevant in this case.

var_comp_opts.delta = 0:0.0001:1.0;
var_comp_opts.pix_2_mod = 'sift';
var_comp_opts.vl_setup_script_path = '../vlfeat/toolbox/';
var_comp_opts.kernel = gp_training_opts.kernel;
var_comp_opts.maxIterations_sigma = 20;
 %For variance analysis I give the tollerance in terms of percentage value
 %from the original one. This is done because the order of magnitude of
 %variance changes significantly dependning on the number of training
 %samples used. Also the plot compares value in terms of percentage value
 %wrt to test set, so it is the natural way to give tollerance here.
var_comp_opts.tollerance_sigma_percentage = 0.5;
var_comp_opts.constrain_2_one = true;
var_comp_opts.output_mode = gp_training_opts.output_mode;
var_comp_opts.max_feats = 1; %just using first feature.
var_comp_opts.max_size = 11;


nameTempl = 'mnist_var_analysis';
num_of_threads = 1;

var_cell_param_analysis = cell(length(training_set_samples),length(depths));

for ii = 1:length(training_set_samples)
    gp_training_opts.num_ts = training_set_samples(ii);
    
    
    for jj = 1:length(depths)
        
        gp_training_opts.kernel_params.depth = depths(jj);
        var_comp_opts.kernel_params = gp_training_opts.kernel_params;

        [test_data,test_norm,test_labels,mu_test,...
            Kstar,training_data,trainedSystem,K_inv,R,gp_training_opts,sigma_test] = main_trainGP(gp_training_opts);
     
        
        if strcmp(gp_training_opts.dataSet,'mnist')
            test_labels = test_labels + 1;
        end
             
        var_comp_opts.ref_value = sigma_test(testPointIdx,testPointIdx);
        var_cell_param_analysis{ii,jj} = main_variance_computation(var_comp_opts,...
            num_of_threads,testPointIdx,testPointIdx,...
            training_data,test_data,test_norm,test_labels,mu_test,Kstar,trainedSystem,K_inv,R);
    end
end


