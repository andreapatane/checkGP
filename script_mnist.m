%%%
%%Analysis and plotting for Figure 3.

clear all
close all
clc
addpath(genpath('utils/'))

%%%some plotting parameters
x_lim_max_plot = 0.5;
colors = {[255 192 203]/255, [173 216 230]/255, [170, 227, 60]/255, [255,165,0]/255,[147, 112, 216]/255 };
%%%

gp_training_opts.kernel = 'ReLU'; %kernel to use
gp_training_opts.dataSet = 'mnist'; %dataset mnist
gp_training_opts.kernel_params.depth = 2; % depth of ReLU kernel 
gp_training_opts.num_ts = 1000;  %number of training samples
gp_training_opts.num_test_samples = 4; %number of test samples used
gp_training_opts.output_mode = 'class_via_regress'; %task is that of classification via regression

%gp_training_opts.data_folder = get_data_folder();
gp_training_opts.data_folder = '../data/flattened_inputs/mnist/'; %directory where mnist data are stored
gp_training_opts.x_train_file = 'x_train2000.csv'; % file name for domain values of training images. This is expected to be a Nx784 matrix, where N is the number of samples.
gp_training_opts.y_train_file = 'y_train2000.csv'; % file name for labels of training images. This is expected to be a N-dimensional vector of integer values between 0 and 9  
gp_training_opts.x_test_file = 'x_test.csv'; % file name for domain values of test images
gp_training_opts.y_test_file = 'y_test.csv'; % file name for labels of testing images.

%hyper-parameter values used for ReLU kernel
gp_training_opts.kernel_params.sigma_w_2 = 3.19; 
gp_training_opts.kernel_params.sigma_b_2 = 0.00;
gp_training_opts.kernel_params.sig_eps = 1e-10;


%%%%%%
%%%Performing training of the GP
[test_data,test_norm,test_labels,mu_test,...
    Kstar,training_data,trainedSystem,K_inv,R,gp_training_opts] = main_trainGP(gp_training_opts);




%%%%%
%Setting up options used for the analysis of \hat{\phi}
upper_bound_comp_opts.epsilons = [0.05,0.15]; %performing analysis for hyper-rectangle values of 0.05 and 0.15
upper_bound_comp_opts.invariance = false; %true: variance property, false: safety property

upper_bound_comp_opts.delta = 0:0.0001:1.0; %values of delta used 
upper_bound_comp_opts.pix_2_mod = 'sift'; %asking for sift features to be computed (this is possbile only if sift has been installed ofc...)
upper_bound_comp_opts.vl_setup_script_path = '../vlfeat/toolbox/'; %relative path to vlfeat toolbox (used for SIFT)
upper_bound_comp_opts.kernel = gp_training_opts.kernel; %updating kernel hyper-parameters
upper_bound_comp_opts.kernel_params = gp_training_opts.kernel_params; %updating kernel hyper-parameters
upper_bound_comp_opts.maxIterations_mu = 5000; %setting maximum number of iterations for bnb on \mu
upper_bound_comp_opts.maxIterations_xi = 25; %setting maximum number of iterations for bnb on \xi
upper_bound_comp_opts.max_feats = 5; %maximum number of features per image
upper_bound_comp_opts.max_size = 11; %maximum size of feature

upper_bound_comp_opts.tollerance_mu = 0.0005; %tollerance on mu optimal value
upper_bound_comp_opts.tollerance_xi = 0.0002; %tollerance on xi optimal value
upper_bound_comp_opts.constrain_2_one = true;

upper_bound_comp_opts.output_mode = gp_training_opts.output_mode; 
nameTempl = 'mnist_1'; %string name to store analysis results into disk
num_of_threads = 1; %number of concurrent thread for quaprog

if strcmp(gp_training_opts.dataSet,'mnist')
    test_labels = test_labels + 1; %just doing some useful transformation on mnist labels (so that they start from 1 and can be used for matrix indexing)
end


%%%%
%Computation of upper-bound on first test point
%
[p_ub,~,~,~] = main_phi_hat_computation(upper_bound_comp_opts,...
    num_of_threads,1,1,...
    training_data,test_data,test_norm,test_labels,mu_test,Kstar,trainedSystem,K_inv,R);
%%%%
%Plotting results
%
figure 
hold on
grid on
title('\gamma =  0.05');
for ii = 1:length(p_ub{1})
    plot(upper_bound_comp_opts.delta,p_ub{1}{ii}(:,1),'Color',colors{ii},'LineWidth',3.0)
end
xlim([0,x_lim_max_plot])
hold off

figure 
hold on
grid on
title('\gamma =  0.15');
for ii = 1:length(p_ub{1})
    plot(upper_bound_comp_opts.delta,p_ub{1}{ii}(:,2),'Color',colors{ii},'LineWidth',3.0)
end
xlim([0,x_lim_max_plot])
hold off

%%%%
%Computation of upper-bound on second test point
%
[p_ub,~,~,~] = main_phi_hat_computation(upper_bound_comp_opts,...
    num_of_threads,2,2,...
    training_data,test_data,test_norm,test_labels,mu_test,Kstar,trainedSystem,K_inv,R);
figure 
hold on
grid on
%%%%
%Plotting results
%
for ii = 1:length(p_ub{1})
    plot(upper_bound_comp_opts.delta,p_ub{1}{ii}(:,1),'Color',colors{ii},'LineWidth',3.0)
end
xlim([0,x_lim_max_plot])
hold off
figure 
hold on
grid on
title('\gamma =  0.15');
for ii = 1:length(p_ub{1})
    plot(upper_bound_comp_opts.delta,p_ub{1}{ii}(:,2),'Color',colors{ii},'LineWidth',3.0)
end
xlim([0,x_lim_max_plot])
hold off

%%%%
%Computation of upper-bound on fourth test point [I skip the third because SIFT finds only 3 significant features on it...]
%

[p_ub,~,~,~] = main_phi_hat_computation(upper_bound_comp_opts,...
    num_of_threads,4,4,...
    training_data,test_data,test_norm,test_labels,mu_test,Kstar,trainedSystem,K_inv,R);
figure 
hold on
grid on

for ii = 1:length(p_ub{1})
    plot(upper_bound_comp_opts.delta,p_ub{1}{ii}(:,1),'Color',colors{ii},'LineWidth',3.0)
end
xlim([0,x_lim_max_plot])
hold off
figure 
hold on
grid on
title('\gamma =  0.15');
for ii = 1:length(p_ub{1})
    plot(upper_bound_comp_opts.delta,p_ub{1}{ii}(:,2),'Color',colors{ii},'LineWidth',3.0)
end
xlim([0,x_lim_max_plot])
hold off