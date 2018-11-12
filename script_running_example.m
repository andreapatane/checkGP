%%%%%%%%%%
%Script for analysis and plotting of Figure 2 in the paper.
%%%%%%%%%%%%%

clear all
close all
clc
addpath(genpath('utils/'))

%setting up gp options
gp_training_opts.kernel = 'sqe'; %squared exponential kernel
gp_training_opts.mle = true; %using mle for hyper-parameters estimation
gp_training_opts.dataSet = 'x_products'; %string name for data set used
gp_training_opts.num_ts = 128; %num of training samples
gp_training_opts.num_test_samples = 2; %num of test samples [just using x^o = (0,0,0) and x^*=(3,3,3) defined in the paper]
gp_training_opts.output_mode = 'regression'; %regression task
gp_training_opts.synthetic_test_set = false; %setting grid test points flag to false


result_folder = 'results/';
nsamples = 100; %number of samples for property empirical approximation


%performing gp training
[test_data,test_norm,test_labels,mu_test,...
    Kstar,training_data,trainedSystem,K_inv,R,gp_training_opts] = main_trainGP(gp_training_opts);


%specifiying analysis for both x^o and x^* (that is index 1 and 2 of the test set)
startIdx = 1; 
endIdx = 2;
%setting options for computation of \hat{\phi}
upper_bound_comp_opts.epsilons = 0.1; %setting radius of hyper-cube
upper_bound_comp_opts.invariance = false; %true: invariance property, false: safety property

upper_bound_comp_opts.delta = 0:0.0001:1.0; %delta values for property evaluation
upper_bound_comp_opts.pix_2_mod = {[1,2]}; %performing modification on the two significant value of the problem
upper_bound_comp_opts.kernel = gp_training_opts.kernel; %updating kernel parameters 
upper_bound_comp_opts.kernel_params = gp_training_opts.kernel_params; %updating kernel parameters 

upper_bound_comp_opts.maxIterations_mu = 30000; %maximum number of iterations to be performed in the bnb for \mu computation
upper_bound_comp_opts.maxIterations_xi = 30000; %maximum number of iterations to be performed in the bnb for \xi computation


%upper_bound_comp_opts.tollerance_mu = 0.000001; %tollerace value used as stopping criteria for \mu bnb. Increase to make processing faster
%upper_bound_comp_opts.tollerance_xi = 0.00003; %tollerace value used as stopping criteria for \xi bnb. Increase to make processing faster

upper_bound_comp_opts.tollerance_mu = 0.001;
upper_bound_comp_opts.tollerance_xi = 0.003;

upper_bound_comp_opts.constrain_2_one = false; %if true, than allowed input perturbation are constrained to belong in [0,1]
upper_bound_comp_opts.output_mode = gp_training_opts.output_mode;
nameTempl = 'sqe_x_prod_xo_and_xc'; %just an identifier for the results

num_of_threads = 1; %num of threads used by quadprog if needed

%%%
%Finally performing the actual analysis.
[p_ub,sup_mu,sup_xi,pixel2Modify_cell] = main_phi_hat_computation(upper_bound_comp_opts,...
    num_of_threads,startIdx,endIdx,...
    training_data,test_data,test_norm,test_labels,mu_test,Kstar,trainedSystem,K_inv,R);

%saving results of the analysis to file
saveResults(p_ub,sup_mu,sup_xi,pixel2Modify_cell,upper_bound_comp_opts,startIdx,endIdx,nameTempl,'formal',result_folder);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$Getting empirical approximation%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Origin test point
testPoint = [0,0];

gp_training_opts.extremes = [testPoint(1) - upper_bound_comp_opts.epsilons,testPoint(1) + upper_bound_comp_opts.epsilons;
    testPoint(2) - upper_bound_comp_opts.epsilons,testPoint(2) + upper_bound_comp_opts.epsilons];

gp_training_opts.mle = false;
gp_training_opts.num_ts = 128;
gp_training_opts.num_test_samples = 2025;
gp_training_opts.synthetic_test_set = true;


[test_data,~,~,mu_test,...
    ~,~,~,~,~,gp_training_opts,sigma_test] = main_trainGP(gp_training_opts);

stat_probs_x_o = main_sampling_approximation(nsamples,test_data,sigma_test,gp_training_opts,...
    upper_bound_comp_opts.delta,upper_bound_comp_opts.invariance,mu_test);

saveResults(stat_probs_x_o,sup_mu,sup_xi,pixel2Modify_cell{1},upper_bound_comp_opts,1,1,'sqe_x_prod_xo','approximated',result_folder)

%Corner test point
testPoint = [3,3];
gp_training_opts.extremes = [testPoint(1) - upper_bound_comp_opts.epsilons,testPoint(1) + upper_bound_comp_opts.epsilons;
    testPoint(2) - upper_bound_comp_opts.epsilons,testPoint(2) + upper_bound_comp_opts.epsilons];

[test_data,~,~,mu_test,...
    ~,~,~,~,~,gp_training_opts,sigma_test] = main_trainGP(gp_training_opts);
stat_probs_x_c = main_sampling_approximation(nsamples,test_data,sigma_test,gp_training_opts,...
    upper_bound_comp_opts.delta,upper_bound_comp_opts.invariance,mu_test);

saveResults(stat_probs_x_c,[],[],pixel2Modify_cell{2},upper_bound_comp_opts,2,2,'sqe_x_prod_xc','approximated',result_folder)


script_plot_example