%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Script for Figure 1 in the paper. Train a GP regression model using a sqe
%kernel on a quadratic regression problem. It then plots a-posteriori mean
%and variance in a grid over [-3,3]x[-3,3] along with location of training
%samples used.
%%%%%%%%%%%%%


clear all
close all
clc
addpath(genpath('utils/'))


%Setting up gp training options thru a matlab structure
gp_training_opts.kernel = 'sqe'; %selecting squared-exponential kernel
gp_training_opts.mle = true; %setting MLE for hyper-parameters selection (otherwise you need to specify hyper-parameters explicitely)
gp_training_opts.dataSet = 'x_products'; %x_products dataset to be used (This is just a quadratic polynomial of two variables)
gp_training_opts.num_ts = 128; % number of training samples
gp_training_opts.num_test_samples = 625; %number of test samples
gp_training_opts.synthetic_test_set = true; %option used to generate a uniform testing grid for the plot
gp_training_opts.extremes = [-3,3;-3,3]; %extremes of the plotting grid (see just above)



gp_training_opts.output_mode = 'regression'; %GP is used for a standard regression problem


%%%%%
%%Performing Training of the GP
[~,~,~,mu_test,...
    ~,training_data,~,~,~,gp_training_opts,sigma_test] = main_trainGP(gp_training_opts);
%%%%

%%%%
%%Producing plots
tt_x = linspace(gp_training_opts.extremes(1,1),gp_training_opts.extremes(1,2),sqrt(gp_training_opts.num_test_samples));
tt_y = linspace(gp_training_opts.extremes(2,1),gp_training_opts.extremes(2,2),sqrt(gp_training_opts.num_test_samples));

mu_4_plot = zeros(length(tt_x),length(tt_y));
var_4_plot = zeros(length(tt_x),length(tt_y));


for ii = 1:length(tt_x)
    for jj = 1:length(tt_y)
        mu_4_plot(ii,jj) = mu_test(  (ii-1)*length(tt_y) + jj  );
    end
end

for ii = 1:length(tt_x)
    for jj = 1:length(tt_y)
        var_4_plot(ii,jj) = sigma_test(  (ii-1)*length(tt_y) + jj, (ii-1)*length(tt_y) + jj );
    end
end

figure
hold on
h = surf(tt_x,tt_y,mu_4_plot);
xlim([gp_training_opts.extremes(1,1),gp_training_opts.extremes(1,2)])
ylim([gp_training_opts.extremes(2,1),gp_training_opts.extremes(2,2)])
colormap(coolwarm(64))
colorbar
scatter3(training_data(:,1),training_data(:,2),ones(size(training_data,1),1),'k','x','LineWidth',1.5)
set(gca,'FontSize',16)
xlabel('x_1')
ylabel('x_2')
hold off

figure
hold on
surf(tt_x,tt_y,var_4_plot);
xlim([gp_training_opts.extremes(1,1),gp_training_opts.extremes(1,2)])
ylim([gp_training_opts.extremes(2,1),gp_training_opts.extremes(2,2)])
colormap(coolwarm(64))
scatter3(training_data(:,1),training_data(:,2),ones(size(training_data,1),1),'k','x','LineWidth',1.5)
colorbar
xlabel('x_1')
ylabel('x_2')
set(gca,'FontSize',16)

hold off
