function [x_train,y_train,x_test,y_test] =  x_product_data(gp_training_opts)
%Data generation for the running example. The script is adapted from the
%one provided in the Hierarchical kernel learning (HKL)
%https://www.di.ens.fr/~fbach/hkl/#matlab and described in "Francis Bach, Exploring Large Feature Spaces
%with Hierarchical Multiple Kernel Learning, NIPS 2008."

% fixing the seed of the random generators to the one showed in the paper
seed=1;
rng(seed)


ntrain = gp_training_opts.num_ts;
ntest = gp_training_opts.num_test_samples;


p = 1024;          
psub = 3;          
n = ntrain + ntest;            
s = 2;              % number of relevant variables
noise_std = 0.5;		% standard deviation of noise

% generate random covariance matrix from a Wishart distribution
Sigma_sqrt = randn(p,p);
Sigma = Sigma_sqrt' * Sigma_sqrt;


% normalize to unit trace and sample
diagonal = diag(Sigma);
Sigma_sqrt =   Sigma_sqrt * diag( 1./diagonal.^.5);
X = randn(n,p) * Sigma_sqrt;
X = X(:,1:psub);

rp = randperm(n);
trainset = rp(1:ntrain);
if ~gp_training_opts.synthetic_test_set
    testset  = rp(ntrain+1: (end -2) );
end

x_train = X(trainset,:);
[y_train,std_Y] = generate_Y_values(x_train,s,ntrain,noise_std);

if gp_training_opts.synthetic_test_set
    tt_x = linspace(gp_training_opts.extremes(1,1),gp_training_opts.extremes(1,2),sqrt(ntest));
    tt_y = linspace(gp_training_opts.extremes(2,1),gp_training_opts.extremes(2,2),sqrt(ntest));
    x_test = zeros(length(tt_x)*length(tt_y),3);
    for ii = 1:length(tt_x)
        for jj = 1:length(tt_y)
            x_test((ii - 1)*length(tt_y) + jj, 1  ) = tt_x(ii);
            x_test((ii - 1)*length(tt_y) + jj, 2  ) = tt_y(jj);
            x_test((ii - 1)*length(tt_y) + jj, 3  ) = X(randi(size(X,1)),3);
        end
    end
else
    x_0 = zeros(size(x_train(1,:)));
    x_c = 3*ones(size(x_train(1,:)));
    x_test = [x_0;x_c;X(testset,:)];
end

[y_test] = generate_Y_values(x_test,s,ntest,0,std_Y);


end
