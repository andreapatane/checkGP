function [x_train,y_train,x_test,y_test,psub,n_classes,std_Y] =  x_product_data(ntrain,ntest,plotPosteriori)
% fixing the seed of the random generators
seed=1;
%randn('state',seed);
%rand('state',seed);
rng(seed)
% toy example characteristics
p = 1024;           % total number of variables (used to generate a Wishart distribution)
psub = 3;          % kept number of variables = dimension of the problem
n = ntrain + ntest;            % number of observations
s = 2;              % number of relevant variables
noise_std = 0.5;		% standard deviation of noise
n_classes = 1;

% generate random covariance matrix from a Wishart distribution
Sigma_sqrt = randn(p,p);
Sigma = Sigma_sqrt' * Sigma_sqrt;


% normalize to unit trace and sample
diagonal = diag(Sigma);
%Sigma = diag( 1./diagonal.^.5) * Sigma * diag( 1./diagonal.^.5);
Sigma_sqrt =   Sigma_sqrt * diag( 1./diagonal.^.5);
X = randn(n,p) * Sigma_sqrt;
X = X(:,1:psub);
p=psub;

rp = randperm(n);
trainset = rp(1:ntrain);
testset  = rp(ntrain+1:end);

x_train = X(trainset,:);
[y_train,std_Y] = generate_Y_values(x_train,s,ntrain,noise_std);

%x_test = X(testset,:);
if plotPosteriori
    tt = linspace(-3.5,3.5,sqrt(ntest));
    x_test = zeros(length(tt)*length(tt),3);
    for ii = 1:length(tt)
        for jj = 1:length(tt)
            x_test((ii - 1)*length(tt) + jj, 1  ) = tt(ii);
            x_test((ii - 1)*length(tt) + jj, 2  ) = tt(jj);
            x_test((ii - 1)*length(tt) + jj, 3  ) = X(randi(size(X,1)),3);
        end
    end
else
    x_0 = zeros(size(x_train(1,:)));
    x_c = 3*ones(size(x_train(1,:)));
    x_test = [x_0;x_c];
end

[y_test] = generate_Y_values(x_test,s,ntest,0,std_Y);





end
