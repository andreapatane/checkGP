function [trainedSystem,R,K_inv,gp_training_opts] = trainGP(training_data,training_labels,gp_training_opts)
%Implements training and inference for GP model


if strcmp(gp_training_opts.output_mode,'class_via_regress')
    one_hot_labels = one_hot_encoding(training_labels); %I do classification via one-hot encoding.
elseif strcmp(gp_training_opts.output_mode,'regression')
    one_hot_labels = training_labels;
end

%training
[K,gp_training_opts] = getKernel(training_data,training_data,gp_training_opts,true,training_labels);


%getting cholesky decomposition of K
while true
    try
        if isfield(gp_training_opts.kernel_params,'sigma_p') %if the error was estimated via MLE than I use it
            K = K + gp_training_opts.kernel_params.sigma_p*eye(size(K));
        else
            K = K + gp_training_opts.kernel_params.sig_eps*eye(size(K)); %otherwise I add a pre-specified va;ue
        end
        R = chol(K,'upper');
        break
    catch
        disp('Matrix K seems ill-conditioned, I am jittering it a bit more..')
    end
end
clear K

trainedSystem = R\(R'\one_hot_labels);
%getting the inverse of K as well for computation of the bounds on the
%variance.
K_inv = inv_chol(R);

end

function Y = inv_chol(R)
% Matrix Inversion using Cholesky Decomposition
%
% Finds the inverse of the matrix X, given its (lower triangular) Cholesky
% Decomposition; i.e. X = LL', according to the paper 'Matrix Inversion
% Using Cholesky Decomposition', Aravindh Krishnamoorthy, Deepak Menon,
% arXiv:1111.4144.
%

% Version 0.1, 2013-05-25, Aravindh Krishnamoorthy
% e-mail: aravindh.k@ieee.org

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(R, 1) ;
%Y = zeros(N, N) ;
% Construct the auxillary diagonal matrix S = 1/rii
Y = inv(diag(diag(R))) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=N:-1:1
    for i=j:-1:1
        Y(i,j) = Y(i,j) - R(i,i+1:end)*Y(i+1:end,j) ;
        Y(i,j) = Y(i,j)/R(i,i) ;
        % Write out the symmetric element
        Y(j,i) = Y(i,j);
    end
end

end