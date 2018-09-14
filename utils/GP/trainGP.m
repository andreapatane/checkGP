function [trainedSystem,R,diags_K,K_inv,sigma_w_2,sigma_b_2] = trainGP(training_data,training_labels,sigma_w_2,sigma_b_2,depth,kernel,sig_eps,cache_K_inv,dataSet,mle_FLag)


% getting normalised one-hot encoding
if strcmp(dataSet,'mnist')
    one_hot_labels = one_hot_encoding(training_labels);
elseif strcmp(dataSet,'x_products')
    one_hot_labels = training_labels;
end
% computing Covariance matrix K
[K,diags_K,sigma_w_2,sigma_b_2,sigma_p] = getKernel(training_data,training_data,kernel,depth,sigma_w_2,sigma_b_2,{},{},mle_FLag,training_labels);


%getting cholesky decomposition of K
while true
    try
        if ~isempty(sigma_p)
            K = K + sigma_p*eye(size(K));
        else
            K = K + sig_eps*eye(size(K));
        end
        R = chol(K,'upper');
        break
    catch
        disp('Matrix K seems ill-conditioned, I am jittering it a bit more..')
    end
end
clear K
trainedSystem = R\(R'\one_hot_labels);

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