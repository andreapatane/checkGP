function [K,diagK,sig_w_2,sig_b_2,sigma_p] = getKernel(X1,X2,kernel,Network_Depth,sig_w_2,sig_b_2,diags_K,diags_Kstarstar,mle_FLag,training_labels)
if nargin < 10
    training_labels = [];
    if nargin < 9
        mle_FLag = false;
        if nargin < 7
            diags_K = {};
            diags_Kstarstar = {};
        end
    end
end

sigma_p = [];
if nargout > 1
    diagK = cell(Network_Depth+1,1);
end

if strcmp(kernel,'sqe')
    
    
    if mle_FLag
        
        gprMdl = fitrgp(X1,training_labels,'KernelFunction','ardsquaredexponential');
        
        sig_w_2 = (gprMdl.KernelInformation.KernelParameters(end))^2;
        sig_b_2 = 1./(2*(gprMdl.KernelInformation.KernelParameters(1:(end-1))).^2);
        if any(sig_b_2 > 1)
            disp('Oh no! U need to update the implementation...')
        end

        K =  get_K_sqe(X1,X2,sig_w_2,sig_b_2);
        sigma_p = gprMdl.Sigma;
        disp(sig_w_2)
        disp(sig_b_2)
        disp(sigma_p)
    else
        
        K = get_K_sqe(X1,X2,sig_w_2,sig_b_2);
        
    end
elseif strcmp(kernel,'1N-ReLU')
    [~, d_in] = size(X1);
    K = zeros(size(X1,1),size(X2,1));
    
    c = sig_w_2 * (sig_b_2 + (sig_w_2/d_in) );
    c = c/(2*pi);
    for ii = 1:size(K,1)
        for jj = 1:size(K,2)
            beta = sig_b_2 + (sig_w_2/d_in) * X1(ii,:) * X2(jj,:)';
            beta = beta/(sig_b_2 + (sig_w_2/d_in));
            beta = max(min(beta,1),-1);
            gamma = acos(beta);
            K(ii,jj) = sig_b_2 + c * ( sin(gamma) + (pi - gamma) * cos(gamma)  );
            
        end
    end
    
elseif strcmp(kernel,'N-ReLU')
    [~, d_in] = size(X1);
    K = zeros(size(X1,1),size(X2,1));
    

    K_hat_prev = sig_b_2 + (sig_w_2/d_in);
    for c_l_i = 1:Network_Depth
        c = (sig_w_2/(2*pi) ) * K_hat_prev;
        for ii = 1:size(K,1)
            for jj = 1:size(K,2)
                if c_l_i == 1
                    beta = sig_b_2 + (sig_w_2/d_in) * X1(ii,:) * X2(jj,:)';
                else
                    beta =  K(ii,jj);
                end
                beta = beta/K_hat_prev;

                beta = max(min(beta,1),-1);
                gamma = acos(beta);
                K(ii,jj) = sig_b_2 + c * ( sin(gamma) + (pi - gamma) * cos(gamma)  );
                
            end
        end
        K_hat_prev = sig_b_2 + (sig_w_2/2)*K_hat_prev;
    end
    
end



end

function logLike = get_logLikelihood(Kdata,whiteData,noiseMatrix,X)


K = Kdata(abs(X)) + noiseMatrix;


logLike = whiteData' * ( K  \ whiteData ) + 2*sum(log(diag(chol(K)  )) )      ;



end

function K = get_K_sqe(X1,X2,sig_w_2,sig_b_2, flagSym)
if nargin < 5
    flagSym = false;
end
K = zeros(size(X1,1),size(X2,1));

if flagSym
    for ii = 1:size(K,1)
        for jj = 1:ii
            K(ii,jj) = sig_w_2 * exp( - dot(sig_b_2, (X1(ii,:) - X2(jj,:)).^2) )   ;
            K(jj,ii) = K(ii,jj);
        end
    end
else
    for ii = 1:size(K,1)
        for jj = 1:size(K,2)
            K(ii,jj) = sig_w_2 * exp( - dot(sig_b_2, (X1(ii,:) - X2(jj,:)).^2) )   ;
        end
    end
end
end






