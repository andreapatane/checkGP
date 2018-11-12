function [K,gp_training_opts] = getKernel(X1,X2,gp_training_opts,trainFlag,training_labels)
if nargin < 5
    training_labels = [];
    if nargin < 4
        trainFlag = false;
    end
end



if strcmp(gp_training_opts.kernel,'sqe')
    
    
    if gp_training_opts.mle && trainFlag
        
        %To perform MLE on sqe hyper-parameters I rely on the GP matlab
        %toolbox implementation.
        disp('Going for MLE of hyper-parameters...')
        gprMdl = fitrgp(X1,training_labels,'KernelFunction','ardsquaredexponential');
        
        %updating the GP option structure with new hyper-parameters extimation
        gp_training_opts.kernel_params.sigma = (gprMdl.KernelInformation.KernelParameters(end))^2;
        gp_training_opts.kernel_params.theta_vec = 1./(2*(gprMdl.KernelInformation.KernelParameters(1:(end-1))).^2);
        gp_training_opts.kernel_params.sigma_p = gprMdl.Sigma;
        if any(gp_training_opts.kernel_params.theta_vec > 1)
            % This case is not currently implemented in the code, as it was
            % not needed yet.
            % Briefly if any theta > 1, than some inequalities used for the
            % computation of the bound, do not hold true anymore. This can
            % be solved by simply rescaling the input features so that
            % every theta <= 1
            error('Case not currently accounted for by the code. Still need to update it...')
        end
        K =  get_K_sqe(X1,X2,gp_training_opts.kernel_params.sigma,gp_training_opts.kernel_params.theta_vec);
    else
        K =  get_K_sqe(X1,X2,gp_training_opts.kernel_params.sigma,gp_training_opts.kernel_params.theta_vec);
    end
elseif strcmp(gp_training_opts.kernel,'ReLU')
    %implements ReLU kernel. For performance reasone the code assumes that
    %every input vector lies in the unit sphere. That is ||x|| = 1, for
    %every x in the training and test set. This is guaranteed by first sending the input set into the to_norm1 function 
    %(ideally this should have Already done in some parent function of this). The code cannot be used if this
    %assumption is not true.
    
    [~, d_in] = size(X1);
    K = zeros(size(X1,1),size(X2,1));
    
    
    
    K_hat_prev = gp_training_opts.kernel_params.sigma_b_2 + (gp_training_opts.kernel_params.sigma_w_2/d_in);
    for c_l_i = 1:gp_training_opts.kernel_params.depth
        c = (gp_training_opts.kernel_params.sigma_w_2/(2*pi) ) * K_hat_prev;
        for ii = 1:size(K,1)
            for jj = 1:size(K,2)
                if c_l_i == 1
                    beta = gp_training_opts.kernel_params.sigma_b_2 + (gp_training_opts.kernel_params.sigma_w_2/d_in) * X1(ii,:) * X2(jj,:)';
                else
                    beta =  K(ii,jj);
                end
                beta = beta/K_hat_prev;
                if (beta > 1.1) || (beta < - 1.1) %checking for numerical errors...
                    error('Numerical errors are building up too quickly. Not sure what is going wrong...')
                end
                beta = max(min(beta,1),-1);
                gamma = acos(beta);
                K(ii,jj) = gp_training_opts.kernel_params.sigma_b_2 + c * ( sin(gamma) + (pi - gamma) * cos(gamma)  );
                
            end
        end
        K_hat_prev = gp_training_opts.kernel_params.sigma_b_2 + (gp_training_opts.kernel_params.sigma_w_2/2)*K_hat_prev;
    end
    
end



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






