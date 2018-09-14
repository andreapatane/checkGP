function [Kstar,Kstarstar] = compute_kernel_on_test(test_data,training_data,kernel,depth,sigma_w_2,sigma_b_2,diags_K)

[Kstarstar, diags_Kstarstar] = getKernel(test_data,test_data,kernel,depth,sigma_w_2,sigma_b_2);
Kstar = getKernel(test_data,training_data,kernel,depth,sigma_w_2,sigma_b_2,diags_Kstarstar,diags_K);



end