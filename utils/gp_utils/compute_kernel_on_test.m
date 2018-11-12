function [Kstar,Kstarstar] = compute_kernel_on_test(test_data,training_data,gp_training_opts)

[Kstarstar] = getKernel(test_data,test_data,gp_training_opts);
Kstar = getKernel(test_data,training_data,gp_training_opts);



end