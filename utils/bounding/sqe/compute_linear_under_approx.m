function [a_i,b_i] = compute_linear_under_approx(y_i,z_i_L,z_i_U)
%getting linear under approximation of sqe kernel

if y_i >= 0
    z_i_M = 0.5*(z_i_L + z_i_U);
    a_i = (1 + z_i_M)*y_i*exp(-z_i_M);
    b_i = -y_i*exp(-z_i_M); 

else
    a_i = y_i * ( exp(-z_i_L)  - z_i_L * (  exp(-z_i_L) -  exp(-z_i_U) )/( z_i_L - z_i_U   )  );
    b_i = y_i * (  exp(-z_i_L) -  exp(-z_i_U) )/( z_i_L - z_i_U   );
end



end