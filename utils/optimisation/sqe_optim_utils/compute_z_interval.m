function [z_i_L, z_i_U] = compute_z_interval(x_i,x_L,x_U,theta_vec)

z_i_L = 0;

for jj = 1:length(x_i)
    
        if (x_i(jj) <= x_U(jj)) && (x_i(jj) >= x_L(jj))
            z_i_L_j = 0;
        else
            z_i_L_j = theta_vec(jj) * min( (x_i(jj) - x_L(jj))^2  , (x_i(jj) - x_U(jj))^2  );
        end       
        
    z_i_L = z_i_L + z_i_L_j;
    
end



z_i_U = dot(theta_vec,max([(x_i - x_L).^2;(x_i - x_U).^2]));



end