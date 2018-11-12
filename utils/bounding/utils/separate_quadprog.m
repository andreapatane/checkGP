function [x_star, f_val] = separate_quadprog(H,f,x_L,x_U)


x_star = x_L';
f_val = 0;
for jj = 1:length(x_L)
    if x_L(jj) < x_U(jj)
        if H(jj) >= 0
            x_critic = -f(jj)/H(jj);
            if (x_critic <= x_U(jj)) && (x_critic >= x_L(jj))
                x_star(jj) = x_critic;
                f_val_par = 0.5*H(jj)*x_critic^2 + f(jj)*x_critic;
            else
                [f_val_par, idx] = min([0.5*H(jj)*x_L(jj)^2 + f(jj)*x_L(jj), 0.5*H(jj)*x_U(jj)^2 + f(jj)*x_U(jj)   ]) ;
                if idx == 2 
                    x_star(jj) = x_U(jj);
                end
                
            end
        else
            [f_val_par, idx] = min([0.5*H(jj)*x_L(jj)^2 + f(jj)*x_L(jj), 0.5*H(jj)*x_U(jj)^2 + f(jj)*x_U(jj)   ]) ;
            if idx == 2
                x_star(jj) = x_U(jj);
            end
        end
    else
        f_val_par = 0.5*H(jj)*x_L(jj)^2 + f(jj)*x_L(jj);
    end
     f_val = f_val + f_val_par;
end




end