function [f_val_min, f_val_max] = simult_separate_quadprog_Hg0(H,f,x_L,x_U)



f_val_min = 0;
f_val_max = 0;
for jj = 1:length(x_L)
    if x_L(jj) < x_U(jj)
        x_critic = -f(jj);
        if (x_critic <= x_U(jj)) && (x_critic >= x_L(jj))
            f_val_min_par = 0;
            f_val_max_par = max(x_L(jj)^2 + f(jj)*x_L(jj),x_U(jj)^2 + f(jj)*x_U(jj));
        else
            fl = x_L(jj)^2 + f(jj)*x_L(jj);
            fu = x_U(jj)^2 + f(jj)*x_U(jj);
            if fl <= fu
                f_val_min_par = fl;
                f_val_max_par = fu;
            else
                f_val_min_par = fu;
                f_val_max_par = fl;
            end
        end      
    else
        f_val_min_par = x_L(jj)^2 + f(jj)*x_L(jj);
        f_val_max_par = f_val_min_par;
    end
    f_val_min = f_val_min + H(jj)*f_val_min_par;
    f_val_max = f_val_max + H(jj)*f_val_max_par;
end

        
end
