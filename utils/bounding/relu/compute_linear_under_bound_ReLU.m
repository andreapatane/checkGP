function [a_i,b_i,a_lower,b_lower,a_upper,b_upper] = compute_linear_under_bound_ReLU(y_i,gamma_l,gamma_u,beta_l,beta_u,k1,k2,a_lower,b_lower,a_upper,b_upper)
%lower bound on ReLU kernel
if nargin < 11
    [a_lower, b_lower, a_upper, b_upper] = linear_bounds_on_K_of_gamma(gamma_l,gamma_u);
    
    
    if 0
        kerCore = @(x) (sin(x) + (pi-x).*cos(x));
        gamma_x = linspace(gamma_l,gamma_u,10);
        figure
        hold on
        plot(gamma_x,a_lower + b_lower*gamma_x)
        plot(gamma_x,a_upper + b_upper*gamma_x)
        plot(gamma_x,kerCore(gamma_x))
    end
    a_lower = k1 + k2 * a_lower;
    b_lower = k2 * b_lower;    
    a_upper = k1 + k2 * a_upper;
    b_upper = k2 * b_upper;
    if 0
        kerCore = @(x) (sin(x) + (pi-x).*cos(x));
        gamma_x = linspace(gamma_l,gamma_u,10);
        kernel_on_gamma = @(x,k1,k2) (k1+k2*kerCore(x));
        figure
        hold on
        plot(gamma_x,a_lower + b_lower*gamma_x)
        plot(gamma_x,a_upper + b_upper*gamma_x)
        plot(gamma_x,kernel_on_gamma(gamma_x,k1,k2))
    end
    
end

%taking into account the sign of y_i
if y_i>=0
    a_gamma = a_lower;
    b_gamma = b_lower;
else
    a_gamma = a_upper;
    b_gamma = b_upper;
end
%and multipling for it
a_i = a_gamma * y_i;
b_i = b_gamma * y_i;

if 0
    kerCore = @(x) (sin(x) + (pi-x).*cos(x));
    gamma_x = linspace(gamma_l,gamma_u,10);
    kernel_on_gamma = @(x,k1,k2) (y_i*(k1+k2*kerCore(x)));
    figure
    hold on
    plot(gamma_x,a_i + b_i*gamma_x)
    plot(gamma_x,kernel_on_gamma(gamma_x,k1,k2))
end

%getting linear bounds on acos
[a_lower_acos, b_lower_acos, a_upper_acos, b_upper_acos] = linear_bounds_on_acos_of_beta(beta_l,beta_u);
%taking into account b_gamma sign
if b_i>=0
    a_acos = a_lower_acos;
    b_acos = b_lower_acos;
else
    a_acos = a_upper_acos;
    b_acos = b_upper_acos;
end


a_i = a_i + b_i * a_acos;
b_i = b_i * b_acos;

if 0
    kerCore = @(x) (sin(x) + (pi-x).*cos(x));
    beta_x = linspace(beta_l,beta_u,10);
    kernel_on_beta = @(x,k1,k2) (y_i*(k1+k2*kerCore(acos(x))));
    figure
    hold on
    plot(beta_x,a_i + b_i*beta_x)
    plot(beta_x,kernel_on_beta(beta_x,k1,k2))
end



end