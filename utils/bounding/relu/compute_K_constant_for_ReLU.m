function [K,alpha_L,alpha_U] =  compute_K_constant_for_ReLU(phi_L,phi_U,sigma_b_2,sigma_w_2,m,depth)

[alpha_L,alpha_U] = min_max_of_alpha_two_variables(phi_L,phi_U);

c1 = sigma_b_2;
c2 = sigma_w_2/m;
c0 = (sigma_w_2^2)/(pi*m);

d_fun_1 = @(x)(c0 * sin(x) .*  (pi - acos( (c1+c2*cos(x))/(c1+c2)   )    )      );

k1 = sigma_b_2 + (sigma_w_2/2) * (sigma_b_2 + (sigma_w_2/m));

fun_1 = @(x)( c1 + c0*(  sin( acos( (c1+c2*cos(x))/(c1+c2)  )    ) + (pi - acos( (c1+c2*cos(x))/(c1+c2) )  ).*((c1+c2*cos(x))/(c1+c2))              )         );   

%c02 = (sigma_w_2/(2*pi) )*k1;
d_fun = @(x,c) (c * (pi - acos(x)));


fun = @(c,x)( c1 + c*(  sin( acos( x  )    ) + (pi - acos( x )  ).*x              )         );   



k_diag = sigma_b_2 + sigma_w_2/m;

dfun2OptCell = {};
fun2OptCell = {};
for ii = 1:depth
    if ii == 1
        dfun2OptCell{ii} = @(x)(d_fun_1(x));
        fun2OptCell{ii} = @(x)(fun_1(x));
    else
        c = k_diag*sigma_w_2/(2*pi);
        fun2OptCell{ii} = @(x)(fun(c,  1/(k_diag).*fun2OptCell{ii-1}(x)   ));
        dfun2OptCell{ii} = @(x)(d_fun( fun2OptCell{ii-1}(x)    ,c) .* (1/k_diag).*dfun2OptCell{ii-1}(x) );
    end
    k_diag = sigma_b_2 + (sigma_w_2/2) * k_diag;
end


%if depth == 1
%    fun2Opt = @(x)(abs(d_fun_1(x)));
%elseif depth == 2
%    fun2Opt = @(x)(  abs(d_fun((1/k1)*fun1(x)   ) .*(1/k1).* d_fun_1(x))   );
%else
%    error('not implemented yet...')
%end

fun2Opt = @(x)(abs(dfun2OptCell{end}(x)));
K1 = fun2Opt(alpha_U);
[~,K2] = fminbnd(@(x)(-fun2Opt(x)),alpha_L,alpha_U);
K2 = -K2;
K = max(K1,K2);


end




