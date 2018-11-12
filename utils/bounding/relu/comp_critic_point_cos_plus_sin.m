function x_c = comp_critic_point_cos_plus_sin(a,b,c)
% compute the critic point of a*cos(x) + b*c*sin(x)

x_c =  atan((c*b)/a);

end