function [m_out,M_out] = max_min_of_coscos_plus_sinsin(m,M,x_L,x_U,y_L,y_U)
%finds minimum and maximum of cos(x)cos(y) + f(t) sin(x)sin(y)
%for x in [x_L,x_U], y in [y_L,y_U], given m and M maximum and minimum of
%f(t).

%following assumptions make code faster to run. Are based on the fact that
%input represents valid pixels
assert(x_L >= 0)
assert(x_U <= pi/2)
assert(y_L >= 0)
assert(y_U <= pi/2)


[m1,~,M1,~] = max_min_of_cos_plus_sin(cos(y_L),sin(y_L),m,M,x_L,x_U);
[m2,~,M2,~] = max_min_of_cos_plus_sin(cos(y_U),sin(y_U),m,M,x_L,x_U);

[m3,~,M3,~] = max_min_of_cos_plus_sin(cos(x_L),sin(x_L),m,M,y_L,y_U);
[m4,~,M4,~] = max_min_of_cos_plus_sin(cos(x_U),sin(x_U),m,M,y_L,y_U);

m_out = min([m1,m2,m3,m4]);
M_out = min([M1,M2,M3,M4]);
end