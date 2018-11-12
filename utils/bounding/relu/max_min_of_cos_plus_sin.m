function [m_out,x_m_out,M_out,x_M_out] = max_min_of_cos_plus_sin(a,b,m,M,x_L,x_U)
%compute maximum and minimum of a*cos(x) + b*g(y)*sin(x), given m and M
%maximum and minimum of g(y)


%fix for b<0: I invert the minimum and the maximum of g(y).
if b < 0
    temp = m;
    m = M;
    M = temp;
end

%x_m_c = comp_critic_point_cos_plus_sin(a,b,m);
%x_M_c = comp_critic_point_cos_plus_sin(a,b,M);

x_m_c = atan((m*b)/a);
x_M_c = atan((M*b)/a);

x_0_L = max(x_L,0);
x_pi_U = min(x_U,pi);


%computations for sin x >= 0
%[m_out,x_m_out] =  get_minimum(a,b,m,x_0_L,x_pi_U,x_m_c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shamelessly inlining get_minimum for performance reasons...
m_out = a*cos(x_0_L) + b*m*sin(x_0_L);
x_m_out = x_0_L;
val2 = a*cos(x_pi_U) + b*m*sin(x_pi_U);
if val2 < m_out
    m_out = val2;
    x_m_out = x_pi_U;
end

if (x_0_L < x_m_c) && (x_m_c < x_pi_U)
    val3 =  a*cos(x_m_c) + b*m*sin(x_m_c);
    if val3 < m_out
        m_out = val3;
        x_m_out = x_m_c;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[M_out,x_M_out] =  get_maximum(a,b,M,x_0_L,x_pi_U,x_M_c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shamelessly inlining get_maximum for performance reasons...
M_out = a*cos(x_0_L) + b*M*sin(x_0_L);
x_M_out = x_0_L;
val2 = a*cos(x_pi_U) + b*M*sin(x_pi_U);
if val2 > M_out
    M_out = val2;
    x_M_out = x_pi_U;
end

if (x_0_L < x_M_c) && (x_M_c < x_pi_U)
    val3 =  a*cos(x_M_c) + b*M*sin(x_M_c);
    if val3 > M_out
        M_out = val3;
        x_M_out = x_M_c;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%computations for the left part of the inverval (if any)
if x_L < x_0_L
    [m_2,x_m_2] =  get_minimum(a,b,M,x_L,0,x_M_c);
    [M_2,x_M_2] =  get_maximum(a,b,m,x_L,0,x_m_c);
    if m_2 < m_out
        m_out = m_2;
        x_m_out = x_m_2;
    end
    if M_2 > M_out
        M_out = M_2;
        x_M_out = x_M_2;
    end
    % else
    %     m_2 = inf;
    %     M_2 = -inf;
    %     x_m_2 = nan;
    %     x_M_2 = nan;
end


if x_U < x_0_L
    [m_3,x_m_3] =  get_minimum(a,b,M,pi,x_U,x_M_c);
    [M_3,x_M_3] =  get_maximum(a,b,m,pi,x_U,x_m_c);
    if m_3 < m_out
        m_out = m_3;
        x_m_out = x_m_3;
    end
    if M_3 > M_out
        M_out = M_3;
        x_M_out = x_M_3;
    end
    % else
    %     m_3 = inf;
    %     M_3 = -inf;
    %     x_m_3 = nan;
    %     x_M_3 = nan;
    
end

% [m_out,m_ind] = min([m_1,m_2,m_3]);
% xms = [x_m_1,x_m_2,x_m_3];
% x_m_out = xms(m_ind);
%
% [M_out,M_ind] = max([M_1,M_2,M_3]);
% xMs = [x_M_1,x_M_2,x_M_3];
% x_M_out = xMs(M_ind);

end

function [m, x_m] = get_minimum(a,b,c,x_L,x_U,x_c)


m = a*cos(x_L) + b*c*sin(x_L);
x_m = x_L;

val2 = a*cos(x_U) + b*c*sin(x_U);
if val2 < m
    m = val2;
    x_m = x_U;
end

if (x_L < x_c) && (x_c < x_U)
    val3 =  a*cos(x_c) + b*c*sin(x_c);
    if val3 < m
        m = val3;
        x_m = x_c;
    end
% else
%     val3 = -inf;
end

% [m,ind] = min([val1,val2,val3]);
% xs = [x_L,x_U,x_c];
% x_m = xs(ind);
end

function [m, x_m] = get_maximum(a,b,c,x_L,x_U,x_c)

m = a*cos(x_L) + b*c*sin(x_L);
x_m = x_L;

val2 = a*cos(x_U) + b*c*sin(x_U);
if val2 > m
    m = val2;
    x_m = x_U;
end

if (x_L < x_c) && (x_c < x_U)
    val3 =  a*cos(x_c) + b*c*sin(x_c);
    if val3 > m
        m = val3;
        x_m = x_c;
    end
% else
%     val3 = -inf;
end

% [m,ind] = max([val1,val2,val3]);
% xs = [x_L,x_U,x_c];
% x_m = xs(ind);
end
