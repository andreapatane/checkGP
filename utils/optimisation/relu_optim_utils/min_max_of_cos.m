function [m,M] = min_max_of_cos(m,M,x_L,x_U)
% find minimum and maximum of f(y) cos(x), given minimum and maximum values
% of f(y).

if (x_L >= 0 && x_U <= pi/2)
    m = m*cos(x_U);
    M = M*cos(x_L);
else
    error('The code should never get here. In case check the code below, which should be correct anyway.')
    %first possible negative part for the cosine
    if x_L < -0.5*pi
        current_upper = min(-0.5*pi,x_U);
        m1 = get_minimum(M,x_L,current_upper);
        M1 = get_maximum(m,x_L,current_upper);
        x_L = -0.5*pi;
    else
        m1 = inf;
        M1 = -inf;
    end
    
    %first possible positive part for the cosine
    if (x_L < 0.5*pi) && (x_L <= x_U)
        current_upper = min(0.5*pi,x_U);
        m2 = get_minimum(m,x_L,current_upper);
        M2 = get_maximum(M,x_L,current_upper);
        x_L = 0.5*pi;
    else
        m2 = inf;
        M2 = -inf;
    end
    
    
    %second possible negative part for the cosine
    if (x_L < (3/2)*pi) && (x_L <= x_U)
        current_upper = min((3/2)*pi,x_U);
        m3 = get_minimum(M,x_L,current_upper);
        M3 = get_maximum(m,x_L,current_upper);
        x_L = (3/2)*pi;
    else
        m3 = inf;
        M3 = -inf;
    end
    
    %second possible positive part for the cosine
    if (x_L >= (3/2)*pi) && (x_L <= x_U)
        m4 = get_minimum(m,x_L,x_U);
        M4 = get_maximum(M,x_L,x_U);
    else
        m4 = inf;
        M4 = -inf;
    end
    m = min([m1,m2,m3,m4]);
    M = max([M1,M2,M3,M4]);
end
end


function m = get_minimum(c,x_L,x_U)
inps = [x_L,x_U,-pi,0,pi,2*pi];
inps = inps( (inps >= x_L) &  (inps <= x_U));

m = c*cos(inps);
m = min(m);
if m < eps
    m = 0;
end
end

function M = get_maximum(c,x_L,x_U)
inps = [x_L,x_U,-pi,0,pi,2*pi];
inps = inps( (inps >= x_L) &  (inps <= x_U));
M = c*cos(inps);
M = max(M);
if M < eps
    M = 0;
end
end