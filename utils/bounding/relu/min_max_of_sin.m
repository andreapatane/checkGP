function [m,M] = min_max_of_sin(m,M,x_L,x_U)
% find minimum and maximum of f(y) sin(x), given minimum and maximum values
% of f(y).
if (x_L >= 0) && (x_U<= pi/2)
    m = m*sin(x_L);
    M = M*sin(x_U);
else
    disp('Not too sure of what is happening here...')
    %first possible negative part for the cosine
    if x_L < 0
        current_upper = min(0,x_U);
        m1 = get_minimum(M,x_L,current_upper);
        M1 = get_maximum(m,x_L,current_upper);
        x_L = 0;
    else
        m1 = inf;
        M1 = -inf;
    end
    
    %first possible positive part for the cosine
    if (x_L < pi) && (x_L <= x_U)
        current_upper = min(pi,x_U);
        m2 = get_minimum(m,x_L,current_upper);
        M2 = get_maximum(M,x_L,current_upper);
        x_L = pi;
    else
        m2 = inf;
        M2 = -inf;
    end
    
    
    %second possible negative part for the cosine
    if (x_L < 2*pi) && (x_L <= x_U)
        current_upper = min((3/2)*pi,x_U);
        m3 = get_minimum(M,x_L,current_upper);
        M3 = get_maximum(m,x_L,current_upper);
        x_L = 2*pi;
    else
        m3 = inf;
        M3 = -inf;
    end
    
    %second possible positive part for the cosine
    if (x_L >= 2*pi) && (x_L <= x_U)
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
inps = [x_L,x_U,-pi/2,pi/2,(3/2)*pi];
inps = inps( (inps >= x_L) &  (inps <= x_U));

m = c*sin(inps);
m = min(m);
end

function M = get_maximum(c,x_L,x_U)
inps = [x_L,x_U,-pi/2,pi/2,(3/2)*pi];
inps = inps( (inps >= x_L) &  (inps <= x_U));
M = c*sin(inps);
M = max(M);
end