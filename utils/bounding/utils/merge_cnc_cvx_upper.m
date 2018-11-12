function [a,b] = merge_cnc_cvx_upper(a_cnc,b_cnc,a_cvx,b_cvx,x_l,x_u,cvxFirstFlag)
%merging the upper bounds from the concave and convex region together in a
%single linear upperbound


if cvxFirstFlag
    g_u = a_cnc + b_cnc * x_u;
    
    g1_l = a_cnc + b_cnc * x_l;
    g2_l = a_cvx + b_cvx * x_l;
    g_l = max(g1_l,g2_l);
else
    g_l =  a_cnc + b_cnc * x_l;
   
    g1_u = a_cnc + b_cnc * x_u;
    g2_u = a_cvx + b_cvx * x_u;
    
    g_u = max(g1_u,g2_u);
end

[a,b] = line_thru_points(x_l,g_l,x_u,g_u);

end




