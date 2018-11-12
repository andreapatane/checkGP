function [a,b] = merge_cnc_cvx_lower(a_cnc,b_cnc,a_cvx,b_cvx,x_l,x_u,cvxFirstFlag)
%merging the lower bounds from the concave and convex region together in a
%single linear lowerbound


if cvxFirstFlag
    f_l = a_cvx + b_cvx * x_l;
    
    f1_u = a_cnc + b_cnc * x_u;
    f2_u = a_cvx + b_cvx * x_u;
    f_u = min(f1_u,f2_u);
else
    f1_l = a_cnc + b_cnc * x_l;
    f2_l = a_cvx + b_cvx * x_l;

    f_l = min(f1_l,f2_l);
    f_u = a_cvx + b_cvx * x_u;
end
[a,b] = line_thru_points(x_l,f_l,x_u,f_u);

end