function [phi_L,phi_U] = compute_angle_space_hyper_rectangle(x_L,x_U)
% I over-approximate the projection of a cartesian hyper-rectangle in
% spherical coordinates. As a result I get an hyper-rectangle in the
% spherical space such that:
%
% [phi_L,phi_U] contains {phi(x) | x \in [x_L,x_U]}

phi_L = zeros(1,length(x_L)-1);
phi_U = zeros(size(phi_L));


%precomputing maximum and minimum of single component squares
min_x_i_square = zeros(size(x_L));
max_x_i_square = zeros(size(x_L));
for ii = 1:length(x_L)
    x_L_square = x_L(ii)^2;
    x_U_square = x_U(ii)^2;
    if (x_L(ii) < 0) && (x_U(ii) > 0 )
        min_x_i_square(ii) = 0;
    else
        min_x_i_square(ii) = min(x_L_square,x_U_square);
    end
    max_x_i_square(ii) = max(x_L_square,x_U_square);

end

%I add some jitter to the last values, to avoid complicating the code too
%much...
if min_x_i_square(end) == 0 
    min_x_i_square(end) = min_x_i_square(end)+eps;
end
if max_x_i_square(end) == 0 
    max_x_i_square(end) = max_x_i_square(end)+eps;
end

%computing actual extremes of hyper-rectangle in spherical coordinates
for ii = 1:length(phi_L)
    a_L = sum(min_x_i_square((ii+1):end));
    a_U = sum(max_x_i_square((ii+1):end));

    phi_L(ii) = acos(   x_U(ii) / sqrt(a_L + x_U(ii))   );
    phi_U(ii) = acos(   x_L(ii) / sqrt(a_U + x_L(ii))   );

end




end