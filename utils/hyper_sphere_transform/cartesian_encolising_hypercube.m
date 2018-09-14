function [x_L,x_U] = cartesian_encolising_hypercube(phi_L,phi_U)
% find an hypercube in cartesian coordiante that encloses the hypercube
% [phi_L,phi_U] given in angular coordinate.

m_prev = 1;
M_prev = 1;
x_L = zeros(length(phi_L)+1,1);
x_U = zeros(size(x_L));
for ii = 1:(length(x_L) - 1)
    [x_L(ii),x_U(ii)] = min_max_of_cos(m_prev,M_prev,phi_L(ii),phi_U(ii));
    [m_prev,M_prev] = min_max_of_sin(m_prev,M_prev,phi_L(ii),phi_U(ii));
end
x_L(end) = m_prev;
x_U(end) = M_prev;


end