function [m_prev,phi_m,M_prev,phi_M] = min_max_of_dot_product_hypersphere(psi,phi_L,phi_U)
%found minimum and maximum of the dot product in the surface of an hypersphere using angular coordiantes:
% max (min) x_1 \dot x_2
% subject to ||x_1||_2 = 1
%     phi_j^L   <=  phi_j <= phi_j^U
% where phi_j is the jth angle of x_1 in spherical coordinates.
%%%%
%
% 

% 
% 
% [~,idxs] = max([psi - phi_L; phi_U - psi]);
% phi_m = phi_L;
% phi_m(logical(idxs - 1)) = phi_U(logical(idxs - 1));
% phi_M = psi;
% M_prev = 1;
% m_prev = dot(get_cartesian_coordinates(phi_m,[]),get_cartesian_coordinates(psi,[]));



phi_m = zeros(size(psi));
phi_M = zeros(size(psi));
m_prev = 1;
M_prev = 1;
m = length(psi);
for ii = 0:(m-1)
   [m_prev,phi_m(m-ii),M_prev,phi_M(m-ii)] = max_min_of_cos_plus_sin( cos(psi(m-ii)) , ...
                                             sin(psi(m-ii)) , m_prev,M_prev,phi_L(m-ii),phi_U(m-ii));
end






end