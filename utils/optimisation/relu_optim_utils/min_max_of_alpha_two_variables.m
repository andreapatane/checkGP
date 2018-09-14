function [alpha_L,alpha_U] = min_max_of_alpha_two_variables(phi_L,phi_U)
%found minimum and maximum of the angle between two vectors of an hypersphere using angular coordiantes:
% max (min) acos(x_1 \dot x_2)
% subject to ||x_1||_2 = 1 , ||x_2||_2 = 1
%     phi_j^L   <=  phi_1j <= phi_j^U
%     phi_j^L   <=  phi_2j <= phi_j^U
% where phi_1j is the jth angle of x_1 in spherical coordinates
% and phi_2j is the jth angle of x_2 in spherical coordinates.
%%%%
%
% 





m_prev = 1;
M_prev = 1;
m = length(phi_L);
for ii = 0:(m-1)
   [m_prev,M_prev] = max_min_of_coscos_plus_sinsin(m_prev,M_prev,phi_L(m-ii),phi_U(m-ii),phi_L(m-ii),phi_U(m-ii));
end

alpha_U = acos(m_prev);
alpha_L = acos(M_prev);




end