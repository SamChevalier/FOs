function [Ys] = MAP_Build_Y(F_vec,G_vec,X,Uv)
% MAP_BUILD_Y: Use Function vectors and variable vectors to build an
%              analytical function for Y
%
% Inputs:
% 1) F_vec     Differential Equations
% 2) G_vec     Algebraic Equations (for Im and Ia)
% 3) X         State Variables
% 4) Uv        Algebraic Variables (V and Theta)
%
% Output:
% 1) Ys        Symbolic MATLAB function for Yf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms Omega_a
j = sqrt(-1);

% Build Jacobians
JacFX  = jacobian(F_vec,X);
JacFUv = jacobian(F_vec,Uv);
JacGX  = jacobian(G_vec,X);
JacGUv = jacobian(G_vec,Uv);

% Admittance Matrix Function
Ys = (JacGX/(j*Omega_a*eye(length(X)) - JacFX))*JacFUv + JacGUv;

end

