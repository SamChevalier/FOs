function Hess_R = MAP_HessianR(x)
% MAP_GRADIENTP: Build the MAP Regularization Hessian
%
% Inputs:
% 1) x            Vector of anonymous generator parameters
%
% Outputs:
% 1) Hess_R       Regularization Hessian
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hessian
n_vars = length(x);
Hess_R = sparse(n_vars,n_vars);

end

