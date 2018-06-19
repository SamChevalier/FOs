function Hess_P = MAP_HessianP(x,var_params)
% MAP_GRADIENTP: Build the MAP Prior Hessian
%
% Inputs:
% 1) x            Vector of anonymous generator parameters
% 2) var_params   A Vector of generator parameter variance values
%
% Outputs:
% 1) Hess_P       Prior Hessian
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hessian Vector: Second Derivatives
Hess_Pvals = 2./var_params;

% Hessian
n_vars = length(x);
rcs    = 1:length(var_params);
Hess_P = sparse(rcs,rcs,Hess_Pvals,n_vars,n_vars);

end

