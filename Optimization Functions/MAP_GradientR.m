function Grad_R = MAP_GradientR(x,n_params,f_range,lambda)
% MAP_GRADIENTR: Build the MAP Regularization Gradient (never used during
%                stage 1, so no stage test)
%
% Inputs:
% 1) x            Vector of anonymous generator parameters + injections
% 2) n_params     How many generator parameters
% 3) f_range      Indices associated with forced oscillations: the
%                 frequencies f_vec(f_range) have current injection vars
% 4) lambda       Regularization penalty parameter
%
% Outputs:
% 1) Grad_R       Regularization Gradient
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The regularization penalty only consists of slack variables
n_vars  = n_params + 4*length(f_range);
n_slack = length(x) - n_vars;
Grad_R  = lambda*(ones(1,n_slack));

% Zero pad at the front
Grad_R = [zeros(1,n_vars) Grad_R];

end

