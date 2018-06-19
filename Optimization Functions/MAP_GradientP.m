function Grad_P = MAP_GradientP(x,GP_vec_P,var_params)
% MAP_GRADIENTP: Build the MAP Prior Gradient
%
% Inputs:
% 1) x            Vector of anonymous generator parameters
% 2) GP_vec_P     Vector of generator parameter mean values
% 3) var_params   A Vector of generator parameter mean values
%
% Outputs:
% 1) Prior        Prior function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Grad_P             = zeros(1,length(x));
n_params           = length(GP_vec_P);
gen_params         = x(1:n_params);
Grad_P(1:n_params) = ((2*(gen_params-GP_vec_P))./var_params);

% Rotate to row vector
Grad_P = Grad_P(:)';

end

