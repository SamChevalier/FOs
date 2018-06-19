function Prior = MAP_Prior(x,GP_vec_P,var_params)
% MAP_PRIOR: Build a MATLAB function which is the portion of the cost
%            function associated with the Gaussian Prior
%
% Inputs:
% 1) x            Vector of anonymous generator parameters
% 2) GP_vec_P     Vector of generator parameter mean values
% 3) var_params   This is a vector of true parameter variances
%
% Outputs:
% 1) Prior        Prior function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Prior
Prior = sum((((x-GP_vec_P).^2)./var_params));

end

