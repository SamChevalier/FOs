function Cost = MAP_Cost(x,Likelihood,Prior,lambda,Stage_1,n_params,f_range)
% MAP_PRIOR: Build the cost function which shall be fed to the optimizer
%
% Inputs:
% 1) x            Vector of anonymous generator parameters, injections, and
%                 slack variables which constrain the injections
% 2) Likelihood   Likelihood function
% 3) Prior        Prior function: scaled sum of squares
% 4) Lambda       Lasso regularization parameter (injection penalty)
% 5) Stage_1      Stage_1? If so, neglect current injection variables
% 6) n_params     Number of generator parameters
% 7) f_range      Indices associated with forced oscillations: the
%                 frequencies f_vec(f_range) have current injection vars
%
% Outputs:
% 1) Cost         Cost function to be given to the optimizer
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure column vector
x = x(:);

% Stage 1?
if Stage_1 == 1 
    % In this case, no current injection variables: 
    gen_params = x(1:n_params);                        % Generator parameters
    Cost       = Prior(gen_params) + Likelihood(gen_params);
    
else
    % In this case, consider current injection variables:
    nf            = length(f_range);
    gen_params    = x(1:n_params);                        % Generator parameters
    gen_params_Is = x(1:(n_params+4*nf));                 % Generator parameters + current injections
    slacks        = x((n_params+4*nf+1):(n_params+8*nf)); % Slack variables
    Cost          = Prior(gen_params) + Likelihood(gen_params_Is) + lambda*sum(slacks);
end

end

