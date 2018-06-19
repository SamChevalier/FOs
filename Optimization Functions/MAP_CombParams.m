function [Param_Matrix] = MAP_CombParams(p1,p2,p3)
% MAP_COMBPARAMS:   Build a matrix of parameter values based on combinations
%                   of the vectors p1, p2, and p3
%
% Inputs:
% 1) p1             Vector of low parameter values
% 2) p2             Vector of mid parameter values
% 3) p3             Vector of high parameter values
%
% Outputs:
% 1) Param_Matrix   Matrix of parameter combinations: each column is a full
%                   set of parameter values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_samps  = 3; % Hardcoded
n_params = length(p1);
rows     = n_params;
cols     = n_samps^n_params;

% Shift Mat & Index vetor
SM      = [p1 p2 p3];
idx_vec = ones(rows,1);
% Perform successive cycling
Param_Matrix = zeros(rows,cols);
it_idx       = 1; 
while it_idx <= cols
    for ii = 1:n_params
        Param_Matrix(ii,it_idx) = SM(ii,idx_vec(ii));
    end
    % Increment properaly
    idx_vec(1) = idx_vec(1) + 1;
    if idx_vec(1) > n_samps
        % Use this as a que to take action
        idx_vec(1) = 1;
        % If n_params > 1, enter loop
        if n_params > 1
            for ii = 2:n_params
                idx_vec(ii) = idx_vec(ii) + 1;
                if idx_vec(ii) > n_samps
                    % Continue
                    idx_vec(ii) = 1;
                else
                    break
                end
            end
        end
    end
    % Incremenet index
    it_idx = it_idx + 1;
end
end



