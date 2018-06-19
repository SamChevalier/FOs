function [opt_sol,fval,hist_sol,Hess_fin] = MAP_nestedOpt(x0,A,b,Optns,Stage_1,Cost,Grad_P,Grad_L,Grad_R,Hess_P,Hess_L,Hess_R)
% MAP_NESTEDOPT: Optimization with Gradient + Hessian
%
% Inputs:
% 1)  x0          Variable Input (these are iterated)
% 2)  A           Linear inequality constraint matrix
% 3)  b           Linear inequality constraint vector
% 4)  Optns       Structure with OptTol, FunTol, MaxIts, StpTol & SolMtd
%                    1 => No Gradient, BFGS Hessian
%                    2 => Analytical Gradient, BFGS Hessian
%                    3 => Analytical Gradient, Analytical Hessian
% 5)  Stage_1     Stage 1?
% 6)  Cost        Formal cost function
% 7)  Grad_P      Prior Gradient
% 8)  Grad_L      Likelihood Gradient
% 9)  Grad_R      Regularization Gradient
% 10) Hess_P      Prior Hessian
% 11) Hess_L      Likelihood Hessian
% 12) Hess_R      Regularization Hessian
%
% Outputs:
% 1) opt_sol      Optimal solution
% 2) fval         Final cost
% 3) hist_sol     Iterative solution history
% 4) Hess_fin     Final Hessian
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial Definitions
hist_sol     = [];
options = optimoptions('fmincon','Algorithm','interior-point','Display','off',...
    'OptimalityTolerance',Optns.OptTol,'MaxIterations',Optns.MaxIts,...
    'OutputFcn',@myoutput,'FunctionTolerance',Optns.FunTol,'StepTolerance',Optns.StpTol);

% Solution Method
switch Optns.SolMtd
    case 1
        options.SpecifyObjectiveGradient = false;
        options.HessianFcn               = [];
    case 2
        options.SpecifyObjectiveGradient = true;
        options.HessianFcn               = [];
    case 3
        options.SpecifyObjectiveGradient = true;
        if Stage_1 == 1 
            options.HessianFcn = @hessinteriorS1;
        else
            options.HessianFcn = @hessinteriorS2;
        end
end

% Stage?
if Stage_1 == 1
    opt_fun = @nestedfunS1;
    n       = length(x0);
    lb      = zeros(n,1);    % Lower bound of 0 for all parameters
    options.Display = 'off'; % Custom messages given
else
    opt_fun = @nestedfunS2;
    n       = length(x0)-length(b);
    lb      = [zeros(n,1); -Inf*ones(length(b),1)]; % Lower bound of 0 for all parameters
    options.Display = 'iter';                       % Custom messages given
end

% Optimize
[opt_sol,fval,~,~,~,~,Hess_fin] = fmincon(opt_fun,x0,A,b,[],[],lb,[],[],options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% %%%%%%%%%% Sub Functions %%%%%%%%%% %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sub Function: Save iteration data
    function stop = myoutput(x,~,state)
        stop = false;
        if isequal(state,'iter')
            hist_sol   = [hist_sol x];
        end
    end

% Computes S1 objective function & gradient
    function [c,g] = nestedfunS1(x0)
        % Cost Function & Gradient
        c = Cost(x0);
        g = Grad_P(x0) + Grad_L(x0);
    end

% Computes S2 objective function & gradient
    function [c,g] = nestedfunS2(x0)
        c = Cost(x0);
        g = Grad_P(x0) + Grad_L(x0) + Grad_R(x0);
    end

% Computes S1 Hessian
    function [h] = hessinteriorS1(x0,~)
        h = Hess_P(x0) + Hess_L(x0);
    end

% Computes S2 Hessian
    function [h] = hessinteriorS2(x0,~)
        h = Hess_P(x0) + Hess_L(x0) + Hess_R(x0);
    end

end