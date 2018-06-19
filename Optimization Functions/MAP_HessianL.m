function Hess_R = MAP_HessianL(x,y_Vm,y_Va,f_vec,Lk_Cov,f_range,Stage_1,S1F,S2F,n_params)
% MAP_HessianL:  This function builds the approximate Hessian associated
%                with the likelihood function.
%
% Inputs:
% 1)  x            Vector of anonymous generator parameters
% 2)  y_Vm         Voltage magnitude fft data
% 3)  y_Va         Voltage phase fft data
% 4)  f_vec        Vector of raw frequencies
% 5)  Lk_Cov       Numerical covariance matrix (sparse)
% 6)  f_range      Indices associated with forced oscillations: the
%                  frequencies f_vec(f_range) have current injection vars
% 7)  Stage_1      Stage_1? If so, neglect current injection variables
% 8)  SF1          Stage 1 functions and gradients
% 9)  SF2          Stage 2 functions and gradients
% 10) n_params     Number of generator parameters
%
% Outputs:
% 1) Grad_L        Gradient of the Likelihood Function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Split *data* into real and imag parts
Vmr = real(y_Vm);
Vmi = imag(y_Vm);
Var = real(y_Va);
Vai = imag(y_Va);

% Angular Frequnecy
Omega = 2*pi*f_vec;

% Current injection Functions
if Stage_1 == 1
    
    % In thise case, just evaluate the gradients across all frequnecies
    Grad_MrFn = S1F.Grad_MrF(Omega,x,Vmr,Vmi,Var,Vai);
    Grad_MiFn = S1F.Grad_MiF(Omega,x,Vmr,Vmi,Var,Vai);
    Grad_PrFn = S1F.Grad_PrF(Omega,x,Vmr,Vmi,Var,Vai);
    Grad_PiFn = S1F.Grad_PiF(Omega,x,Vmr,Vmi,Var,Vai);
    
    % Gradient Vector
    GV = [Grad_MrFn; Grad_MiFn; Grad_PrFn; Grad_PiFn];
    
    % Assemble Approximate Hessian: This Hessian has, effectively,
    % neglected 2nd order derivatives entirely, and is only using gradient
    % data for Hessian assembly.
    Hess_R = 2*(GV.')*(Lk_Cov\GV);

else
    
    % Split x into generator parameters and Injection variables
    nf         = length(f_range);
    gen_params = x(1:n_params);
    
    % Now, split the frequnecies into two subsets
    idF  = f_range;            % FO Freqs
    idnF = 1:length(f_vec);
    idnF(f_range) = [];        % non-FO Freqs
    
    % Evaluate the non-FO Gradients
    Grad_MrFn = S1F.Grad_MrF(Omega(idnF),gen_params,Vmr(idnF),Vmi(idnF),Var(idnF),Vai(idnF));
    Grad_MiFn = S1F.Grad_MiF(Omega(idnF),gen_params,Vmr(idnF),Vmi(idnF),Var(idnF),Vai(idnF));
    Grad_PrFn = S1F.Grad_PrF(Omega(idnF),gen_params,Vmr(idnF),Vmi(idnF),Var(idnF),Vai(idnF));
    Grad_PiFn = S1F.Grad_PiF(Omega(idnF),gen_params,Vmr(idnF),Vmi(idnF),Var(idnF),Vai(idnF));
    
    % Add in the derivatives associated with the Injections
    nfv       = length(f_vec) - nf;  % Number of freqs with no injections
    zm        = zeros(nfv,nf);
    Grad_MrFn = [Grad_MrFn zm zm zm zm];
    Grad_MiFn = [Grad_MiFn zm zm zm zm];
    Grad_PrFn = [Grad_PrFn zm zm zm zm];
    Grad_PiFn = [Grad_PiFn zm zm zm zm];
    
    % Evaluate the FO Gradients
    Grad_MrIFn = S2F.Grad_MrIF(Omega(idF),gen_params,Vmr(idF),Vmi(idF),Var(idF),Vai(idF));
    Grad_MiIFn = S2F.Grad_MiIF(Omega(idF),gen_params,Vmr(idF),Vmi(idF),Var(idF),Vai(idF));
    Grad_PrIFn = S2F.Grad_PrIF(Omega(idF),gen_params,Vmr(idF),Vmi(idF),Var(idF),Vai(idF));
    Grad_PiIFn = S2F.Grad_PiIF(Omega(idF),gen_params,Vmr(idF),Vmi(idF),Var(idF),Vai(idF));
    
    % Add in the derivatives associated with the Injections
    zm         = zeros(nf,nf);
    n1m        = -1*eye(nf,nf);
    Grad_MrIFn = [Grad_MrIFn n1m zm zm zm];
    Grad_MiIFn = [Grad_MiIFn zm n1m zm zm];
    Grad_PrIFn = [Grad_PrIFn zm zm n1m zm];
    Grad_PiIFn = [Grad_PiIFn zm zm zm n1m];
    
    % Assemble Gradient Indices
    id1 = 1:(f_range(1)-1);
    
    % Assemble Gradient
    Grad_MrFn0 = [Grad_MrFn(id1,:); Grad_MrIFn; Grad_MrFn((id1(end)+1):end,:)];
    Grad_MiFn0 = [Grad_MiFn(id1,:); Grad_MiIFn; Grad_MiFn((id1(end)+1):end,:)];
    Grad_PrFn0 = [Grad_PrFn(id1,:); Grad_PrIFn; Grad_PrFn((id1(end)+1):end,:)];
    Grad_PiFn0 = [Grad_PiFn(id1,:); Grad_PiIFn; Grad_PiFn((id1(end)+1):end,:)];

    % Gradient Vector
    GV = [Grad_MrFn0; Grad_MiFn0; Grad_PrFn0; Grad_PiFn0];
    
    % Assemble Approximate Hessian
    n_vars = length(x);
    Hess_R = sparse(n_vars,n_vars);
    n_IPs  = n_vars - 4*nf;
    
    % Assemble Approximate Hessian: This Hessian has, effectively,
    % neglected 2nd order derivatives entirely, and is only using gradient
    % data for Hessian assembly.
    Hess_R(1:n_IPs,1:n_IPs) = 2*(GV.')*(Lk_Cov\GV);

end
end
