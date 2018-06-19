function Grad_L = MAP_GradientL(x,y_Vm,y_Va,y_Im,y_Ia,f_vec,Lk_Cov,f_range,Stage_1,S1F,S2F,n_params)
% MAP_GRADIENTL: This function builds the gradient associated with the
%                likelihood function. The gradient element associated with 
%                parameter p has the form (dC/dp) = 2(dx'/dp)*G*x. 
%                Therefore, given the parameters in x, the gradient may 
%                be computed.
%
% Inputs:
% 1)  x            Vector of anonymous generator parameters
% 2)  y_Vm         Voltage magnitude fft data
% 3)  y_Va         Voltage phase fft data
% 4)  y_Im         Current magnitude fft data
% 5)  y_Ia         Current phase fft data
% 6)  f_vec        Vector of raw frequencies
% 7)  Lk_Cov       Numerical covariance matrix (sparse)
% 8)  f_range      Indices associated with forced oscillations: the
%                  frequencies f_vec(f_range) have current injection vars
% 9)  Stage_1      Stage_1? If so, neglect current injection variables
% 10) SF1          Stage 1 functions and gradients
% 11) SF2          Stage 2 functions and gradients
% 12) n_params     Number of generator parameters
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
Imr = real(y_Im);
Imi = imag(y_Im);
Iar = real(y_Ia);
Iai = imag(y_Ia);

% Angular Frequnecy
Omega = 2*pi*f_vec;

% Current injection Functions
if Stage_1 == 1
    
    % Evaluate the functions
    MrFn = S1F.MrF(Omega,x,Imr,Vmr,Vmi,Var,Vai);
    MiFn = S1F.MiF(Omega,x,Imi,Vmr,Vmi,Var,Vai);
    PrFn = S1F.PrF(Omega,x,Iar,Vmr,Vmi,Var,Vai);
    PiFn = S1F.PiF(Omega,x,Iai,Vmr,Vmi,Var,Vai);
    
    % In thise case, just evaluate the gradients across all frequnecies
    Grad_MrFn = S1F.Grad_MrF(Omega,x,Vmr,Vmi,Var,Vai);
    Grad_MiFn = S1F.Grad_MiF(Omega,x,Vmr,Vmi,Var,Vai);
    Grad_PrFn = S1F.Grad_PrF(Omega,x,Vmr,Vmi,Var,Vai);
    Grad_PiFn = S1F.Grad_PiF(Omega,x,Vmr,Vmi,Var,Vai);
    
    % Assemble Gradient
    Grad_L = zeros(1,length(x));
    FV     = [MrFn; MiFn; PrFn; PiFn];
    for ii = 1:n_params
        DV         = [Grad_MrFn(:,ii); Grad_MiFn(:,ii); Grad_PrFn(:,ii); Grad_PiFn(:,ii)];
        Grad_L(ii) = 2*(DV.')*(Lk_Cov\FV);
    end
    
else
    
    % Split x into generator parameters and Injection variables
    nf         = length(f_range);
    nf_rng     = 1:nf;
    gen_params = x(1:n_params);
    Inj_mr     = x(n_params + 0*nf + nf_rng);
    Inj_mi     = x(n_params + 1*nf + nf_rng);
    Inj_pr     = x(n_params + 2*nf + nf_rng);
    Inj_pi     = x(n_params + 3*nf + nf_rng);
    
    % Now, split the frequnecies into two subsets
    idF  = f_range;
    idnF = 1:length(f_vec);
    idnF(f_range) = [];
    
    % Evaluate the functions: non-FO
    MrFn = S1F.MrF(Omega(idnF),gen_params,Imr(idnF),Vmr(idnF),Vmi(idnF),Var(idnF),Vai(idnF));
    MiFn = S1F.MiF(Omega(idnF),gen_params,Imi(idnF),Vmr(idnF),Vmi(idnF),Var(idnF),Vai(idnF));
    PrFn = S1F.PrF(Omega(idnF),gen_params,Iar(idnF),Vmr(idnF),Vmi(idnF),Var(idnF),Vai(idnF));
    PiFn = S1F.PiF(Omega(idnF),gen_params,Iai(idnF),Vmr(idnF),Vmi(idnF),Var(idnF),Vai(idnF));
    
    % Evaluate the functions: FO
    MrIFn = S2F.MrIF(Omega(idF),gen_params,Imr(idF),Vmr(idF),Vmi(idF),Var(idF),Vai(idF),Inj_mr);
    MiIFn = S2F.MiIF(Omega(idF),gen_params,Imi(idF),Vmr(idF),Vmi(idF),Var(idF),Vai(idF),Inj_mi);
    PrIFn = S2F.PrIF(Omega(idF),gen_params,Iar(idF),Vmr(idF),Vmi(idF),Var(idF),Vai(idF),Inj_pr);
    PiIFn = S2F.PiIF(Omega(idF),gen_params,Iai(idF),Vmr(idF),Vmi(idF),Var(idF),Vai(idF),Inj_pi);
    
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
    
    % Function Assembly (static)
    FV = [MrFn(id1); MrIFn; MrFn((id1(end)+1):end);
          MiFn(id1); MiIFn; MiFn((id1(end)+1):end);
          PrFn(id1); PrIFn; PrFn((id1(end)+1):end);
          PiFn(id1); PiIFn; PiFn((id1(end)+1):end)];

    % Assemble Gradient: Nothing to say about the slack variables
    Grad_L = zeros(1,length(x));
    for ii = 1:(n_params+4*nf)
        % Derivative Assembly
        DV = [Grad_MrFn(id1,ii); Grad_MrIFn(:,ii); Grad_MrFn((id1(end)+1):end,ii);
              Grad_MiFn(id1,ii); Grad_MiIFn(:,ii); Grad_MiFn((id1(end)+1):end,ii);
              Grad_PrFn(id1,ii); Grad_PrIFn(:,ii); Grad_PrFn((id1(end)+1):end,ii);
              Grad_PiFn(id1,ii); Grad_PiIFn(:,ii); Grad_PiFn((id1(end)+1):end,ii)];
        
        % Gradient
        Grad_L(ii) = 2*(DV.')*(Lk_Cov\FV);
    end
end

end
