function [S1F,S2F] = MAP_LGS(X_params,Ys_s)
% MAP_LGS:         MAP Likelihood Gradient Structure              
%
% Inputs:
% 1) X_params      Vector of symbolic variables which are used to take the
%                  derivative of the functions
% 2) Ys_s          Symbolic admittance matrix
%
% Outputs:
% 1) S1F           Stage 1 function and gradient structure
% 2) S2F           Stage 2 function and gradient structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a new set of symbolic variables: Voltage/Current Data & Injections
syms Vmrs Vmis Vars Vais Imrs Imis Iars Iais 
syms Inj_mrs Inj_mis Inj_prs Inj_pis

% Evaluate YfC =)
Y11   = Ys_s(1,1);
Y12   = Ys_s(1,2);
Y21   = Ys_s(2,1);
Y22   = Ys_s(2,2);

% Split *admittance* into real and imag components
Y11r = real(Y11);
Y11i = imag(Y11);
Y12r = real(Y12);
Y12i = imag(Y12);
Y21r = real(Y21);
Y21i = imag(Y21);
Y22r = real(Y22);
Y22i = imag(Y22);

% ****** STAGE 1 ****** => Build Mr, Mi, Pr, and Pi
Mr = Imrs - Y11r.*Vmrs + Y11i.*Vmis - Y12r.*Vars + Y12i.*Vais;
Mi = Imis - Y11i.*Vmrs - Y11r.*Vmis - Y12i.*Vars - Y12r.*Vais;
Pr = Iars - Y21r.*Vmrs + Y21i.*Vmis - Y22r.*Vars + Y22i.*Vais;
Pi = Iais - Y21i.*Vmrs - Y21r.*Vmis - Y22i.*Vars - Y22r.*Vais;

% Take Derivatives
Grad_Mr = jacobian(Mr,X_params);
Grad_Mi = jacobian(Mi,X_params);
Grad_Pr = jacobian(Pr,X_params);
Grad_Pi = jacobian(Pi,X_params);

% Build Matlab equations
S1F.MrF      = matlabFunction(Mr,'vars',{'Omega_a',X_params,'Imrs','Vmrs','Vmis','Vars','Vais'});
S1F.MiF      = matlabFunction(Mi,'vars',{'Omega_a',X_params,'Imis','Vmrs','Vmis','Vars','Vais'});
S1F.PrF      = matlabFunction(Pr,'vars',{'Omega_a',X_params,'Iars','Vmrs','Vmis','Vars','Vais'});
S1F.PiF      = matlabFunction(Pi,'vars',{'Omega_a',X_params,'Iais','Vmrs','Vmis','Vars','Vais'});
S1F.Grad_MrF = matlabFunction(Grad_Mr,'vars',{'Omega_a',X_params,'Vmrs','Vmis','Vars','Vais'});
S1F.Grad_MiF = matlabFunction(Grad_Mi,'vars',{'Omega_a',X_params,'Vmrs','Vmis','Vars','Vais'});
S1F.Grad_PrF = matlabFunction(Grad_Pr,'vars',{'Omega_a',X_params,'Vmrs','Vmis','Vars','Vais'});
S1F.Grad_PiF = matlabFunction(Grad_Pi,'vars',{'Omega_a',X_params,'Vmrs','Vmis','Vars','Vais'});

% ****** STAGE 2 ****** => Reformulate
MrI = Imrs - Y11r.*Vmrs + Y11i.*Vmis - Y12r.*Vars + Y12i.*Vais - Inj_mrs;
MiI = Imis - Y11i.*Vmrs - Y11r.*Vmis - Y12i.*Vars - Y12r.*Vais - Inj_mis;
PrI = Iars - Y21r.*Vmrs + Y21i.*Vmis - Y22r.*Vars + Y22i.*Vais - Inj_prs;
PiI = Iais - Y21i.*Vmrs - Y21r.*Vmis - Y22i.*Vars - Y22r.*Vais - Inj_pis;

% Build Matlab functions
S2F.MrIF      = matlabFunction(MrI,'vars',{'Omega_a',X_params,'Imrs','Vmrs','Vmis','Vars','Vais','Inj_mrs'});
S2F.MiIF      = matlabFunction(MiI,'vars',{'Omega_a',X_params,'Imis','Vmrs','Vmis','Vars','Vais','Inj_mis'});
S2F.PrIF      = matlabFunction(PrI,'vars',{'Omega_a',X_params,'Iars','Vmrs','Vmis','Vars','Vais','Inj_prs'});
S2F.PiIF      = matlabFunction(PiI,'vars',{'Omega_a',X_params,'Iais','Vmrs','Vmis','Vars','Vais','Inj_pis'});
S2F.Grad_MrIF = S1F.Grad_MrF; % These are effectively the same, so no need to reformulate (df/dinj = -1 or 0)
S2F.Grad_MiIF = S1F.Grad_MiF;
S2F.Grad_PrIF = S1F.Grad_PrF;
S2F.Grad_PiIF = S1F.Grad_PiF;

end
