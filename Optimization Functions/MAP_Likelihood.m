function [Likelihood] = MAP_Likelihood(x,y_Vm,y_Va,y_Im,y_Ia,f_vec,YfC,Stage_1,f_range,Lk_Cov,n_params)
% MAP_LIKELIHOOD: Build the Likelihood Function
%
% Inputs:
% 1)  x            Vector of anonymous generator parameters, current
%                  injection variables, and slack variables
% 2)  y_Vm         Voltage magnitude fft data
% 3)  y_Va         Voltage phase fft data
% 4)  y_Im         Current magnitude fft data
% 5)  y_Ia         Current phase fft data
% 6)  f_vec        Vector of raw frequencies
% 7)  YfC          Admittance matrix function (cell matrix)
% 8)  Stage_1      Stage_1? If so, neglect current injection variables
% 9)  f_range      Indices associated with forced oscillations: the
%                  frequencies f_vec(f_range) have current injection vars
% 10) Lk_Cov       Numerical covariance matrix (sparse)
% 11) n_params     Unknown parameters in Y
%
% Outputs:
% 1) Likelihood   Likelihood function (actually, it is twice the negative 
%                 log of the true Gaussian Likelihood function)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate YfC =)
Omega = 2*pi*f_vec;
Y11   = YfC{1,1}(Omega,x(1:n_params));
Y12   = YfC{1,2}(Omega,x(1:n_params));
Y21   = YfC{2,1}(Omega,x(1:n_params));
Y22   = YfC{2,2}(Omega,x(1:n_params));

% Split *admittance* into real and imag components
Y11r = real(Y11);
Y11i = imag(Y11);
Y12r = real(Y12);
Y12i = imag(Y12);
Y21r = real(Y21);
Y21i = imag(Y21);
Y22r = real(Y22);
Y22i = imag(Y22);

% Split *data* into real and imag parts
Vmr = real(y_Vm);
Vmi = imag(y_Vm);
Var = real(y_Va);
Vai = imag(y_Va);
Imr = real(y_Im);
Imi = imag(y_Im);
Iar = real(y_Ia);
Iai = imag(y_Ia);

% Build Mr, Mi, Pr, and Pi
Mr = Imr - Y11r.*Vmr + Y11i.*Vmi - Y12r.*Var + Y12i.*Vai;
Mi = Imi - Y11i.*Vmr - Y11r.*Vmi - Y12i.*Var - Y12r.*Vai;
Pr = Iar - Y21r.*Vmr + Y21i.*Vmi - Y22r.*Var + Y22i.*Vai;
Pi = Iai - Y21i.*Vmr - Y21r.*Vmi - Y22i.*Var - Y22r.*Vai;

% Stage 1?
if Stage_1 == 1 
    % In this case, no current injection variables: 
    %  => Don't alter Mr, Mi, Pr, and Pi
    
else
    % In this case, current injection variables are added in f_range
    nf     = length(f_range);
    nf_rng = 1:nf;
    
    % Add in current injections
    Mr(f_range) = Mr(f_range) - x(n_params + 0*nf + nf_rng);
    Mi(f_range) = Mi(f_range) - x(n_params + 1*nf + nf_rng);
    Pr(f_range) = Pr(f_range) - x(n_params + 2*nf + nf_rng);
    Pi(f_range) = Pi(f_range) - x(n_params + 3*nf + nf_rng);
    
end

% Concatenate
L = [Mr.' Mi.' Pr.' Pi.'].';

% Now, build the Likelihood function
Likelihood = L.'*(Lk_Cov\L);

end

