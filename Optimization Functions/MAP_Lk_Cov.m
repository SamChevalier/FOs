function [Lk_Cov] = MAP_Lk_Cov(f_vec,YfC,std_w,GP_vec_P)
% MAP_LL_COV:      Build the Likelihood covariance matrix using
%                  probabilistic methods based on sampling from normal
%                  distributions associated with the parameters.
%
% Inputs:
% 1) f_vec        Vector of raw frequencies
% 2) YfC          Admittance matrix: cell array of functions
% 3) std_w        Frequency domain noise strength: structure
% 4) GP_vec_P     Mean prior values (for evaluating Yf)
%
% Outputs:
% 1) Lk_Cov        Constant numerical covariance matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define sparse Covariance matrix
n      = length(f_vec);
rng    = 1:n;
Lk_Cov = sparse(4*n,4*n);

% Evaluate the Admittance Entries
Omega = 2*pi*f_vec;
Yn11  = YfC{1,1}(Omega,GP_vec_P);
Yn12  = YfC{1,2}(Omega,GP_vec_P);
Yn21  = YfC{2,1}(Omega,GP_vec_P);
Yn22  = YfC{2,2}(Omega,GP_vec_P);

% AWGN variance across frequnecies (real and imaginary)
var_w.Vm = std_w.Vm^2;
var_w.Va = std_w.Va^2;
var_w.Im = std_w.Im^2;
var_w.Ia = std_w.Ia^2;

% Build the Covariance Matrix's diagonal matrices
Nr = real(Yn11).^2*var_w.Vm + imag(Yn11).^2*var_w.Vm + real(Yn12).^2*var_w.Va + imag(Yn12).^2*var_w.Va + var_w.Im;
Ni = imag(Yn11).^2*var_w.Vm + real(Yn11).^2*var_w.Vm + imag(Yn12).^2*var_w.Va + real(Yn12).^2*var_w.Va + var_w.Im;
Qr = real(Yn21).^2*var_w.Vm + imag(Yn21).^2*var_w.Vm + real(Yn22).^2*var_w.Va + imag(Yn22).^2*var_w.Va + var_w.Ia;
Qi = imag(Yn21).^2*var_w.Vm + real(Yn21).^2*var_w.Vm + imag(Yn22).^2*var_w.Va + real(Yn22).^2*var_w.Va + var_w.Ia;

% Update Covariance Matrix
Lk_Cov(0*n+rng,0*n+rng) = diag(sparse(Nr));
Lk_Cov(1*n+rng,1*n+rng) = diag(sparse(Ni));
Lk_Cov(2*n+rng,2*n+rng) = diag(sparse(Qr));
Lk_Cov(3*n+rng,3*n+rng) = diag(sparse(Qi));

% Build the Covariance Matrix's off-diagonal matrices: Start with Zero Mats
NiNr = zeros(n,1);
QiQr = zeros(n,1);

% Update Covariance Matrix
Lk_Cov(0*n+rng,1*n+rng) = diag(sparse(NiNr));
Lk_Cov(1*n+rng,0*n+rng) = diag(sparse(NiNr));

Lk_Cov(2*n+rng,3*n+rng) = diag(sparse(QiQr));
Lk_Cov(3*n+rng,2*n+rng) = diag(sparse(QiQr));

% Build the Covariance Matrix's off-diagonal matrices: Non-Zero Mats
QrNr = real(Yn11)*var_w.Vm.*real(Yn21) + imag(Yn11)*var_w.Vm.*imag(Yn21) + real(Yn12)*var_w.Va.*real(Yn22) + imag(Yn12)*var_w.Va.*imag(Yn22);
QrNi = imag(Yn11)*var_w.Vm.*real(Yn21) - real(Yn11)*var_w.Vm.*imag(Yn21) + imag(Yn12)*var_w.Va.*real(Yn22) - real(Yn12)*var_w.Va.*imag(Yn22);
QiNr = real(Yn11)*var_w.Vm.*imag(Yn21) - imag(Yn11)*var_w.Vm.*real(Yn21) + real(Yn12)*var_w.Va.*imag(Yn22) - imag(Yn12)*var_w.Va.*real(Yn22);
QiNi = imag(Yn11)*var_w.Vm.*imag(Yn21) + real(Yn11)*var_w.Vm.*real(Yn21) + imag(Yn12)*var_w.Va.*imag(Yn22) + real(Yn12)*var_w.Va.*real(Yn22);

% Update Covariance Matrix
Lk_Cov(0*n+rng,2*n+rng) = diag(sparse(QrNr));
Lk_Cov(2*n+rng,0*n+rng) = diag(sparse(QrNr));

Lk_Cov(1*n+rng,2*n+rng) = diag(sparse(QrNi));
Lk_Cov(2*n+rng,1*n+rng) = diag(sparse(QrNi));

Lk_Cov(0*n+rng,3*n+rng) = diag(sparse(QiNr));
Lk_Cov(3*n+rng,0*n+rng) = diag(sparse(QiNr));

Lk_Cov(1*n+rng,3*n+rng) = diag(sparse(QiNi));
Lk_Cov(3*n+rng,1*n+rng) = diag(sparse(QiNi));

end

