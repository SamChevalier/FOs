function [sol_S1,sol_S2,y_data] = MAP_Opt(data,Y,STD_ns,freq_data,prior_data,Optns,S1F,S2F,S2_lambda)
% MAP_OPT:        Solve the MAP problem for a single generator
% 
% Inputs:
% 1) data         This the time series data associated with the voltage and
%                 current: data.Vm, data.Va, data.Im, & data.Ia
% 2) Y            Admittance matrix structure with three elements:
%                     1) Y.YfC      => Cell matrix of admittance matrix functions
%                     2) Y.Ys_s     => Symbolic admittance matrix
%                     3) Y.X_params => Vector of symbolic variables which are 
%                                      used to take derivatives of the functions
% 3) STD_ns       Structure of standard deviations of measurement noise:
%                     1) STD_ns.Vm
%                     2) STD_ns.Va
%                     3) STD_ns.Im
%                     4) STD_ns.Ia
% 4) freq_data    This structure contains the upper and lower bounds of the
%                 FO frequnecy range, and it contains the max frequency of
%                 concern: LowerFB, UpperFB, Max, & dt
% 5) prior_data   This is a structure which contains prior information:
%                     1) prior_data.prior_mean => Stage 1 initial conditions
%                     2) prior_data.prior_std  => Stage 1 prior std
%
% 6) Optns        Structure with OptTol, FunTol, MaxIts, StpTol & SolMtd
%                     1 => No Gradient, BFGS Hessian
%                     2 => Analytical Gradient, BFGS Hessian
%                     3 => Analytical Gradient, Analytical Hessian
%                 Also includes loop stopping criteria: op_its & op_SCP
% 7) S1F          Stage 1 function and gradient structure (hefty)
% 8) S2F          Stage 2 function and gradient structure (hefty)
% 9) S2_lambda    (Optional parameter) Stage 2 regularization parameter
% 
% Outputs:
% 1) sol_S1       History of Stage 1 Solutions
% 2) sol_S2       History of Stage 2 Solutions
% 3) y_data       The generator's fft data
%                 1) y_data.Vm
%                 2) y_data.Va
%                 3) y_data.Im
%                 4) y_data.Ia
%                 5) y_data.f       (*** Skip DC ***)
%                 6) y_data.f_range (range of FO's in y_data.f)
%                 7) y_data.prior_var
%                 8) y_data.std_w
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0) *** Pre-Stage Processing ***

% 0.1) Upack Data
Vm_PMUn = data.Vm;
Va_PMUn = data.Va;
Im_PMUn = data.Im;
Ia_PMUn = data.Ia;

% 0.2) Apply FFTs to raw data => yields raw (R) fft data 
N           = length(Vm_PMUn);
Hann_Wind   = hann(N);
[~,y_VmR]   = MAP_FFT(Hann_Wind.*detrend(Vm_PMUn,'constant'),freq_data.dt,N);
[~,y_VaR]   = MAP_FFT(Hann_Wind.*detrend(Va_PMUn,'constant'),freq_data.dt,N);
[~,y_ImR]   = MAP_FFT(Hann_Wind.*detrend(Im_PMUn,'constant'),freq_data.dt,N);
[fwR,y_IaR] = MAP_FFT(Hann_Wind.*detrend(Ia_PMUn,'constant'),freq_data.dt,N);
fwR         = fwR(:);

% 0.3) Immediately trim data: exclude DC and everything after the max frequency
[~,ind_MaxFreq] = min(abs(fwR - freq_data.Max));      % Maximum Frequnecy
ind_trm = 2:ind_MaxFreq;                              % *** Trim DC ***
y_Vm    = y_VmR(ind_trm);
y_Va    = y_VaR(ind_trm);
y_Im    = y_ImR(ind_trm);
y_Ia    = y_IaR(ind_trm);
fw      = fwR(ind_trm);

% 0.4) Determine standard deviation of the measurement noise at each frequnecy
std_w = MAP_MeasNoise(N,STD_ns);

% 0.5) Parse admittance matrix
X_params = Y.X_params;
YfC      = Y.YfC;

% 0.6) Unpack prior information
prior_mean = prior_data.prior_mean; % Initial conditions
prior_std  = prior_data.prior_std;  % Stage 1 prior std's (large)

%% 1) *** Stage 1 ***
% Stage 1 Decision Variable Structure: In this stage, x only consists of
% generator prameters, and the prior is extremely weak.
%

Stage    = 1;
n_params = length(X_params);

% Set the regularization parameter: not relevant in stage 1
lambda = 0;

% 1.1) Define data in f_range
f_range  = [];
n_FOs    = length(freq_data.LowerFB);                        % How many FOs?
for ii = 1:n_FOs
    [~,ind_LowerFB] = min(abs(fw - freq_data.LowerFB(ii)));  % Lower Injection Freq
    [~,ind_UpperFB] = min(abs(fw - freq_data.UpperFB(ii)));  % Upper Injection Freq
    f_range         = [f_range; (ind_LowerFB:ind_UpperFB)']; % FO index range
end

% Set the frequnecy range & remove frequencies in the FO range(s)
S1_range          = 1:length(fw);
S1_range(f_range) = [];

% 1.2) Define S1 data vectors
y_VmS1   = y_Vm(S1_range);
y_VaS1   = y_Va(S1_range);
y_ImS1   = y_Im(S1_range);
y_IaS1   = y_Ia(S1_range);
fw_S1    = fw(S1_range);

% 1.3) Define Prior
prior_var = (prior_mean.*prior_std).^2;
Prior     = @(x)MAP_Prior(x,prior_mean,prior_var);

% Initialize loop
PT        = 100; % Percent
fval      = 0;   % Cost
fval_prev = 0;   % Previous Cost
its       = 0;   % Iteration Count 
sol_S1    = prior_mean;

% Loop => Update Covariance Matrix
while (its <= Optns.op_its) && (PT > Optns.op_SCP)
    
    % Incremenet & Update
    its = its + 1;
    if its > 1
        fval_prev = fval;
    end
        
    % 1.4) Build Likelihood Covariance Matrix
    Lk_Cov = MAP_Lk_Cov(fw_S1,YfC,std_w,prior_mean);

    % 1.5) Likelihood Function
    Likelihood = @(x)MAP_Likelihood(x,y_VmS1,y_VaS1,y_ImS1,y_IaS1,fw_S1,YfC,Stage,f_range,Lk_Cov,n_params);
    
    % 1.6) Cost Function
    Cost = @(x)MAP_Cost(x,Likelihood,Prior,lambda,Stage,n_params,f_range);

    % 1.7) Build Gradients (Prior, Likelihood, and Regularization)
    Grad_P    = @(x)MAP_GradientP(x,prior_mean,prior_var);
    Grad_L    = @(x)MAP_GradientL(x,y_VmS1,y_VaS1,y_ImS1,y_IaS1,fw_S1,Lk_Cov,f_range,Stage,S1F,S2F,n_params);
    Grad_R    = @(x)MAP_GradientR(x,n_params,f_range,lambda);
    % Quick note: a Jacobian takes the form J = [df1/dx1...df1/dxn
    %                                                .        .
    %                                            dfn/dx1...dfn/dxn]
    % Matlab agrees, even for a scalar function f. Therefore, all gradients
    % will be defined as *row* vectors. Deal with it. Also, S1F & S2F are the
    % Stage 1 & 2 function and gradient structures

    % 1.8) Build Hessian
    Hess_P = @(x)MAP_HessianP(x,prior_var);
    Hess_L = @(x)MAP_HessianL(x,y_VmS1,y_VaS1,fw_S1,Lk_Cov,f_range,Stage,S1F,S2F,n_params);
    Hess_R = @(x)MAP_HessianR(x);
    
    % 1.9) Set up and Optimize
    x0 = prior_mean;
    A  = [];
    b  = [];
    [opt_sol,fval,~,Hess_fin] = MAP_nestedOpt(x0,A,b,Optns,Stage,Cost,Grad_P,Grad_L,Grad_R,Hess_P,Hess_L,Hess_R);
    
    % Print a message :)
    fprintf('After iteration %d, the cost is %d\n',[its fval]);
    
    % 1.10) Reset the prior mean, etc.
    sol_S1     = [sol_S1 opt_sol];
    prior_mean = opt_sol;
    PT         = 100*abs((fval - fval_prev)/fval_prev);
end

% Send message to user:
if (PT < Optns.op_SCP)
    fprintf('Stage 1 optimization loop terminated because percent change is %d.\n\n',PT);
end

%% 2) *** Stage 2 ***
% Stage 2 Decision Variable Structure
%
% What is the structure of x? Define the following local variables:
% n = # generator parameters
% m = # number of frequencies over which FO(s) is/are active
%
% x has length (n+2*4*m) => 2*(real and imag of current mag & phase)
% This is multiplied by 2 because of the slack variables
%
% x structure:
%   x = [n "gen params";
%        m "Ir_mag vars"
%        m "Ii_mag vars"
%        m "Ir_phs vars"
%        m "Ii_phs vars"
%        m "Ir_mag slacks"
%        m "Ii_mag slacks"
%        m "Ir_phs slacks"
%        m "Ii_phs slacks"
Stage = 2;

% 2.1) Update x0 and prior-mean
x0         = opt_sol;
prior_mean = x0;
x0         = [x0; zeros(8*length(f_range),1)];

% Update Stage 2 lambda (only if specified)
if nargin == 9
    lambda = S2_lambda;
else
    % Just set it equal to the value of the previous cost function
    lambda = fval;
end

% 2.2) Re-define the frequency domain data
y_VmS2   = y_Vm;
y_VaS2   = y_Va;
y_ImS2   = y_Im;
y_IaS2   = y_Ia;
fw_S2    = fw;

% 2.3) Update Prior with Stage 1 inverse Hessian
prior_var = diag(inv(Hess_fin));
Prior     = @(x)MAP_Prior(x,prior_mean,prior_var);

% 2.4) Build Likelihood Covariance Matrix
Lk_Cov = MAP_Lk_Cov(fw_S2,YfC,std_w,prior_mean);

% 2.5) Re-run Likelihood Function
Likelihood = @(x)MAP_Likelihood(x,y_VmS2,y_VaS2,y_ImS2,y_IaS2,fw_S2,YfC,Stage,f_range,Lk_Cov,n_params);

% 2.6) Re-compute Cost Function
Cost = @(x)MAP_Cost(x,Likelihood,Prior,lambda,Stage,n_params,f_range);

% 2.7) Re-build Gradients (Prior, Likelihood, and Regularization)
Grad_P    = @(x)MAP_GradientP(x,prior_mean,prior_var);
Grad_L    = @(x)MAP_GradientL(x,y_VmS2,y_VaS2,y_ImS2,y_IaS2,fw_S2,Lk_Cov,f_range,Stage,S1F,S2F,n_params);
Grad_R    = @(x)MAP_GradientR(x,n_params,f_range,lambda);

% 2.8) Re-build Hessians
Hess_P = @(x)MAP_HessianP(x,prior_var);
Hess_L = @(x)MAP_HessianL(x,y_VmS2,y_VaS2,fw_S2,Lk_Cov,f_range,Stage,S1F,S2F,n_params);
Hess_R = @(x)MAP_HessianR(x);

% 2.9) Build inequality constraints A and b associated with I injs & slacks:
% -s < I < s
n       = 4*length(f_range);
b       = sparse(2*n,1);
rws     = 1:2:2*n;
cls     = 1:n;
rws_all = [rws rws+1 rws rws+1];
cls_all = [cls cls cls+n cls+n];
n1v     = -1*ones(length(cls),1);
nv      = [n1v -n1v n1v n1v];
A       = sparse(rws_all,cls_all,nv,2*n,2*n);
A       = [sparse(2*n,n_params) A];

% 2.10) Optimize
Optns.MaxIts   = 100;
[~,~,sol_S2,~] = MAP_nestedOpt(x0,A,b,Optns,Stage,Cost,Grad_P,Grad_L,Grad_R,Hess_P,Hess_L,Hess_R);

%% 3) Define Output Data
% We want the output data to be defined from index = 2, to index = max
% freq. See above. Again, DC is *not* included!
y_data.Vm        = y_Vm;
y_data.Va        = y_Va;
y_data.Im        = y_Im;
y_data.Ia        = y_Ia;
y_data.f         = fw;
y_data.f_range   = f_range;
y_data.prior_var = prior_var;
y_data.std_w     = std_w;

end

