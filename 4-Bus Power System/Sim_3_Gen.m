%% 4 Bus System Simulation: Inf Bus and 3 3rd Order + AVR Gens
%  Bus 1: Gen
%  Bus 2: Gen (FO)
%  Bus 3: Gen
%  Bus 4: Infinite Bus
%
% Gen1 --- Gen2 --- Gen3 --- Infinte_Bus ->
%
%  Use PF to model power flows
clear variables;

%% EE's rule
j = sqrt(-1);

% Simulation Length
dt       = 0.05;
tf       = 120;

% Oscillation (Percent Torque) & Perturbation Data
Osc_Amp  = 0.05;
Osc_Frq  = 0.5;
Rnd_Amp  = 0.01;

%% Steady State Power Flow: PV, PV, PV and Swing
%  Voltage
IC.V1 = 1;
IC.V2 = 1;
IC.V3 = 1;
IC.V4 = 1;

% Power: Electrical and Mechanical (Stator Resistance Ignored)
IC.P1 = 2;   IC.Pm1 = IC.P1;
IC.P2 = 3;   IC.Pm2 = IC.P2;
IC.P3 = 1.5; IC.Pm3 = IC.P3;

% Swing Bus Angle
IC.T4 = 0;

% X Data :)
X12  = 0.02;
R12  = 0.002;
X23  = 0.03;
R23  = 0.003;
X34  = 0.015;
R34  = 0.0015;

% Solve Power Flow: 3 eqs (P_inj) and 3 unknowns (T1, T2, T3)
syms    T1 T2 T3
vars = [T1 T2 T3];
Vc1  = IC.V1*exp(j*T1);
Vc2  = IC.V2*exp(j*T2);
Vc3  = IC.V3*exp(j*T3);
Vc4  = IC.V4*exp(j*IC.T4);
eqs  = [IC.P1 == real(Vc1*conj((Vc1-Vc2)/(R12+j*X12)));
        IC.P2 == real(Vc2*conj(((Vc2-Vc1)/(R12+j*X12) + (Vc2-Vc3)/(R23+j*X23))));
        IC.P3 == real(Vc3*conj(((Vc3-Vc2)/(R23+j*X23) + (Vc3-Vc4)/(R34+j*X34))))];
Soln = vpasolve(eqs,vars,[0 0 0]);

% Parse Solution
IC.T1 = double(Soln.T1);
IC.T2 = double(Soln.T2);
IC.T3 = double(Soln.T3);

% Determine Reactive Power
Vc1   = IC.V1*exp(j*IC.T1);
Vc2   = IC.V2*exp(j*IC.T2);
Vc3   = IC.V3*exp(j*IC.T3);
Vc4   = IC.V4*exp(j*IC.T4);
IC.Q1 = imag(Vc1*conj((Vc1-Vc2)/(R12+j*X12)));
IC.Q2 = imag(Vc2*conj(((Vc2-Vc1)/(R12+j*X12) + (Vc2-Vc3)/(R23+j*X23))));
IC.Q3 = imag(Vc3*conj(((Vc3-Vc2)/(R23+j*X23) + (Vc3-Vc4)/(R34+j*X34))));

%% Assign Dynamic Parameters and Solve for IC's

% Gen 1
M1    = 1;
D1    = 0.1;
Xd1   = 0.6;
Xdp1  = 0.1;
Xq1   = 0.4;
Td0p1 = 0.1;
KA1   = 15;
TA1   = 0.05;

% Gen 2
M2    = 0.5;
D2    = 0.05;
Xd2   = 0.5;
Xdp2  = 0.1;
Xq2   = 0.5;
Td0p2 = 0.1;
KA2   = 10;
TA2   = 0.1;

% Gen 3
M3    = .75;
D3    = .05;
Xd3   = 0.4;
Xdp3  = 0.1;
Xq3   = 0.3;
Td0p3 = 0.18;
KA3   = 20;
TA3   = 0.1;

% Solve ICs
syms    d1 d2 d3 eqp1 eqp2 eqp3 Ef1 Ef2 Ef3 Vref1 Vref2 Vref3
vars = [d1 d2 d3 eqp1 eqp2 eqp3 Ef1 Ef2 Ef3 Vref1 Vref2 Vref3];

% Define Alg constraints: V and I
ed1 = IC.V1*sin(d1 - IC.T1);
eq1 = IC.V1*cos(d1 - IC.T1);
ed2 = IC.V2*sin(d2 - IC.T2);
eq2 = IC.V2*cos(d2 - IC.T2);
ed3 = IC.V3*sin(d3 - IC.T3);
eq3 = IC.V3*cos(d3 - IC.T3);
id1 = (eqp1 - eq1)/Xdp1;
iq1 = ed1/Xq1;
id2 = (eqp2 - eq2)/Xdp2;
iq2 = ed2/Xq2;
id3 = (eqp3 - eq3)/Xdp3;
iq3 = ed3/Xq3;

eqs = [... Gen 1
       0     == (Ef1 - (Xd1-Xdp1)*id1-eqp1)/Td0p1;
       0     == (KA1*(Vref1-IC.V1)-Ef1)/TA1;
       IC.P1 == ed1*id1 + eq1*iq1;
       IC.Q1 == id1*eq1 - iq1*ed1;
       ... Gen 2
       0     == (Ef2 - (Xd2-Xdp2)*id2-eqp2)/Td0p2;
       0     == (KA2*(Vref2-IC.V2)-Ef2)/TA2;
       IC.P2 == ed2*id2 + eq2*iq2;
       IC.Q2 == id2*eq2 - iq2*ed2;
       ... Gen 3
       0     == (Ef3 - (Xd3-Xdp3)*id3-eqp3)/Td0p3;
       0     == (KA3*(Vref3-IC.V3)-Ef3)/TA3;
       IC.P3 == ed3*id3 + eq3*iq3;
       IC.Q3 == id3*eq3 - iq3*ed3];

% Solve numerically: Make sure voltages and references are positive
Soln = vpasolve(eqs,vars,[NaN NaN; NaN NaN; NaN NaN; 0 4; 0 4; 0 4; 0 4; 0 4; 0 4; 0 4; 0 4; 0 4]);

% Parse the Solution
IC.d1    = double(Soln.d1);
IC.d2    = double(Soln.d2);
IC.d3    = double(Soln.d3);
IC.eqp1  = double(Soln.eqp1);
IC.eqp2  = double(Soln.eqp2);
IC.eqp3  = double(Soln.eqp3);
IC.Ef1   = double(Soln.Ef1);
IC.Ef2   = double(Soln.Ef2);
IC.Ef3   = double(Soln.Ef3);
IC.Vref1 = double(Soln.Vref1);
IC.Vref2 = double(Soln.Vref2);
IC.Vref3 = double(Soln.Vref3);

%% Set up Dynamic Equations to Simulate
syms w1(t) d1(t) eqp1(t) Ef1(t) V1(t) T1(t)...
     w2(t) d2(t) eqp2(t) Ef2(t) V2(t) T2(t)...
     w3(t) d3(t) eqp3(t) Ef3(t) V3(t) T3(t)...
     V4(t) T4(t)                           ...
     Pm2(t)

% DAE Variables
DAEvars = [w1(t) d1(t) w2(t) d2(t) w3(t) d3(t) eqp1(t) Ef1(t) eqp2(t) Ef2(t)...
           eqp3(t) Ef3(t) V1(t) T1(t) V2(t) T2(t) V3(t) T3(t)];

% Voltages and Currents
ed1 = V1(t)*sin(d1(t) - T1(t));
eq1 = V1(t)*cos(d1(t) - T1(t));
id1 = (eqp1(t) - eq1)/Xdp1;
iq1 = ed1/Xq1;

ed2 = V2(t)*sin(d2(t) - T2(t));
eq2 = V2(t)*cos(d2(t) - T2(t));
id2 = (eqp2(t) - eq2)/Xdp2;
iq2 = ed2/Xq2;

ed3 = V3(t)*sin(d3(t) - T3(t));
eq3 = V3(t)*cos(d3(t) - T3(t));
id3 = (eqp3(t) - eq3)/Xdp3;
iq3 = ed3/Xq3;

% Electrical Powers
Pe1 = ed1*id1+eq1*iq1;
Pe2 = ed2*id2+eq2*iq2;
Pe3 = ed3*id3+eq3*iq3;

% Define Complex Voltages
Vc1 = V1(t)*exp(j*T1(t));
Vc2 = V2(t)*exp(j*T2(t));
Vc3 = V3(t)*exp(j*T3(t));
Vc4 = V4(t)*exp(j*T4(t));

% DAEs
DAEs = [... Dynamic Equations
        diff(w1(t))    == (IC.Pm1 - Pe1 - D1*w1(t))/M1, ...
        diff(d1(t))    == w1(t), ...
        diff(w2(t))    == (IC.Pm2 + Pm2(t) - Pe2 - D2*w2(t))/M2, ...
        diff(d2(t))    == w2(t), ...
        diff(w3(t))    == (IC.Pm3 - Pe3 - D3*w3(t))/M3, ...
        diff(d3(t))    == w3(t), ...
        diff(eqp1(t))  == (Ef1-(Xd1-Xdp1)*id1-eqp1)/Td0p1, ...
        diff(Ef1(t))   == (KA1*(IC.Vref1-V1(t))-Ef1)/TA1, ...
        diff(eqp2(t))  == (Ef2-(Xd2-Xdp2)*id2-eqp2)/Td0p2, ...
        diff(Ef2(t))   == (KA2*(IC.Vref2-V2(t))-Ef2)/TA2, ...
        diff(eqp3(t))  == (Ef3-(Xd3-Xdp3)*id3-eqp3)/Td0p3, ...
        diff(Ef3(t))   == (KA3*(IC.Vref3-V3(t))-Ef3)/TA3, ...
        ... Power Injections
        ed1*id1+eq1*iq1 == real(Vc1*conj( (Vc1-Vc2)/(R12+j*X12))), ...
        eq1*id1-ed1*iq1 == imag(Vc1*conj( (Vc1-Vc2)/(R12+j*X12))), ...
        ed2*id2+eq2*iq2 == real(Vc2*conj( (Vc2-Vc1)/(R12+j*X12) + (Vc2-Vc3)/(R23+j*X23))), ...
        eq2*id2-ed2*iq2 == imag(Vc2*conj( (Vc2-Vc1)/(R12+j*X12) + (Vc2-Vc3)/(R23+j*X23))), ...
        ed3*id3+eq3*iq3 == real(Vc3*conj( (Vc3-Vc2)/(R23+j*X23) + (Vc3-Vc4)/(R34+j*X34))), ...
        eq3*id3-ed3*iq3 == imag(Vc3*conj( (Vc3-Vc2)/(R23+j*X23) + (Vc3-Vc4)/(R34+j*X34)))];

% Set Up System
fF  = daeFunction(DAEs, DAEvars,V4(t),T4(t),Pm2(t));

% Build a vector of Random Noise: Independent Mag and Phase Noise
VN_vec = [0; Rnd_Amp*randn(length(0:dt:tf)-1,1)];
TN_vec = [0; Rnd_Amp*randn(length(0:dt:tf)-1,1)];

% Use a spline add noise
t   = linspace(0,tf,length(VN_vec));
pV4 = spline(t,VN_vec);
pT4 = spline(t,TN_vec);

% Voltage Magnitude and Phase Ocillation
V4  = @(t) (IC.V4 + ppval(pV4,t));
T4  = @(t) (IC.T4 + ppval(pT4,t));
Pm2 = @(t) (IC.Pm2*Osc_Amp*sin(2*pi*Osc_Frq*t));

% Redefine System of Equations
FF   = @(t, Y, YP) fF(t, Y, YP, V4(t), T4(t), Pm2(t));

% Assign the initial conditions
yp0est = zeros(18,1);
y0est  = [0; 
          IC.d1;
          0;
          IC.d2;
          0;
          IC.d3;
          IC.eqp1;
          IC.Ef1;
          IC.eqp2;
          IC.Ef2;
          IC.eqp3;
          IC.Ef3;
          IC.V1;
          IC.T1;
          IC.V2;
          IC.T2;
          IC.V3;
          IC.T3];

% Determine True Initial Conditions
[y0, yp0] = decic(FF, 0, y0est, [], yp0est, []);

% Test Results - Are my initial conditions correct?
error_y  = y0 - y0est;
error_yp = yp0 - yp0est;
error_in = [error_y; error_yp];
if max(error_in) > 10e-8
    disp('*** Initial Conditions do no match ***')
    return
end

%% Simulate the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan  = 0:dt:tf;
[t_out,y_out]  = ode15i(FF,tspan,y0,yp0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! Parse Output !!! %
r.t    = t_out;
r.w1   = y_out(:,1);
r.d1   = y_out(:,2);
r.w2   = y_out(:,3);
r.d2   = y_out(:,4);
r.w3   = y_out(:,5);
r.d3   = y_out(:,6);
r.eqp1 = y_out(:,7);
r.Ef1  = y_out(:,8);
r.eqp2 = y_out(:,9);
r.Ef2  = y_out(:,10);
r.eqp3 = y_out(:,11);
r.Ef3  = y_out(:,12);
r.V1   = y_out(:,13);
r.T1   = y_out(:,14);
r.V2   = y_out(:,15);
r.T2   = y_out(:,16);
r.V3   = y_out(:,17);
r.T3   = y_out(:,18);

% Voltage
r.V4  = IC.V4 + ppval(pV4,r.t);
r.T4  = IC.T4 + ppval(pT4,r.t);

% Terminal Voltage and Current
r.ed1 = r.V1.*sin(r.d1 - r.T1);
r.eq1 = r.V1.*cos(r.d1 - r.T1);
r.id1 = (r.eqp1 - r.eq1)/Xdp1;
r.iq1 = r.ed1/Xq1;
r.ed2 = r.V2.*sin(r.d2 - r.T2);
r.eq2 = r.V2.*cos(r.d2 - r.T2);
r.id2 = (r.eqp2 - r.eq2)/Xdp2;
r.iq2 = r.ed2/Xq2;
r.ed3 = r.V3.*sin(r.d3 - r.T3);
r.eq3 = r.V3.*cos(r.d3 - r.T3);
r.id3 = (r.eqp3 - r.eq3)/Xdp3;
r.iq3 = r.ed3/Xq3;

% Determine Negative Current Injections: Use Transformation
r.Ir1 = -sin(r.d1).*r.id1 - cos(r.d1).*r.iq1;
r.Ii1 =  cos(r.d1).*r.id1 - sin(r.d1).*r.iq1;
r.Ir2 = -sin(r.d2).*r.id2 - cos(r.d2).*r.iq2;
r.Ii2 =  cos(r.d2).*r.id2 - sin(r.d2).*r.iq2;
r.Ir3 = -sin(r.d3).*r.id3 - cos(r.d3).*r.iq3;
r.Ii3 =  cos(r.d3).*r.id3 - sin(r.d3).*r.iq3;

% Convert to Polar
r.Im1 = sqrt((r.Ir1).^2 + (r.Ii1).^2);
r.Ia1 = unwrap(angle(r.Ir1 + j*r.Ii1));
r.Im2 = sqrt((r.Ir2).^2 + (r.Ii2).^2);
r.Ia2 = unwrap(angle(r.Ir2 + j*r.Ii2));
r.Im3 = sqrt((r.Ir3).^2 + (r.Ii3).^2);
r.Ia3 = unwrap(angle(r.Ir3 + j*r.Ii3));

%% Save Simulation Data
  symObj = syms;
  cellfun(@clear,symObj);
  save('Sim_Data')

%% Apply FFT to Voltage and (Negative) Currents

% Define and Apply Hanning Window
L            = length(r.V1);  
Hann_Wind    = hann(L);

% Apply FFT to current and votlage data
N         = length(r.V1);
[~,y_Vm1] = Apply_FFT_N(Hann_Wind.*detrend(r.V1,'constant'),dt,N);
[~,y_Va1] = Apply_FFT_N(Hann_Wind.*detrend(r.T1,'constant'),dt,N);
[~,y_Vm2] = Apply_FFT_N(Hann_Wind.*detrend(r.V2,'constant'),dt,N);
[~,y_Va2] = Apply_FFT_N(Hann_Wind.*detrend(r.T2,'constant'),dt,N);
[~,y_Vm3] = Apply_FFT_N(Hann_Wind.*detrend(r.V3,'constant'),dt,N);
[~,y_Va3] = Apply_FFT_N(Hann_Wind.*detrend(r.T3,'constant'),dt,N);
[~,y_Im1] = Apply_FFT_N(Hann_Wind.*detrend(r.Im1,'constant'),dt,N);
[~,y_Ia1] = Apply_FFT_N(Hann_Wind.*detrend(r.Ia1,'constant'),dt,N);
[~,y_Im2] = Apply_FFT_N(Hann_Wind.*detrend(r.Im2,'constant'),dt,N);
[~,y_Ia2] = Apply_FFT_N(Hann_Wind.*detrend(r.Ia2,'constant'),dt,N);
[~,y_Im3] = Apply_FFT_N(Hann_Wind.*detrend(r.Im3,'constant'),dt,N);
[f,y_Ia3] = Apply_FFT_N(Hann_Wind.*detrend(r.Ia3,'constant'),dt,N);
f_vec     = f;

% Only save 40% of the data
f_ind = round(length(f_vec)/3.3333);
rng_i = 1:f_ind;

% Reset
y_Vm1 = y_Vm1(rng_i);
y_Va1 = y_Va1(rng_i);
y_Vm2 = y_Vm2(rng_i);
y_Va2 = y_Va2(rng_i);
y_Vm3 = y_Vm3(rng_i);
y_Va3 = y_Va3(rng_i);
y_Im1 = y_Im1(rng_i);
y_Ia1 = y_Ia1(rng_i);
y_Im2 = y_Im2(rng_i);
y_Ia2 = y_Ia2(rng_i);
y_Im3 = y_Im3(rng_i);
y_Ia3 = y_Ia3(rng_i);
f_vec = f_vec(rng_i);

%% Predict Currents - Compute Admittance Matrices
syms d_a w_a eqp_a Ef_a V_a T_a Omega_a Pm_a Xq_a Xd_a Xdp_a KA_a TA_a D_a M_a Td0p_a Vref_a
X   = [d_a w_a eqp_a Ef_a];
Uv  = [V_a T_a];
ed = V_a*sin(d_a-T_a);
eq = V_a*cos(d_a-T_a);
id = (eqp_a-eq)/Xdp_a;
iq = ed/Xq_a;
Pe = ed*id + eq*iq;

% System Model
F_vec = [w_a;
    (Pm_a - Pe - D_a*w_a)/M_a;
    (Ef_a - (Xd_a-Xdp_a)*id - eqp_a)/Td0p_a;
    (KA_a*(Vref_a-V_a) - Ef_a)/TA_a];

% Output
Ir    = -sin(d_a)*id - cos(d_a)*iq;
Ii    =  cos(d_a)*id - sin(d_a)*iq;
G_vec = [sqrt(id^2+iq^2); angle(Ir+j*Ii)];
    
% Build Jacobians
JacFX  = jacobian(F_vec,X);
JacFUv = jacobian(F_vec,Uv);
JacGX  = jacobian(G_vec,X);
JacGUv = jacobian(G_vec,Uv);

% Transformation Matrix and Admittance Matrix
Y_a = (JacGX/(j*Omega_a*eye(4) - JacFX))*JacFUv + JacGUv;
Y_f = matlabFunction(Y_a);

% Initialize Empty Predictions
y_Imp1 = zeros(length(f_vec),1);
y_Iap1 = zeros(length(f_vec),1);
y_Imp2 = zeros(length(f_vec),1);
y_Iap2 = zeros(length(f_vec),1);
y_Imp3 = zeros(length(f_vec),1);
y_Iap3 = zeros(length(f_vec),1);

%% Loop
for ii = 1:length(f_vec)
    % Choose Frequency
    Omega_d = f_vec(ii)*2*pi;
    
    % Evaluate Y for Generators 1 and 2
    Y_n1 = Y_f(D1,KA1,M1,Omega_d,TA1,IC.T1,Td0p1,IC.V1,Xd1,Xdp1,Xq1,IC.d1,IC.eqp1);
    Y_n2 = Y_f(D2,KA2,M2,Omega_d,TA2,IC.T2,Td0p2,IC.V2,Xd2,Xdp2,Xq2,IC.d2,IC.eqp2);
    Y_n3 = Y_f(D3,KA3,M3,Omega_d,TA3,IC.T3,Td0p3,IC.V3,Xd3,Xdp3,Xq3,IC.d3,IC.eqp3);
    
    % Test Currents
    I_pred1 = Y_n1*[y_Vm1(ii); y_Va1(ii)];
    I_pred2 = Y_n2*[y_Vm2(ii); y_Va2(ii)];
    I_pred3 = Y_n3*[y_Vm3(ii); y_Va3(ii)];
    
    % Save
    y_Imp1(ii) = I_pred1(1);
    y_Iap1(ii) = I_pred1(2);
    y_Imp2(ii) = I_pred2(1);
    y_Iap2(ii) = I_pred2(2);
    y_Imp3(ii) = I_pred3(1);
    y_Iap3(ii) = I_pred3(2);
end

% Plots
clf

% Generator 1
subplot(3,2,1)
semilogy(f_vec,abs(y_Im1))
hold on
semilogy(f_vec,abs(y_Imp1))
ylim([10^-6 10^0])

subplot(3,2,2)
semilogy(f_vec,abs(y_Ia1))
hold on
semilogy(f_vec,abs(y_Iap1))
ylim([10^-6 10^0])

% Generator 2
subplot(3,2,3)
semilogy(f_vec,abs(y_Im2))
hold on
semilogy(f_vec,abs(y_Imp2))
ylim([10^-6 10^0])

subplot(3,2,4)
semilogy(f_vec,abs(y_Ia2))
hold on
semilogy(f_vec,abs(y_Iap2))
ylim([10^-6 10^0])

% Generator 3
subplot(3,2,5)
semilogy(f_vec,abs(y_Im3))
hold on
semilogy(f_vec,abs(y_Imp3))
ylim([10^-6 10^0])

subplot(3,2,6)
semilogy(f_vec,abs(y_Ia3))
hold on
semilogy(f_vec,abs(y_Iap3))
ylim([10^-6 10^0])
