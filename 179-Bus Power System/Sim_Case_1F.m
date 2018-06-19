%% Use PSAT to Simulate Case 1 (Source = Bus 4, 65)
%  Sam Chevalier
clear variables; clc;
global Source_V Amplitude_V Frequency_V vref0   ...
       Source_P1 Amplitude_P1 Frequency_P1 pm01 ...
       intstep VP_nos VQ_nos tcorr P0 Q0 PQ_std

%%% Careful! All Bus indices are shifted by 1. Therfore, "Bus 4" %%%
%%% (source Bus) is actually Bus #3 when calling its variables.  %%%

%% Initial Data
WECC_179
npq = size(PQ.con,1);           % Number of load buses
P0  = PQ.con(:,4);              % Initial value of loads' active power
Q0  = PQ.con(:,5);              % Initial value of loads' reactive power
    
%% Solve Power Flow
initpsat;
datafile = 'WECC_179';              
runpsat(datafile,'data');           
Settings.freq   = 60;   
Settings.maxvar = 5000;    % Increase # of Variables
runpsat('pf');             % Run the Almighty Power Flow %

% Parse Power Flow Results Data
n        = Bus.n;
Angles   = DAE.y(Bus.a);
Voltages = DAE.y(Bus.v);
EMFs     = Syn.vf0;

%% Store 3rd Order Generator Initialization Data
%  Save 2 items: eqp and Ef
i1        = Syn.e1q;
Syn3.eqp  = DAE.x(i1(1));
Syn3.Ef   = DAE.x(Exc.vf);

%% Choose Excitation and Mechanical Oscillation Source Buses
%  Sources are generator indices

% Exciter Reference
Source_V    = 1;
Amplitude_V = 0.015;
Frequency_V = 0.86;
vref0       = Exc.vref0(Source_V);

% Mechanical Power 2 (Generator Bus 65)
Source_P1    = 15;
Amplitude_P1 = 0.05;
Frequency_P1 = 0.7;
pm01         = Syn.pm0(Source_P1);

%% Initialize Time Domain Simulation
intstep = 0.0333333;
tcorr   = 1;
PQ_std  = 0.01;
tbegin  = 0;
tfinal  = 120;
runpsat('Vref_perturb','pert');     % Perturbation file
Settings.freq   = 60;               % Change System Freq from default to 60
clpsat.readfile = 1;                % Read data file before running power flow
VP_nos  = zeros(npq,1);             % Vector of noise for load P
VQ_nos  = zeros(npq,1);             % Vector of noise for load Q

% SETTINGS FOR TIME DOMAIN SIMULATION
Settings.coi   = 0;                % Do *NOT* use center of inertia for synchronous machines
Settings.t0    = tbegin;           % Initial simulation time
Settings.tf    = tfinal;           % Final simulation time
Settings.pq2z  = 0;                % Do not convert PQ loads to constant impedance (Z) loads after power flow
Settings.fixt  = 1;                % Enable fixed time-step solver
Settings.tstep = intstep;          % Fixed time step value
nL             = Line.n + Ltc.n + Phs.n + Hvdc.n + Lines.n; % Num circuit elements
Varname.idx    = 1:DAE.n + DAE.m + 2*Bus.n + 6*nL + 2;      % Plot Variable Indexes (ask for current to be included)

% Run Time Domain Simulation
runpsat('td');                     

% Define Output Variables
ix_Va  = DAE.n+1:DAE.n+Bus.n;                                              % Index of voltage angles
ix_Vm  = DAE.n+Bus.n+1:DAE.n+2*Bus.n;                                      % Index of voltage magnitudes
ix_Iij = DAE.n + DAE.m + 2*Bus.n + 4*nL+1:DAE.n + DAE.m + 2*Bus.n + 5*nL;  % Index of line currents
ix_Iji = DAE.n + DAE.m + 2*Bus.n + 5*nL+1:DAE.n + DAE.m + 2*Bus.n + 6*nL;  % Index of line currents
ix_P   = DAE.n + DAE.m + 2*Bus.n + 1:DAE.n + DAE.m + 2*Bus.n+nL;           % Index of line active power flows (i -> j)
ix_Q   = DAE.n + DAE.m + 2*Bus.n +2*nL +1:DAE.n + DAE.m + 2*Bus.n+3*nL;    % Index of line reactive powers flows (i -> j)
ix_de  = Syn.delta;                                                        % Index of rotor angles
ix_om  = Syn.omega;                                                        % Index of generator speed
ix_Pi  = DAE.n + DAE.m + 1:DAE.n + DAE.m + Bus.n;                          % Index of Active Power Injections
ix_Pm  = DAE.n + Syn.pm;  

% Call all of the algebraic variables saved during the TD Simulation
state_ind  = 1:DAE.n;
alg_ind    = (DAE.n + 1):(DAE.n + DAE.m);
Alg_Vars   = Varout.vars(:,alg_ind);
State_Vars = Varout.vars(:,state_ind);

% Output variable values
Va        = Varout.vars(:,ix_Va);   % Va: Bus voltage angles
Vm        = Varout.vars(:,ix_Vm);   % Vm: Bus voltage magnitudes
time      = Varout.t;
deltas    = Varout.vars(:,ix_de);
omegas    = Varout.vars(:,ix_om);
Il_ij     = Varout.vars(:,ix_Iij);
Il_ji     = Varout.vars(:,ix_Iji);
Pl        = Varout.vars(:,ix_P);    % Active Power Flows
Ql        = Varout.vars(:,ix_Q);    % Reactive Power Flows
P         = Varout.vars(:,ix_Pi);   % Active Power Injections
Pm        = Varout.vars(:,ix_Pm);   % Mechanical Powers

%% Other Data
V_ref   = Alg_Vars(:,Exc.vref);
eqp_ind = Syn.e1q;
eqp     = State_Vars(:,eqp_ind(1));

%% Save all Voltage and Current Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  save('Sim_Data','Vm','Va','deltas','time','V_ref','eqp','EMFs','intstep','Syn3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% What are the Active Power Injections at Each Generator Terminal?
Gen_inds = Syn.con(:,1);
P_Gens   = zeros(size(P,1),length(Gen_inds));
for ii = 1:length(Gen_inds)
    Gen_bus      = Gen_inds(ii) - 1; % Because Bus #4 is the 3rd bus
    P_Gens(:,ii) = detrend(P(:,Gen_bus),'constant');
end
