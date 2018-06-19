%% Analyze "Sim_3_Gen" via MAP_Opt
j = sqrt(-1);

% Load Relevant Data ;)
load('Sim_Data');
load('Admittance_Data');
load('Prior_Data');

%% 1. Initialize Solver Constants

% Frequency settings
freq_data.LowerFB    = 0.48;  % Set lower frequnecy (band edge) for current injections
freq_data.UpperFB    = 0.52;  % Set upper frequnecy (band edge) for current injections
freq_data.Max        = 3;     % Max frequnecy
freq_data.dt         = dt;    % Time Step

% Optimizer Settings: fmincon
Optns.OptTol = 1e-6;
Optns.FunTol = 1e-6;
Optns.StpTol = 1e-6;
Optns.MaxIts = 1;    % One step, typically
Optns.SolMtd = 3;    % 1 => No Gradient, BFGS Hessian
                     % 2 => Analytical Gradient, BFGS Hessian
                     % 3 => Analytical Gradient, Analytical Hessian

% Optimizer Settings: loop
Optns.op_its = 100;   % Max number of optimization loops
Optns.op_SCP = 0.0025; % Stopping Criteria Percentage: if fval changes less 
                      % than this percentage, stop the optimization loop

%% 2. Perturb System Parameters & Add Measurement Noise
param_dist = 1.5; % Uniform Distribution (-PD/2  0  +PD/2)
GP_vec1    = [D1 KA1 M1 TA1 Td0p1 Xd1 Xdp1 Xq1]';
GP_vec2    = [D2 KA2 M2 TA2 Td0p2 Xd2 Xdp2 Xq2]';
GP_vec3    = [D3 KA3 M3 TA3 Td0p3 Xd3 Xdp3 Xq3]';
rnd_vec1   = param_dist*rand(length(GP_vec1),1) - param_dist/2;
rnd_vec2   = param_dist*rand(length(GP_vec2),1) - param_dist/2;
rnd_vec3   = param_dist*rand(length(GP_vec3),1) - param_dist/2;

% What is the prior parameter spread? Take, as an input, the standard
% deviation of the parameter as though its mean value has been normalized
% to mu = 1.
%
% sigma = 0.1  => Almost all of the PDF exists between +-25%
% sigma = 0.25 => Almost all of the PDF exists between +-60%
%
%
% ~~ The following code is loaded from memory for plotting convenience ~~
%         prior_data1.prior_mean = GP_vec1.*(1+rnd_vec1);
%         prior_data1.prior_std  = ones(8,1);
% 
%         prior_data2.prior_mean = GP_vec2.*(1+rnd_vec2);
%         prior_data2.prior_std  = ones(8,1);
% 
%         prior_data3.prior_mean = GP_vec3.*(1+rnd_vec3);
%         prior_data3.prior_std  = ones(8,1);

% Noise Strength
SNR   = 45;
d_coi = 0;

% Apply Noise: G1
data_in.Vm = r.V1; data_in.Va = unwrap(r.T1); data_in.Im = r.Im1; data_in.Ia = unwrap(r.Ia1);
[data1,STD_ns1] = MAP_SNR(data_in,SNR,d_coi);

% Apply Noise: G2
data_in.Vm = r.V2; data_in.Va = unwrap(r.T2); data_in.Im = r.Im2; data_in.Ia = unwrap(r.Ia2);
[data2,STD_ns2] = MAP_SNR(data_in,SNR,d_coi);

% Apply Noise: G3
data_in.Vm = r.V3; data_in.Va = unwrap(r.T3); data_in.Im = r.Im3; data_in.Ia = unwrap(r.Ia3);
[data3,STD_ns3] = MAP_SNR(data_in,SNR,d_coi);

%% 3. Build Admittance Matrix
if ~exist('Y1','var')
    
    % Generator varaibles and DAE model
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
    
    % Output
    Ys = MAP_Build_Y(F_vec,G_vec,X,Uv);
    
    % Generator 1 Substitution
    d_a   = IC.d1;
    eqp_a = IC.eqp1;
    V_a   = IC.V1;
    T_a   = IC.T1;
    Y1s_s = simplify(subs(Ys));
    
    % Generator 2 Substitution
    d_a   = IC.d2;
    eqp_a = IC.eqp2;
    V_a   = IC.V2;
    T_a   = IC.T2;
    Y2s_s = simplify(subs(Ys));
    
    % Generator 2 Substitution
    d_a   = IC.d3;
    eqp_a = IC.eqp3;
    V_a   = IC.V3;
    T_a   = IC.T3;
    Y3s_s = simplify(subs(Ys));
    
    % Vector of Uncertain Generator Paramerers
    X_params = [D_a; KA_a; M_a; TA_a; Td0p_a; Xd_a; Xdp_a; Xq_a];
    n_params = length(X_params);
    
    %*% Now Define Admittance Matrix Function: G1 %*%
    YfC1      = cell(2,2);
    YfC1{1,1} = matlabFunction(Y1s_s(1,1),'vars',{'Omega_a', X_params});
    YfC1{1,2} = matlabFunction(Y1s_s(1,2),'vars',{'Omega_a', X_params});
    YfC1{2,1} = matlabFunction(Y1s_s(2,1),'vars',{'Omega_a', X_params});
    YfC1{2,2} = matlabFunction(Y1s_s(2,2),'vars',{'Omega_a', X_params});
    
    %*% Now Define Admittance Matrix Function: G2 %*%
    YfC2      = cell(2,2);
    YfC2{1,1} = matlabFunction(Y2s_s(1,1),'vars',{'Omega_a', X_params});
    YfC2{1,2} = matlabFunction(Y2s_s(1,2),'vars',{'Omega_a', X_params});
    YfC2{2,1} = matlabFunction(Y2s_s(2,1),'vars',{'Omega_a', X_params});
    YfC2{2,2} = matlabFunction(Y2s_s(2,2),'vars',{'Omega_a', X_params});
    
    %*% Now Define Admittance Matrix Function: G3 %*%
    YfC3      = cell(2,2);
    YfC3{1,1} = matlabFunction(Y3s_s(1,1),'vars',{'Omega_a', X_params});
    YfC3{1,2} = matlabFunction(Y3s_s(1,2),'vars',{'Omega_a', X_params});
    YfC3{2,1} = matlabFunction(Y3s_s(2,1),'vars',{'Omega_a', X_params});
    YfC3{2,2} = matlabFunction(Y3s_s(2,2),'vars',{'Omega_a', X_params});
    
    % Package it all up into structure Y
    Y1.X_params = X_params;
    Y1.Ys_s     = Y1s_s;
    Y1.YfC      = YfC1;
    Y2.X_params = X_params;
    Y2.Ys_s     = Y2s_s;
    Y2.YfC      = YfC2;
    Y3.X_params = X_params;
    Y3.Ys_s     = Y3s_s;
    Y3.YfC      = YfC3;
end

%% 4. Get Function and Gradient Structures - Update if System Changes!
if exist('S1F_1') == 1
    % It is loaded already
elseif exist('SF.mat','file')
    load('SF');
else
    [S1F_1,S2F_1] = MAP_LGS(Y1.X_params,Y1.Ys_s);
    [S1F_2,S2F_2] = MAP_LGS(Y2.X_params,Y2.Ys_s);
    [S1F_3,S2F_3] = MAP_LGS(Y3.X_params,Y3.Ys_s);
end

%% 5. Optimize Generator 1
lamda_S2                  = 5e6; % From XVal results
[sol1_S1,sol1_S2,y1_data] = MAP_Opt(data1,Y1,STD_ns1,freq_data,prior_data1,Optns,S1F_1,S2F_1,lamda_S2);

% Video
n_params                = length(prior_data1.prior_mean);
[y1_pred0S1,y1_predxS1] = MAP_S1Video(sol1_S1,Y1,y1_data,0.2);
[y1_pred0S2,y1_predxS2] = MAP_S2Video(sol1_S2,Y1,y1_data,0.2,n_params);

%% 6. Optimize Generator 2
[sol2_S1,sol2_S2,y2_data] = MAP_Opt(data2,Y2,STD_ns2,freq_data,prior_data2,Optns,S1F_2,S2F_2,lamda_S2);

% Video
n_params                = length(prior_data1.prior_mean);
[y2_pred0S1,y2_predxS1] = MAP_S1Video(sol2_S1,Y2,y2_data,0.1);
[y2_pred0S2,y2_predxS2] = MAP_S2Video(sol2_S2,Y2,y2_data,0.1,n_params);

%% 7. Optimize Generator 3
[sol3_S1,sol3_S2,y3_data] = MAP_Opt(data3,Y3,STD_ns3,freq_data,prior_data3,Optns,S1F_3,S2F_3,lamda_S2);

% Video
n_params                = length(prior_data1.prior_mean);
[y3_pred0S1,y3_predxS1] = MAP_S1Video(sol3_S1,Y3,y3_data,0.2);
[y3_pred0S2,y3_predxS2] = MAP_S2Video(sol3_S2,Y3,y3_data,0.2,n_params);

%% 8a. Quick Plot: Initial Mismatch
clf

% Generator 1
subplot(3,2,1)
semilogy(y1_data.f,abs(y1_data.Im))
hold on
semilogy(y1_data.f,abs(y1_pred0S1.y_Imp))
ylim([10^-6 10^0])

subplot(3,2,2)
semilogy(y1_data.f,abs(y1_data.Ia))
hold on
semilogy(y1_data.f,abs(y1_pred0S1.y_Iap))
ylim([10^-6 10^0])

% Generator 2
subplot(3,2,3)
semilogy(y1_data.f,abs(y2_data.Im))
hold on
semilogy(y1_data.f,abs(y2_pred0S1.y_Imp))
ylim([10^-6 10^0])

subplot(3,2,4)
semilogy(y1_data.f,abs(y2_data.Ia))
hold on
semilogy(y1_data.f,abs(y2_pred0S1.y_Iap))
ylim([10^-6 10^0])

% Generator 3
subplot(3,2,5)
semilogy(y1_data.f,abs(y3_data.Im))
hold on
semilogy(y1_data.f,abs(y3_pred0S1.y_Imp))
ylim([10^-6 10^0])

subplot(3,2,6)
semilogy(y1_data.f,abs(y3_data.Ia))
hold on
semilogy(y1_data.f,abs(y3_pred0S1.y_Iap))
ylim([10^-6 10^0])

%% 8b. Quick Plot: Final S1 Mismatch
clf

% Generator 1
subplot(3,2,1)
semilogy(y1_data.f,abs(y1_data.Im))
hold on
semilogy(y1_data.f,abs(y1_predxS1.y_Imp))
ylim([10^-6 10^0])

subplot(3,2,2)
semilogy(y1_data.f,abs(y1_data.Ia))
hold on
semilogy(y1_data.f,abs(y1_predxS1.y_Iap))
ylim([10^-6 10^0])

% Generator 2
subplot(3,2,3)
semilogy(y1_data.f,abs(y2_data.Im))
hold on
semilogy(y1_data.f,abs(y2_predxS1.y_Imp))
ylim([10^-6 10^0])

subplot(3,2,4)
semilogy(y1_data.f,abs(y2_data.Ia))
hold on
semilogy(y1_data.f,abs(y2_predxS1.y_Iap))
ylim([10^-6 10^0])

% Generator 3
subplot(3,2,5)
semilogy(y1_data.f,abs(y3_data.Im))
hold on
semilogy(y1_data.f,abs(y3_predxS1.y_Imp))
ylim([10^-6 10^0])

subplot(3,2,6)
semilogy(y1_data.f,abs(y3_data.Ia))
hold on
semilogy(y1_data.f,abs(y3_predxS1.y_Iap))
ylim([10^-6 10^0])

%% 9. Paper Plot 1: Initial Mismatch
clf
f_vec = y1_data.f;

% Plot Colors
bl = [0         0.4470    0.7410];
rd = [0.8500    0.3250    0.0980];
c1 = [0 0.4470 0.7410];

%%%%%% ----------- Generator 1 ----------- %%%%%%
subplot(3,2,1)
semilogy(f_vec,abs(y1_pred0S1.y_Imp).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y1_data.Im).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde {\rm I}} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
ylim([10^-9 10^-1])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 1$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-2.1,'$({\bf a})$','Interpreter','latex','FontSize',14)

subplot(3,2,2)
semilogy(f_vec,abs(y1_pred0S1.y_Iap).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y1_data.Ia).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde \phi} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
ylim([10^-9 10^-1.5])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 1$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-2.6,'$({\bf b})$','Interpreter','latex','FontSize',14)

%%%%%% ----------- Generator 2 ----------- %%%%%%
subplot(3,2,3)
semilogy(f_vec,abs(y2_pred0S1.y_Imp).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y2_data.Im).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde {\rm I}} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
ylim([10^-10 10^-1])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 2$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-2.2,'$({\bf c})$','Interpreter','latex','FontSize',14)

subplot(3,2,4)
semilogy(f_vec,abs(y2_pred0S1.y_Iap).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y2_data.Ia).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde \phi} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
ylim([10^-10 10^-2])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 2$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-3,'$({\bf d})$','Interpreter','latex','FontSize',14)

%%%%%% ----------- Generator 3 ----------- %%%%%%
subplot(3,2,5)
semilogy(f_vec,abs(y3_pred0S1.y_Imp).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y3_data.Im).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde {\rm I}} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
xlabel('$\rm{Frequency \; (Hz)}$','Interpreter','latex','FontSize',14)
ylim([10^-7 10^-0.5])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 3$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-1.3,'$({\bf e})$','Interpreter','latex','FontSize',14)

subplot(3,2,6)
semilogy(f_vec,abs(y3_pred0S1.y_Iap).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y3_data.Ia).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde \phi} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
xlabel('$\rm{Frequency \; (Hz)}$','Interpreter','latex','FontSize',14)
ylim([10^-7 10^-2])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 3$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-2.65,'$({\bf f})$','Interpreter','latex','FontSize',14)

% Size
set(gcf,'Units','inches','Position',[0 0 10 5.5])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% 10. Paper Plot 2: End of Stage 1 Mismatch
clf
f_vec = y1_data.f;

% Plot Colors
bl = [0         0.4470    0.7410];
rd = [0.8500    0.3250    0.0980];

%%%%%% ----------- Generator 1 ----------- %%%%%%
subplot(3,2,1)
semilogy(f_vec,abs(y1_predxS1.y_Imp).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y1_data.Im).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde {\rm I}} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
ylim([10^-9 10^-1])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 1$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-2.1,'$({\bf a})$','Interpreter','latex','FontSize',14)

subplot(3,2,2)
semilogy(f_vec,abs(y1_predxS1.y_Iap).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y1_data.Ia).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde \phi} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
ylim([10^-9 10^-1.5])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 1$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-2.6,'$({\bf b})$','Interpreter','latex','FontSize',14)

%%%%%% ----------- Generator 2 ----------- %%%%%%
subplot(3,2,3)
semilogy(f_vec,abs(y2_predxS1.y_Imp).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y2_data.Im).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde {\rm I}} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
ylim([10^-10 10^-1])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 2$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-2.2,'$({\bf c})$','Interpreter','latex','FontSize',14)

subplot(3,2,4)
semilogy(f_vec,abs(y2_predxS1.y_Iap).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y2_data.Ia).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde \phi} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
ylim([10^-10 10^-2])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 2$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-3,'$({\bf d})$','Interpreter','latex','FontSize',14)

%%%%%% ----------- Generator 3 ----------- %%%%%%
subplot(3,2,5)
semilogy(f_vec,abs(y3_predxS1.y_Imp).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y3_data.Im).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde {\rm I}} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
xlabel('$\rm{Frequency \; (Hz)}$','Interpreter','latex','FontSize',14)
ylim([10^-7 10^-0.5])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 3$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-1.3,'$({\bf e})$','Interpreter','latex','FontSize',14)

subplot(3,2,6)
semilogy(f_vec,abs(y3_predxS1.y_Iap).^2,'Linewidth',1.5,'color',[c1 0.5])
hold on
semilogy(f_vec,abs(y3_data.Ia).^2,'Linewidth',0.5,'color','black')
rectangle('Position',[0.48 10^-15 0.04 1],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
ylabel({'${\tilde \phi} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',14)
xlabel('$\rm{Frequency \; (Hz)}$','Interpreter','latex','FontSize',14)
ylim([10^-7 10^-2])
xlim([0.3 0.7])
set(gca,'FontName','Times','FontSize',13)
set(gca,'ytick',[])
title({'$\rm{\bf Generator} \;\; \bf 3$'},'Interpreter','latex','FontSize',14)
text(0.305,10^-2.65,'$({\bf f})$','Interpreter','latex','FontSize',14)

% Size
set(gcf,'Units','inches','Position',[0 0 10 5.5])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% 11. Paper Plot 3: Stage 2 Injections
clf

% Parse Currents
I1  = sol1_S2(9:end,end);
I2  = sol2_S2(9:end,end);
I3  = sol3_S2(9:end,end);
I1n = sqrt(I1(1:5).^2 + I1(6:10).^2 + I1(11:15).^2 + I1(16:20).^2);
I2n = sqrt(I2(1:5).^2 + I2(6:10).^2 + I2(11:15).^2 + I2(16:20).^2);
I3n = sqrt(I3(1:5).^2 + I3(6:10).^2 + I3(11:15).^2 + I3(16:20).^2);

% Frequnecies
freq_vals = y1_data.f(y1_data.f_range);

clf
subplot(1,3,1)
semilogy(freq_vals,I1n,'*','markersize',6,'LineWidth',1.1)
hold on
for ii = 1:length(freq_vals)
    line([freq_vals(ii) freq_vals(ii)],[abs(I1n(ii)) 1e-5],'color','black')
end
set(gca,'FontName','Times','FontSize',14)
ylim([2.8e-5 1e-1])
xlim([0.47 0.53])
ylabel({'${{\rm Injection \;\; Mag. \;}\left\Vert{\mathcal I}\right\Vert}$'},'Interpreter','latex','FontSize',14)
xlabel('$\rm{Frequency \; (Hz)}$','Interpreter','latex','FontSize',14)
title({'$\rm{\bf Generator} \;\; \bf 1$'},'Interpreter','latex','FontSize',14)
set(gca,'YTick',[1.0000e-04 1.0000e-03 0.0100])
set(gca,'YGrid','on');

subplot(1,3,2)
semilogy(freq_vals,I2n,'*','markersize',6,'LineWidth',1.1)
hold on
for ii = 1:length(freq_vals)
    line([freq_vals(ii) freq_vals(ii)],[abs(I2n(ii)) 1e-5],'color','black')
end
set(gca,'FontName','Times','FontSize',14)
ylim([2.8e-5 1e-1])
xlim([0.47 0.53])
set(gca,'YTick',[1.0000e-04 1.0000e-03 0.0100])
xlabel('$\rm{Frequency \; (Hz)}$','Interpreter','latex','FontSize',14)
title({'$\rm{\bf Generator} \;\; \bf 2$'},'Interpreter','latex','FontSize',14)
set(gca,'YGrid','on','YTickLabel',{});

subplot(1,3,3)
semilogy(freq_vals,I3n,'*','markersize',6,'LineWidth',1.1)
for ii = 1:length(freq_vals)
    line([freq_vals(ii) freq_vals(ii)],[abs(I3n(ii)) 1e-5],'color','black')
end
set(gca,'FontName','Times','FontSize',14)
xlabel('$\rm{Frequency \; (Hz)}$','Interpreter','latex','FontSize',14)
ylim([2.8e-5 1e-1])
xlim([0.47 0.53])
set(gca,'YTick',[1.0000e-04 1.0000e-03 0.0100])
title({'$\rm{\bf Generator} \;\; \bf 3$'},'Interpreter','latex','FontSize',14)
set(gca,'YGrid','on','YTickLabel',{});
set(gca,'YMinorGrid','on')
set(gcf,'Units','inches','Position',[0 0 9 2])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

