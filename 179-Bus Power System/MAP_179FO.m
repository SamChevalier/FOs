%% Analyze "Sim_Case_1F" via MAP_Opt
j = sqrt(-1);

% Load Relevant Data ;)
load('SF');
load('Sim_Data');
load('Param_Perts');
load('Admittance_3rd');

% The value of lambda (the regularization parameter) is set equal to
% "lambda_S" which is sufficiently large to find sparse current injections.

%% 1. Initialize Solver Constants

% Frequency settings
freq_data.LowerFB    = [0.68 0.84];  % Set lower frequnecy (band edge) for current injections
freq_data.UpperFB    = [0.72 0.88];  % Set upper frequnecy (band edge) for current injections
freq_data.Max        = 3;            % Max frequnecy
freq_data.dt         = 1/30;         % Time Step

% Optimizer Settings: fmincon
Optns.OptTol = 1e-6;
Optns.FunTol = 1e-6;
Optns.StpTol = 1e-6;
Optns.MaxIts = 1;    % One step, typically
Optns.SolMtd = 3;    % 1 => No Gradient, BFGS Hessian
                     % 2 => Analytical Gradient, BFGS Hessian
                     % 3 => Analytical Gradient, Analytical Hessian

% Optimizer Settings: loop
Optns.op_its = 100;    % Max number of optimization loops
Optns.op_SCP = 0.01;   % Stopping Criteria Percentage: if fval changes less 
                       % than this percentage, stop the optimization loop

%% 2. Perturb System Parameters & Add Measurement Noise
param_dist = 1.5;  % Uniform Distribution (-PD/2  0  +PD/2)

% Call System Data: Call and Convert all Machine Data
WECC_179
initpsat;
datafile = 'WECC_179';              
runpsat(datafile,'data');           
Settings.freq = 60;   
runpsat('pf');

% 2nd Order Machine Parameter Vectors
GIv    = Syn.con(2:end,1);
GIvD   = Syn.con(2:end,1)-1;
Mv     = Syn.con(2:end,18);
Dv     = Syn.con(2:end,19);
Xdpv   = Syn.con(2:end,9);
EMFv   = EMFs(2:end);
n2nd   = length(EMFv);
Vmv    = Vm(:,GIvD);
Vav    = Va(:,GIvD);
Deltav = deltas(:,2:end);

% 3rd Order Generator
G3     = Syn.con(1,1);
GD3    = G3 - 1;
Vm3    = Vm(:,GD3);
Va3    = Va(:,GD3);
Delta3 = deltas(:,1);
M3     = Syn.con(1,18);
D3     = Syn.con(1,19); % 0 - not included in the model
xd3    = Syn.con(1,8);
xdp3   = Syn.con(1,9);
xq3    = Syn.con(1,13);
Td0p3  = Syn.con(1,11);
KA3    = Exc.con(1,5);  % Exciter Gain
TA3    = Exc.con(1,6);  % Exciter TC

% Assemble Data Vectors
GP_Mat2 = [Dv'; EMFv'; Mv'; Xdpv'];
GP_vec3 = [KA3 M3 TA3 Td0p3 xd3 xdp3 xq3]';

% Perturbations: load from memory
%         rnd_Mat2 = param_dist*rand(size(GP_Mat2))     - param_dist/2;
%         rnd_Vec3 = param_dist*rand(length(GP_vec3),1) - param_dist/2;
                     
%% 3. Compute all Currents (Flowing INTO Generators)

% a) Second Order Generator Currents
Im_2nd = zeros(size(Vm,1),n2nd);
Ia_2nd = zeros(size(Vm,1),n2nd);

for ii = 1:n2nd
    I            = (Vmv(:,ii).*exp(j*Vav(:,ii)) - EMFv(ii).*exp(j*Deltav(:,ii)))/(j*Xdpv(ii));
    Im_2nd(:,ii) = abs(I);
    Ia_2nd(:,ii) = unwrap(angle(I));
end

% b) 3rd Order Currents
id = (eqp - Vm3.*cos(Delta3-Va3) )/xdp3;
iq =  (Vm3.*sin(Delta3-Va3))/xq3;

% Convert 3rd Order Gen Currents to Polar
Im_3rd = zeros(size(id,1),1);
Ia_3rd = zeros(size(id,1),1);

for ii = 1:length(id)
    T = -[ sin(Delta3(ii))  cos(Delta3(ii));
          -cos(Delta3(ii))  sin(Delta3(ii))];
    I = T*[id(ii); iq(ii)];
    Im_3rd(ii) = abs(I(1) + j*I(2));
    Ia_3rd(ii) = angle(I(1) + j*I(2));
end

% Unwrap 4th Order Generator
Ia_3rd = unwrap(Ia_3rd);

% Rename Voltages (Ensure an Odd Number)
if mod(size(Vmv,1),2) == 0; mv = 1; else mv = 0; end
    
Vm_2nd = Vmv(1:end-mv,:);
Va_2nd = Vav(1:end-mv,:);
Vm_3rd = Vm3(1:end-mv);
Va_3rd = Va3(1:end-mv);

% Same For Currents
Im_2nd = Im_2nd(1:end-mv,:);
Ia_2nd = Ia_2nd(1:end-mv,:);
Im_3rd = Im_3rd(1:end-mv);
Ia_3rd = Ia_3rd(1:end-mv);

% Misc.
deltas = deltas(1:end-1,:);

% Measurement noise is added to acheive a SNR of 45. All angular data is
% subtracted from the COI, as defined by Milano in PSAT
M_gens = Syn.con(:,18);
d_coi  = (deltas*M_gens)/sum(M_gens);

%% 4. Build 2nd Order Admittance Matrix
syms del_a w_a V_a T_a Ef_a D_Ya M_Ya X_Ya Pm_a Omega
X   = [del_a w_a];
Uv  = [V_a T_a];
w0    = 2*pi*60;
Pe    = V_a*Ef_a*sin(del_a-T_a)/X_Ya;
F_vec = [w0*w_a;
        (Pm_a - Pe - D_Ya*w_a)/M_Ya];
G_vec = [(1/X_Ya)*sqrt(V_a^2 + Ef_a^2 - 2*Ef_a*V_a*cos(T_a-del_a)); 
         atan((V_a*sin(T_a)-Ef_a*sin(del_a))/(V_a*cos(T_a)-Ef_a*cos(del_a)))-pi];

% Output
Y2s = MAP_Build_Y(F_vec,G_vec,X,Uv);

%% 5. Build 3rd Order Admittance Matrix
if ~exist('Y3','var')
    syms del_a w_a eqp_a Ef_a V_a T_a Vr_a Omega Pm_a
    syms KA3_Ya M3_Ya TA3_Ya Td0p3_Ya xd3_Ya xdp3_Ya xq3_Ya Vr3_a
    X   = [del_a w_a eqp_a Ef_a];
    Uv  = [V_a T_a];
    ed = V_a*sin(del_a-T_a);
    eq = V_a*cos(del_a-T_a);
    id = (eqp_a-eq)/xdp3_Ya;
    iq = (ed)/xq3_Ya;
    Pe = ed*id + eq*iq;
    w0 = 2*pi*60;
    F_vec = [w0*w_a;
        (Pm_a - Pe)/M3_Ya;
        (Ef_a - (xd3_Ya-xdp3_Ya)*id - eqp_a)/Td0p3_Ya;
        (KA3_Ya*(Vr3_a-V_a)-Ef_a)/TA3_Ya];
    Ir    = -sin(del_a)*id - cos(del_a)*iq;
    Ii    =  cos(del_a)*id - sin(del_a)*iq;
    G_vec = [sqrt(id^2+iq^2); angle(Ir+j*Ii)];
    
    % Output
    Y3s = MAP_Build_Y(F_vec,G_vec,X,Uv);
    
    % Admittance Matrix: Subs in
    del_a  = Delta3(1);
    w_a    = 0;
    eqp_a  = Syn3.eqp;
    Ef_a   = Syn3.Ef;
    V_a    = Vm3(1);
    T_a    = Va3(1);
    Ys_s   = simplify(subs(Y3s));
    
    % Vector of Uncertain Generator Paramerers
    X_params = [KA3_Ya; M3_Ya; TA3_Ya; Td0p3_Ya; xd3_Ya; xdp3_Ya; xq3_Ya];
    
    %*% Now Define Admittance Matrix Function %*%
    YfC      = cell(2,2);
    YfC{1,1} = matlabFunction(Ys_s(1,1),'vars',{'Omega_a', X_params});
    YfC{1,2} = matlabFunction(Ys_s(1,2),'vars',{'Omega_a', X_params});
    YfC{2,1} = matlabFunction(Ys_s(2,1),'vars',{'Omega_a', X_params});
    YfC{2,2} = matlabFunction(Ys_s(2,2),'vars',{'Omega_a', X_params});
    
    % Package it all up into structure Y
    Y3.X_params = X_params;
    Y3.Ys_s     = Ys_s;
    Y3.YfC      = YfC;
end

% Get Function and Gradient Structures - Update if System Changes!
if exist('S1F3','var') == 1
else
    [S1F3,S2F3] = MAP_LGS(X_params,Ys_s);
end

%% 6. MAP: 2nd Order Generators

for ii = 1:n2nd
    % What is the prior parameter spread? Take, as an input, the standard
    % deviation of the parameter as though its mean value has been normalized
    % to mu = 1.
    %
    % sigma = 0.1  => Almost all of the PDF exists between +-25%
    % sigma = 0.25 => Almost all of the PDF exists between +-60%
    n_params              = size(GP_Mat2,1);
    prior_data.prior_mean = GP_Mat2(:,ii).*(1+rnd_Mat2(:,ii));
    prior_data.prior_std  = 0.01*ones(n_params,1);
    
    % Measurement Noise Strength
    SNR = 45;
    
    % Apply Noise
    data_in.Vm = Vm_2nd(:,ii); data_in.Va = Va_2nd(:,ii); data_in.Im = Im_2nd(:,ii); data_in.Ia = Ia_2nd(:,ii);
    [data,STD_ns] = MAP_SNR(data_in,SNR,d_coi);
    
    % Admittance Matrix: Subs in
    V_a   = Vmv(1,ii);
    del_a = Deltav(1,ii);
    T_a   = Vav(1,ii);
    Ys_s  = simplify(subs(Y2s));
    
    % Vector of Uncertain Generator Paramerers
    X_params = [D_Ya; Ef_a; M_Ya; X_Ya];
    
    %*% Now Define Admittance Matrix Function %*%
    YfC      = cell(2,2);
    YfC{1,1} = matlabFunction(Ys_s(1,1),'vars',{'Omega_a', X_params});
    YfC{1,2} = matlabFunction(Ys_s(1,2),'vars',{'Omega_a', X_params});
    YfC{2,1} = matlabFunction(Ys_s(2,1),'vars',{'Omega_a', X_params});
    YfC{2,2} = matlabFunction(Ys_s(2,2),'vars',{'Omega_a', X_params});
    
    % Package it all up into structure Y
    Y.X_params = X_params;
    Y.Ys_s     = Ys_s;
    Y.YfC      = YfC;
    
    % Get Function and Gradient Structures
    [S1F,S2F] = MAP_LGS(X_params,Ys_s);
    
    % Optimize
    lambda_S2              = 5e5;
    [sol_S1,sol_S2,y_data] = MAP_Opt(data,Y,STD_ns,freq_data,prior_data,Optns,S1F,S2F,lambda_S2);
    
    % Initialize
    if ii == 1
        nf        = length(y_data.f);
        sol_S1M   = zeros(size(sol_S1,1),Optns.op_its,n2nd);
        sol_S2M   = zeros(size(sol_S2,1),Optns.op_its,n2nd);
        Vm_M      = zeros(nf,n2nd);
        Va_M      = zeros(nf,n2nd);
        Im_M      = zeros(nf,n2nd);
        Ia_M      = zeros(nf,n2nd);
        it_termV1 = zeros(n2nd,1);
        it_termV2 = zeros(n2nd,1);
        
        % Predictions
        Im_pM0S1 = zeros(nf,n2nd);
        Ia_pM0S1 = zeros(nf,n2nd);
        Im_pMxS1 = zeros(nf,n2nd);
        Ia_pMxS1 = zeros(nf,n2nd);
        Im_pM0S2 = zeros(nf,n2nd);
        Ia_pM0S2 = zeros(nf,n2nd);
        Im_pMxS2 = zeros(nf,n2nd);
        Ia_pMxS2 = zeros(nf,n2nd);
    end
    
    % Put into 3D Matrix
    it_term1                 = size(sol_S1,2);
    it_term2                 = size(sol_S2,2);
    sol_S1M(:,1:it_term1,ii) = sol_S1;
    sol_S2M(:,1:it_term2,ii) = sol_S2;
    Vm_M(:,ii)      = y_data.Vm;
    Va_M(:,ii)      = y_data.Va;
    Im_M(:,ii)      = y_data.Im;
    Ia_M(:,ii)      = y_data.Ia;
    it_termV1(ii)   = it_term1;
    it_termV2(ii)   = it_term2;
    
    % Convergence Video
    dt = 0.01;
    [y_pred01,y_predx1] = MAP_S1Video(sol_S1,Y,y_data,dt);
    [y_pred02,y_predx2] = MAP_S2Video(sol_S2,Y,y_data,dt,n_params);
    
    % Save all the predictions
    Im_pM0S1(:,ii) = y_pred01.y_Imp;
    Ia_pM0S1(:,ii) = y_pred01.y_Iap;
    Im_pMxS1(:,ii) = y_predx1.y_Imp;
    Ia_pMxS1(:,ii) = y_predx1.y_Iap;
    Im_pM0S2(:,ii) = y_pred01.y_Imp;
    Ia_pM0S2(:,ii) = y_pred01.y_Iap;
    Im_pMxS2(:,ii) = y_predx1.y_Imp;
    Ia_pMxS2(:,ii) = y_predx1.y_Iap;
end

%% 7. MAP: 3rd Order Generator
% What is the prior parameter spread? Take, as an input, the standard
% deviation of the parameter as though its mean value has been normalized
% to mu = 1.
%
% sigma = 0.1  => Almost all of the PDF exists between +-25%
% sigma = 0.25 => Almost all of the PDF exists between +-60%
n_params              = length(GP_vec3);
prior_data.prior_mean = GP_vec3.*(1+rnd_Vec3);
prior_data.prior_std  = 0.01*ones(n_params,1);

% Measurement Noise Strength
SNR = 45;

% Apply Noise
data_in.Vm = Vm_3rd; data_in.Va = Va_3rd; data_in.Im = Im_3rd; data_in.Ia = Ia_3rd;
[data,STD_ns] = MAP_SNR(data_in,SNR,d_coi);

% Optimize
[sol3_S1,sol3_S2,y3_data] = MAP_Opt(data,Y3,STD_ns,freq_data,prior_data,Optns,S1F3,S2F3,lambda_S2);

% Convergence Video
dt = 0.01;
[y3_pred01,y3_predx1] = MAP_S1Video(sol3_S1,Y3,y3_data,dt);
[y3_pred02,y3_predx2] = MAP_S2Video(sol3_S2,Y3,y3_data,dt,n_params);

stop = 1;

%% 8. Plot Initial Mismatches: Percent Difference
clf

fr1 = y_data.f_range(1:5);
fr2 = y_data.f_range(6:11);

% Loop Over 2nd Order Gens
for ii = 1:n2nd

    % Predicted: 4 Values
    P_Imr = real(Im_pM0S1(:,ii));
    P_Imi = imag(Im_pM0S1(:,ii));
    P_Iar = real(Ia_pM0S1(:,ii));
    P_Iai = imag(Ia_pM0S1(:,ii));
    
    % Measured: 4 Values
    M_Imr = real(Im_M(:,ii));
    M_Imi = imag(Im_M(:,ii));
    M_Iar = real(Ia_M(:,ii));
    M_Iai = imag(Ia_M(:,ii));
    
    % Norm of the Difference
    nDiff = sqrt( (P_Imr-M_Imr).^2 + (P_Imi-M_Imi).^2 + (P_Iar-M_Iar).^2 + (P_Iai-M_Iai).^2);
    
    % Norm of their Average
    n1    = sqrt(P_Imr.^2 + P_Imi.^2 + P_Iar.^2 + P_Iai.^2);
    n2    = sqrt(M_Imr.^2 + M_Imi.^2 + M_Iar.^2 + M_Iai.^2);
    nAvg  = 0.5*(n1 + n2);
    
    % Plot
    if ii == 14
        subplot(2,1,1)
        pt = nDiff(fr1(3))./nAvg(fr1(3));
        p1 = semilogy(ii+1,pt,'red*','MarkerSize',10,'LineWidth',1.1);
        line([ii+1 ii+1],[pt 1e-3],'color','black')
        hold on
        
        subplot(2,1,2)
        pt = nDiff(fr2(3))./nAvg(fr2(3));
        semilogy(ii+1,pt,'blackx','MarkerSize',7,'LineWidth',1.1)
        line([ii+1 ii+1],[pt 1e-3],'color','black')
        hold on
    else
        subplot(2,1,1)
        pt = nDiff(fr1(3))./nAvg(fr1(3));
        p2 = semilogy(ii+1,pt,'x','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',7,'LineWidth',1.1);
        line([ii+1 ii+1],[pt 1e-3],'color','black')
        hold on
        
        subplot(2,1,2)
        pt = nDiff(fr2(3))./nAvg(fr2(3));
        p4 = semilogy(ii+1,pt,'x','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',7,'LineWidth',1.1);
        line([ii+1 ii+1],[pt 1e-3],'color','black')
        hold on
    end
end

% 3rd Order Generator
fr = y_data.f_range;

% Predicted: 4 Values
P_Imr = real(y3_pred01.y_Imp);
P_Imi = imag(y3_pred01.y_Imp);
P_Iar = real(y3_pred01.y_Iap);
P_Iai = imag(y3_pred01.y_Iap);

% Measured: 4 Values
M_Imr = real(y3_data.Im);
M_Imi = imag(y3_data.Im);
M_Iar = real(y3_data.Ia);
M_Iai = imag(y3_data.Ia);

% Norm of the Difference
nDiff = sqrt( (P_Imr-M_Imr).^2 + (P_Imi-M_Imi).^2 + (P_Iar-M_Iar).^2 + (P_Iai-M_Iai).^2);

% Norm of their Average
nAvg  = 0.5*sqrt( (P_Imr+M_Imr).^2 + (P_Imi+M_Imi).^2 + (P_Iar+M_Iar).^2 + (P_Iai+M_Iai).^2);

% Plot
subplot(2,1,1)
pt = nDiff(fr1(3))./nAvg(fr1(3));
semilogy(1,pt,'x','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',7,'LineWidth',1.1)
line([1 1],[pt 1e-3],'color','black')
set(gca,'FontName','Times','FontSize',15)
ylabel({'${\rm Prediction \; Error}$','${\rm for} \; f_d = 0.70 \; {\rm Hz}$'},'Interpreter','latex','FontSize',15)
text(0.2,45,'$({\bf a})$','Interpreter','latex','FontSize',15)
legend([p1 p2],{'${\rm Source \; Gen}$','${\rm Non-Source \; Gens}$'},'Interpreter','latex','box','off','Location','NE','FontSize',14)
ylim([1e-1 1e2])

subplot(2,1,2)
pt = nDiff(fr2(3))./nAvg(fr2(3));
p3 = semilogy(1,pt,'red*','MarkerSize',10,'LineWidth',1.1);
line([1 1],[pt 1e-3],'color','black')
set(gca,'FontName','Times','FontSize',15)
xlabel('$\rm{Generator \; Index}$','Interpreter','latex','FontSize',15)
ylabel({'${\rm Prediction \; Error}$','${\rm for} \; f_d = 0.86 \; {\rm Hz}$'},'Interpreter','latex','FontSize',15)
ylim([1e-1 1e2])
legend([p3 p4],{'${\rm Source \; Gen}$','${\rm Non-Source \; Gens}$'},'Interpreter','latex','box','off','Location','NE','FontSize',14)
text(0.2,45,'$({\bf b})$','Interpreter','latex','FontSize',15)
set(gcf,'Units','inches','Position',[0 0 9 4])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% 9. Plot Final Current Injections
clf

fr1 = y_data.f_range(1:5);
fr2 = y_data.f_range(6:11);

% Loop Over 2nd Order Gens
for ii = 1:n2nd
    Inj = sol_S2M(5:end,it_termV2(ii),ii);
    Inj = sqrt(Inj(1:11).^2 + Inj(12:22).^2 + Inj(23:33).^2 + Inj(34:44).^2);
    Inj_70 = Inj(3);
    Inj_86 = Inj(8);
    
    % Plot
    if ii == 14
        subplot(2,1,1)
        pt = Inj_70;
        p1 = semilogy(ii+1,pt,'red*','MarkerSize',10,'LineWidth',1.1);
        line([ii+1 ii+1],[pt 1e-6],'color','black')
        hold on
        
        subplot(2,1,2)
        pt = Inj_86;
        semilogy(ii+1,pt,'blackx','MarkerSize',7,'LineWidth',1.1)
        line([ii+1 ii+1],[pt 1e-6],'color','black')
        hold on
    else
        subplot(2,1,1)
        pt = Inj_70;
        p2 = semilogy(ii+1,pt,'x','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',7,'LineWidth',1.1);
        line([ii+1 ii+1],[pt 1e-6],'color','black')
        hold on
        
        subplot(2,1,2)
        pt = Inj_86;
        p4 = semilogy(ii+1,pt,'x','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',7,'LineWidth',1.1);
        line([ii+1 ii+1],[pt 1e-6],'color','black')
        hold on
    end
end

% 3rd Order Generator
Inj = sol3_S2(8:end,end);
Inj = sqrt(Inj(1:11).^2 + Inj(12:22).^2 + Inj(23:33).^2 + Inj(34:44).^2);
Inj_70 = Inj(3);
Inj_86 = Inj(8);

subplot(2,1,1)
pt = Inj_70;
p2 = semilogy(1,pt,'blackx','MarkerSize',7,'LineWidth',1.1);
line([1 1],[pt 1e-6],'color','black')
ylim([4e-6 1])
set(gca,'FontName','Times','FontSize',15)
ylabel({'${\rm Injection \;}\left\Vert \mathcal{I}\right\Vert$','${\rm at} \; f_d=0.70 \; {\rm Hz}$'},'Interpreter','latex','FontSize',15)
text(-0.75,0.3,'$({\bf a})$','Interpreter','latex','FontSize',15)
legend([p1 p2],{'${\rm Source \; Gen}$','${\rm Non-Source \; Gens}$'},'Interpreter','latex','box','off','Location','NE','FontSize',14)
ylim([1e-4 1])
xlim([-1 30])

subplot(2,1,2)
pt = Inj_86;
p3 = semilogy(1,pt,'red*','MarkerSize',10,'LineWidth',1.1);
line([1 1],[pt 1e-6],'color','black')
ylim([4e-6 1])
set(gca,'FontName','Times','FontSize',15)
xlabel('$\rm{Generator \; Index}$','Interpreter','latex','FontSize',15)
ylabel({'${\rm Injection \;}\left\Vert \mathcal{I}\right\Vert$','${\rm at} \; f_d=0.86 \; {\rm Hz}$'},'Interpreter','latex','FontSize',15)
set(gcf,'Units','inches','Position',[0 0 9 4])
ylim([1e-4 1])
xlim([-1 30])
legend([p3 p4],{'${\rm Source \; Gen}$','${\rm Non-Source \; Gens}$'},'Interpreter','latex','box','off','Location','NE','FontSize',14)
text(-0.75,0.3,'$({\bf b})$','Interpreter','latex','FontSize',15)
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% 10. Plot Initial & Final Mismatches at 2nd Order Generator #2: Magnitude
clf
num_g2 = 2;
c1     = [0 0.4470 0.7410];

% Current Magnitude: Start of S1
subplot(2,1,1)
semilogy(y_data.f,abs(Im_pM0S1(:,num_g2)).^2,'Linewidth',1.5,'color',[c1 0.5]);
hold on
semilogy(y_data.f,abs(Im_M(:,num_g2)).^2,'Linewidth',0.5,'color','black');
set(gca,'FontName','Times','FontSize',15)
rectangle('Position',[0.68 10^-15 0.04 10],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
rectangle('Position',[0.84 10^-15 0.04 10],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',15)
ylabel({'${\tilde {\rm I}} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',15)
text(0.02,0.5e-1,'$({\bf a})$','Interpreter','latex','FontSize',15)
ca = gca;
ca.XTick = [0 0.5 0.7 0.86 1.5000 2 2.5000 3];
ylim([1e-10 1])
xlim([0 2])

% Current Magnitude: End of S1
subplot(2,1,2)
semilogy(y_data.f,abs(Im_pMxS1(:,num_g2)).^2,'Linewidth',1.5,'color',[c1 0.5]);
hold on
semilogy(y_data.f,abs(Im_M(:,num_g2)).^2,'Linewidth',0.5,'color','black');
set(gca,'FontName','Times','FontSize',15)
rectangle('Position',[0.68 10^-15 0.04 10],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
rectangle('Position',[0.84 10^-15 0.04 10],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',15)
ylabel({'${\tilde {\rm I}} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',15)
xlabel('${\rm Frequency \; (Hz)}$','Interpreter','latex','FontSize',15)
text(0.02,0.5e-1,'$({\bf b})$','Interpreter','latex','FontSize',15)
ca = gca;
ca.XTick = [0 0.5 0.7 0.86 1.5000 2 2.5000 3];
set(gcf,'Units','inches','Position',[0 0 9 5])
ylim([1e-10 1])
xlim([0 2])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% 11. Plot Initial & Final Mismatches at 2nd Order Generator #14: Magnitude
clf
num_g2 = 14;
c1     = [0 0.4470 0.7410];

% Current Magnitude: Start of S1
subplot(2,1,1)
semilogy(y_data.f,abs(Im_pM0S1(:,num_g2)).^2,'Linewidth',1.5,'color',[c1 0.5]);
hold on
semilogy(y_data.f,abs(Im_M(:,num_g2)).^2,'Linewidth',0.5,'color','black');
set(gca,'FontName','Times','FontSize',15)
ylabel({'${\tilde {\rm I}} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',15)
text(0.03,1e0,'$({\bf a})$','Interpreter','latex','FontSize',15)
rectangle('Position',[0.68 10^-15 0.04 10],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
rectangle('Position',[0.84 10^-15 0.04 10],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',15)
ca = gca;
ca.XTick = [0 0.5 0.7 0.86 1.5000 2 2.5000 3];
ca.YTick = [1e-5 1e0];
ca.YMinorTick = 'off';
ylim([1e-7 1e1])
xlim([0 2])

% Current Magnitude: End of S1
subplot(2,1,2)
semilogy(y_data.f,abs(Im_pMxS1(:,num_g2)).^2,'Linewidth',1.5,'color',[c1 0.5]);
hold on
semilogy(y_data.f,abs(Im_M(:,num_g2)).^2,'Linewidth',0.5,'color','black');
set(gca,'FontName','Times','FontSize',15)
ylabel({'${\tilde {\rm I}} \;\; {\rm PSD} $'},'Interpreter','latex','FontSize',15)
xlabel('${\rm Frequency \; (Hz)}$','Interpreter','latex','FontSize',15)
text(0.03,1e0,'$({\bf b})$','Interpreter','latex','FontSize',15)
rectangle('Position',[0.68 10^-15 0.04 10],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
rectangle('Position',[0.84 10^-15 0.04 10],'FaceColor',[0.8500    0.3250    0.0980 .3],'EdgeColor','none','LineWidth',0.1)
legend({'${\rm Predicted}$','${\rm Measured}$'},'Interpreter','latex','box','off','Location','NE','FontSize',15)
ca = gca;
ca.XTick = [0 0.5 0.7 0.86 1.5000 2 2.5000 3];
ca.YTick = [1e-5 1e0];
ca.YMinorTick = 'off';
set(gcf,'Units','inches','Position',[0 0 9 5])
ylim([1e-7 1e1])
xlim([0 2])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"
