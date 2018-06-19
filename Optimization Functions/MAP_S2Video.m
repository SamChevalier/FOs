function [y_pred0,y_predx] = MAP_S2Video(sol_S2,Y,y_data,dt,n_params)
% MAP_S1VIDEO: Iteratively plot the optimizer's Stage 2 solution.
%
% Inputs:
% 1) sol_S2       Optimizer's solution matrix
% 2) Y            Admittance structure
% 3) y_data       Frequency domain data
% 4) dt           Time between iterations
% 5) n_params     Number of parameters
%
% Outputs:
% 1) y_pred0      Initial current prediction
% 2) y_predx      Final current prediction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:size(sol_S2,2)
    Omega = y_data.f*2*pi;
    xc    = sol_S2(1:n_params,ii);
    Yn11  = Y.YfC{1,1}(Omega,xc);
    Yn12  = Y.YfC{1,2}(Omega,xc);
    Yn21  = Y.YfC{2,1}(Omega,xc);
    Yn22  = Y.YfC{2,2}(Omega,xc);
    
    % Parse currents
    nf     = length(y_data.f_range);
    Inj_mr = sol_S2(n_params + 0*nf + (1:nf),ii);
    Inj_mi = sol_S2(n_params + 1*nf + (1:nf),ii);
    Inj_pr = sol_S2(n_params + 2*nf + (1:nf),ii);
    Inj_pi = sol_S2(n_params + 3*nf + (1:nf),ii);
    Inj_Im = Inj_mr + j*Inj_mi;
    Inj_Ip = Inj_pr + j*Inj_pi;
    
    % Predict
    y_Imp = Yn11.*y_data.Vm + Yn12.*y_data.Va;
    y_Iap = Yn21.*y_data.Vm + Yn22.*y_data.Va;
    
    % Track first prediction
    if ii == 1
        y_pred0.y_Imp = y_Imp;
        y_pred0.y_Iap = y_Iap;
    end
    
    % Add currents
    y_Imp(y_data.f_range) = y_Imp(y_data.f_range) + Inj_Im;
    y_Iap(y_data.f_range) = y_Iap(y_data.f_range) + Inj_Ip;
    
    % Plot Magnitudes
    clf
    subplot(2,1,1)
    semilogy(y_data.f,abs(y_data.Im))
    hold on
    semilogy(y_data.f,abs(y_Imp))
    
    % Plot Angles
    subplot(2,1,2)
    semilogy(y_data.f,abs(y_data.Ia))
    hold on
    semilogy(y_data.f,abs(y_Iap))
    
    % Pause
    if ii == 1
        pause(0.5)
    else
        pause(dt)
    end
end

% Track final prediction
y_predx.y_Imp = y_Imp;
y_predx.y_Iap = y_Iap;

end

