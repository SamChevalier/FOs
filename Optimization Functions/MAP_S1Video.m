function [y_pred0,y_predx] = MAP_S1Video(sol_S1,Y,y_data,dt)
% MAP_S1VIDEO: Iteratively plot the optimizer's Stage 1 solution.
%
% Inputs:
% 1) sol_S1       Optimizer's solution matrix
% 2) Y            Admittance structure
% 3) y_data       Frequency domain data
% 4) dt           Time between iterations
%
% Outputs:
% 1) y_pred0      Initial current prediction
% 2) y_predx      Final current prediction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:size(sol_S1,2)
    Omega = y_data.f*2*pi;
    xc    = sol_S1(:,ii);
    Yn11  = Y.YfC{1,1}(Omega,xc);
    Yn12  = Y.YfC{1,2}(Omega,xc);
    Yn21  = Y.YfC{2,1}(Omega,xc);
    Yn22  = Y.YfC{2,2}(Omega,xc);
    
    % Predict
    y_Imp = Yn11.*y_data.Vm + Yn12.*y_data.Va;
    y_Iap = Yn21.*y_data.Vm + Yn22.*y_data.Va;
    
    % Track first prediction
    if ii == 1
        y_pred0.y_Imp = y_Imp;
        y_pred0.y_Iap = y_Iap;
    end
    
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
        pause(1)
    else
        pause(dt)
    end
end

% Track final prediction
y_predx.y_Imp = y_Imp;
y_predx.y_Iap = y_Iap;

end

