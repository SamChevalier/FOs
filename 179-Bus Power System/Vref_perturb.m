function Vref_perturb(t)
global Source_V Amplitude_V Frequency_V ...
       Source_P1 Amplitude_P1 Frequency_P1 pm01 ...
       Exc vref0 intstep VP_nos VQ_nos tcorr ...
       P0 Q0 PQ_std PQ Syn 
   
% 1. Cause the Reference to Oscillate on Bus 4 (Index 3)
    Exc.vref0(Source_V) = vref0*(1 + Amplitude_V*sin(2*pi*Frequency_V*t));
    
% 2. Cause the Torque to Oscillate on Bus n (Index n-1)
    Syn.pm0(Source_P1) = pm01*(1 + Amplitude_P1*sin(Frequency_P1*2*pi*t));

% 4. Add Load Perturbations
    gamma    = 1/tcorr;                  % Inverse noise correlation time
    nos_std  = PQ_std*sqrt(2*gamma);     % Should be 0.01*sqrt(2*gamma); 
    rnd_vec  = randn(length(P0),1);
    
    % Load noise
    VP_nos   =  VP_nos*(1 - gamma * intstep) + nos_std * sqrt(intstep) * rnd_vec;
    VQ_nos   =  VQ_nos*(1 - gamma * intstep) + nos_std * sqrt(intstep) * rnd_vec;
    
    % Add the noise
    PQ.con(:,4) = P0.*(1+VP_nos);
    PQ.con(:,5) = Q0.*(1+VQ_nos);

end