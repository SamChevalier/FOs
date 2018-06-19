function [data_out,STD_ns] = MAP_SNR(data_in,SNR,d_coi)
% MAP_SNR: Apply noise to the data in the structure "data" of the given
%          SNR. Output the appropriate noise standard deviation.
%
% Inputs:
% 1) data_in   A structure with the following entires:
%                  1) data_in.Vm
%                  2) data_in.Va
%                  3) data_in.Im
%                  4) data_in.Ia
% 2) SNR       The desired SNR
% 3) d_coi     Center of intertia to subtract from angular data. If there
%              is no COI, just set it to 0
%
% Outputs:
% 1) data_out  A data structure with the following (corrupted) entries
%                  1) data_in.Vm + noise
%                  2) data_in.Va + noise
%                  3) data_in.Im + noise
%                  4) data_in.Ia + noise   
% 2) STD_ns    A data structure with the standard deviation of the
%              associated PMU signals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(data_in.Vm);

% Magnitude Data
data_out.Vm = data_in.Vm + randn(n,1)*std(data_in.Vm)/(10^(SNR/20));
data_out.Im = data_in.Im + randn(n,1)*std(data_in.Im)/(10^(SNR/20));

% Angle Data: Subtract out COI (might be 0)
data_out.Va = data_in.Va + randn(n,1)*std(data_in.Va - d_coi)/(10^(SNR/20));
data_out.Ia = data_in.Ia + randn(n,1)*std(data_in.Ia - d_coi)/(10^(SNR/20));

% Collect standard deviations
STD_ns.Vm = std(data_in.Vm)/(10^(SNR/20));
STD_ns.Im = std(data_in.Im)/(10^(SNR/20));
STD_ns.Va = std(data_in.Va - d_coi)/(10^(SNR/20));
STD_ns.Ia = std(data_in.Ia - d_coi)/(10^(SNR/20));
            
end

