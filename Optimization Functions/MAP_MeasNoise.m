function [std_w] = MAP_MeasNoise(N,STD_ns)
% MAP_MEASNOISE: Compute the standard deviation of the measurement noise
%                assuming an N point FFT has been applied via Apply_FFT_N
%
% Inputs:
% 1) N           Lenth of the fft: N = 2K+1
% 2) STD_ns      Vector of standard deviations of the AWGN (White Noise)
%
% Outputs:
% 1) std_w       Structure of standard deviations of the real and imaginary
%                componenets of the fft'd AWGN at any frequnecy w (omega)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std_w.Vm = STD_ns.Vm*sqrt(2/N); 
std_w.Va = STD_ns.Va*sqrt(2/N); 
std_w.Im = STD_ns.Im*sqrt(2/N); 
std_w.Ia = STD_ns.Ia*sqrt(2/N); 

% Proof =>
% ns_std = 17.4;
% N      = 10001;
% trls   = 1000;
% val    = zeros(trls,1);
% for ii = 1:1000
%     test_ns = randn(N,1)*ns_std;
%     [~,y_fft] = Apply_FFT_N(test_ns,1,N);
%     val(ii) = real(y_fft(1000));
% end
% std_w_TRUE = ns_std*sqrt(2/N);
% std_w_EMPC = std(val);

end

