function [f_vec,y_fft] = MAP_FFT(signal,tstep,N)
% MAP_FFT: Apply an N point fft to signal with delta_t = t_step
%
% Inputs:
% 1) signal  Time domain signal
% 2) tstep   Time step between signal data points
% 3) N       N point fft (dtft) - must be odd
%
% Output:
% 1) f_vec   Frequency vector (Hz)
% 2) y_fft   Raw fft data (scaled and divided properly)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure that N is odd
if mod(N,2) == 0
    N = N-1;
end

% Now, take the N point FFT
y     = fft(signal/N,N);
fs    = 1/tstep;
f_vec = fs*(0:((N-1)/2))/N;

% Analyze Signal
y_fft = y(1:((N+1)/2));

% Double all but DC
y_fft(2:end) = 2*y_fft(2:end);
end

