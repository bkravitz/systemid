function [SNR] = snr_est_regression(y,Fi,Fs)
%SNR_EST_REGRESSION uses regression to estimate the signal to noise ratio
%   Assumes that the signal is a sine/cosine with frequency Fi
% 
% written by: Bethany Sutherland
% last edited: 11/6/2018
%
% Inputs:
%   y       = signal
%   Fi      = frequency of signal you are looking for
%   Fs      = sampling frequency
%
% Outputs:
%   SNR     = estimate of the signal to noise ratio    

yrT = 365;
t = 0:1/Fs:length(y)/Fs-1/Fs;
ptT = Fs/Fi;

A = [sin(2*pi*t/yrT)' cos(2*pi*t/yrT)' sin(2*pi*t/ptT)' cos(2*pi*t/ptT)' ones(size(t'))];

load = A\y';

bg_reg = load(1)*sin(2*pi*t/yrT)+load(2)*cos(2*pi*t/yrT);
pt_reg = load(3)*sin(2*pi*t/ptT)+load(4)*cos(2*pi*t/ptT);

% Calculate amplitude and phase of the recreated signal
SigA = sqrt(load(4)^2+(load(3))^2);
noise = y-pt_reg;
navg = mean(noise);
Sigma_sig = SigA/sqrt(2);
Sigma_noise = sqrt(sum((noise-navg).^2)/length(noise));
SNR = (Sigma_sig)^2/(Sigma_noise)^2;


end

