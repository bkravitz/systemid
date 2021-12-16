function [SNR,gammasq,indx] = SNR_from_coherence(C,F,Fi,Fs)
%SNR_FROM_COHERENCE estimates the signal to noise ration from the magnitude
% squared coherence
% 
% written by: Bethany Sutherland
% last edited: 11/6/2018
%
% Inputs:
%   C       = magnitude squared coherence
%   F       = frequencys at which C is caluclated
%   Fi      = frequency of signal you are looking for
%   Fs      = sampling frequency
%
% Outputs:    
%   SNR     = estimated signal to noise ratio

[tmp, indx] = min(abs(F-Fi/Fs));
gammasq = C(indx);  % use value closest to frequency of perturbation
SNR = gammasq/(1-gammasq); fprintf('using Fay form of the equation!\n')
% SNR = sqrt(gammasq/(1-gammasq)); fprintf('USING EQU 7.69 SWANSON)\n')

end

