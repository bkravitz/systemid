function [C,f,Sxx,Syy,Sxy,P1X,P1Y,X_k,Y_k] = mscest(x,y,nfft,Fs,segSz,window,OL)
% MSCEST estimates the magnitude squared coherence between two
% signals
%
% written by: Bethany Sutherland
% last edited: 11/13/2018
%
% Inputs:
%   x       = input signal
%   y       = output signal
%   nfft    = number of sampling points used to calculate the DFT
%   Fs      = sampling frequency
%   SegSz   = number of points in each individual segment
%   window  = window to apply to each segment
%
% Outputs:
%   Sxx     = Auto-correlation of input signal
%   Syy     = auto-correlation of output signal
%   Sxy     = correlation between input and output signals
%   C       = magnitude squared coherence between input and output signals
%   H       = Transfer function estimate
%   F       = Frequency of coherence estimate bin centers 

if nfft<segSz
    fprintf('warning from mscohere_est: number of points in fft is less than segment size, segments are being truncated to nfft\n')
end

% determine correct frequencies
fprintf('check frequencies, maybe should be like in removingphase.m')
dF = Fs/nfft;
f = (0:(nfft-1))*dF; 
if mod(nfft,2) == 0 % if nfft == even 
    nyq = f(end/2+1);
else
    nyq = f(ceil(end/2))+dF/2;
end
f(f>nyq) = f(f>nyq) - (nyq*2);
F = f(f>0);

% estimate auto and cross correlations through periodograms
overlap = floor(OL*segSz); % number of overlap
if overlap ~= OL*segSz
    fprintf(['cannot do exactly ' num2str(OL*100) '%% overlap, rounding to ' num2str(overlap*100/segSz) '%%\n'])
end
leftover = segSz-overlap;
fst = 1;  si = 1; 
lst = fst + segSz -1;
while lst(si) < length(x)
    si = si+1;
    fst(si) = fst(si-1) + leftover;
    lst(si) = fst(si) + segSz -1;
end

if length(fst)>1
    fst = fst(1:end-1);
    lst = lst(1:end-1);
end

for k = 1:length(fst)
    
    % extract individual segment
    whichtimes=fst(k):lst(k);
    x_k_t = x(whichtimes);            % extract segment,
    x_k_t = x_k_t.*window';                  % apply window
    y_k_t = y(whichtimes);            % extract segment,
    y_k_t = y_k_t.*window';                  % apply window
    
    % take Fourier Transform
    X_k(:,k) = fft(x_k_t,nfft)/segSz;              % compute the fourier transform   
    Y_k(:,k) = fft(y_k_t,nfft)/segSz;              % compute the fourier transform
    
    % compute the two sided power spectrum
    P2X(:,k) = abs(X_k(:,k));
    P2Y(:,k) = abs(Y_k(:,k));
    
    % compute the single sided power spectrum
    P1X(:,k) = P2X(f>=0,k);
    P1X(2:ceil(nfft/2),k) = 2*P1X(2:ceil(nfft/2),k);
    P1Y(:,k) = P2Y(f>=0,k);
    P1Y(2:ceil(nfft/2),k) = 2*P1Y(2:ceil(nfft/2),k);
    
end

Sxy = mean(conj(X_k).*Y_k,2);
Sxx = mean(conj(X_k).*X_k,2);
Syy = mean(conj(Y_k).*Y_k,2);

C = abs(Sxy).^2./(Sxx.*Syy);
end

