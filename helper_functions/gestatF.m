function [P,C,F,err2n,X,Y] = gestatF(x,y,nfft,Fsamp,Fmax,Fat)
% function [P,C,F,err2n,X,Y] = gestben(x,y,nfft,Fsamp,Fmax)
%
% modified from gest.m by Doug MacMartin (original date 9/28/15)
% estimates transfer function and spectral coherence between input and output signals
% all fitting and plotting functionality removed
%
% Inputs:
%   x     = input signal
%   y     = output signal
%   nfft  = number of points to use in fft (more means more resolution, less averaging)
%           possibly a vector (e.g. nfft=[8 16 32 64 128])
%   Fsamp = sample frequency
%   Fmax  = maximum real frequency to use in estimate (e.g., 6/year)
%           if length(Fmax)==2, then first element is minimum frequency
%
% Outputs:
%   P     = estimate of the transfer function between x and y
%   C     = estimate of the magnitude squared coherence
%   F     = frequencies at which P and C are calculated
%   err2n = variance of normalized error
%   X     = power spectral density of X
%   Y     = power spectral density of Y

if nargin<4
    Fsamp=1;
end
if isempty(Fsamp)
    Fsamp=1;
end
if nargin<5
    Fmax=inf;
end
if isempty(Fmax)
    Fmax=inf;
end
if length(Fmax)==2
    Fmin=Fmax(1);
    Fmax=Fmax(2);
else
    Fmin=eps;
end

% compute number of averages (for estimating error in g)
M=length(x);
nfft=sort(nfft);
for knf=1:length(nfft)
    L=nfft(knf);
    noverlap = fix(0.5.*L); % number of overlap, bs:fix rounds down to integer
    LminusOverlap = L-noverlap;
    k = fix((M-noverlap)/LminusOverlap); % or, k=1+(M-L)/(L-O) = (M-O)/(L-O)

    % test if anti-aliasing helps
    %[b,a]=butter(1,0.5);
    %x=filtfilt(b,a,x);
    %y=filtfilt(b,a,y);
    
    win=hamming(nfft(knf));
    [tmp(knf).P,tmp(knf).F]=tfestimate(x,y,win,[],nfft(knf),Fsamp,Fat);
    [tmp(knf).C,tmp(knf).F]=mscohere(x,y,win,[],nfft(knf),Fsamp,Fat);
    [tmp(knf).X,tmp(knf).F]=pwelch(x,win,[],nfft(knf),Fsamp,Fat);
    [tmp(knf).Y,tmp(knf).F]=pwelch(y,win,[],nfft(knf),Fsamp,Fat);
    if knf==1,P=tmp(1).P;F=tmp(1).F;C=tmp(1).C;X=tmp(1).X;Y=tmp(1).Y;
        err2n=(1-C)./C/(2*k);   % variance of normalized error
    else
        kf=2;       % bs: appends lower frequency information to bottom
        II=find(tmp(knf).F<F(kf));  % note F(1) is always zero
        P=[tmp(knf).P(II);P(kf:end)]; 
        F=[tmp(knf).F(II);F(kf:end)];
        C=[tmp(knf).C(II);C(kf:end)];
        err2n=[(1-C(II))./C(II)/(2*k);err2n(kf:end)];   % variance of normalized error
        X=[tmp(knf).X(II);X(kf:end)];
        Y=[tmp(knf).Y(II);Y(kf:end)];
    end
end
