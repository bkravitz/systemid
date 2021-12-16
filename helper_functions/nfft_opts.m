function [nfft] = nfft_opts(N,smp,T,fs)
%Gives you the the larges possible nfft to use which will give you an output
% bin centered on the desired frequency, which still having the desired
% number of samples
%
% Inputs:
%   N       = length of signal
%   smp     = desired number of samples
%   T       = period you care most about
%   fs      = sample frequency
%
% Outputs:
%   nfft     = best nfft to use



% nfft = 2N/(smp+1) assuming overlap of 50%
maxnfft = 2*N/(smp+1);

% f(m) = m*fs/N
% need largest N which creates a integer
% N = m*fs*T
n = 0; m = 0;
while n < maxnfft
    m = m + 1;
    n(m) = m*fs*T;
end
nfft = n(1:m-1);

end

