function [epsilon,snreb]=ceb(C,nd)

% function out=ceb(in)
% Computing magnitude squared coherence error bars using a simple error formula
% Written by Ben Kravitz (ben.kravitz@pnnl.gov or ben.kravitz.work@gmail.com)
% Last updated 21 November 2018
%
% Inputs
%   gamma2 = magnitude-squared coherence (at each frequency)
%   nd = number of independent segments (total series length divided by segment length, rounded down)
%
% Outputs
%   epsilon = random error estimate (at each frequency)
%
% 95 percent confidence interval error bars are [gamma2*(1-2*epsilon),gamma2*(1+2*epsilon)]

epsilon = (1-C).*sqrt(2/nd./C);

%  Added by Bethany Sutherland (11/21/18)
dsnr = sqrt(8*C/nd./(1-C).^2);

snreb = dsnr;%10*log10(dsnr); % 

end