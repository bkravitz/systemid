function [F_range] = F2use(Fapprox, F)
% REGIONS2CALCULATE alters the phase of the output, y, to match the phase of the 
% input, x.
%
% written by: Bethany Sutherland
% last edited: 1/22/2019
%
% Inputs:
%   Fapprox         = [low high] the approximate frequencies to use for filter
%   F               = array of the DFT frequency bins


[~, indxlow] = min(abs(F-Fapprox(1)));
[~, indxhigh] = min(abs(F-Fapprox(2)));

F_range = [F(indxlow) F(indxhigh)];