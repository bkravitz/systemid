function [lagout,Corrout] = lagcorrect(a,b,varargin)

% function  [lagout,Corrout] = lagcorrect(a,b,[maxlag,doplot,tm])
% This function looks at correlation coefficients between two timeseries
% and estimates the amount of lag between the input and output.
%
% Inputs:
%   a       = input signal
%   b       = output signal
%   maxlag  = the maximum allowable lag between the input and output, if
%               not specified maxlag = 1/2 the size of the timeseries. 
%   tm      = time vector
%
% Outputs:
%
%   lagout  = estimated lag between input and output
%               [first local maximum; 
%                first local max outside 1 standard deviation;
%                first local max outside 2 standard deviations]
%   Corrout = correlation between signals after x is shifted and truncated 
%               by lagout timesteps
%
%
% System identification for teleconnections
% Copyright (C) 2020  Bethany Sutherland and Ben Kravitz
% Contact:  ben.kravitz.work@gmail.com
%
% Written by Bethany Sutherland
% Last updated 6 November 2020
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% setting default values for variable arguments
maxlag = floor(length(tm)/2);
doplot = false;
tm = 0:length(a)-1;

if isempty(varargin)==0;
    maxlag = varargin{1};
    if length(varargin)>1;
        doplot = varargin{2};
    end
    if length(varargin)>2;
        tm = varargin{3};
    end
end

if size(a,1) >1; a = a'; end
if size(b,1) >1; b = b'; end

if doplot
    figure
end

Cor = nan(1,maxlag);

for ni = 0:maxlag
    n = ni+1;
    atmp = [nan(1,ni) a(1:end-ni)];
    
    Corno0tmp = corrcoef(atmp(~isnan(atmp)),b(~isnan(atmp)));
    Cor(n) = Corno0tmp(1,2);
    
    if doplot
        subplot(3,1,1)
        plot(tm,atmp,tm,b)
        legend('NI','PI')
        title('timeseries')
        
        subplot(3,1,2)
        plot(tm,Cor)
        title('correlation')
        
        pause(0.05)
    end
end

sig = std(Cor);

[~,pks] = findpeaks(abs([0 Cor])); % add 0 at start to ensure no lag also gets checks as a peak
pks = pks - 1; % adjust for having added the 0

% find first peak outside 1 standard deviation
pk1sig = find(abs(Cor(pks))-sig>0,1,'first');
lag1sig = pks(pk1sig)-1 ; % -1 to account for starting at 0 lag
if isempty(pk1sig); lag1sig = NaN; end

% find first peak outside 1.5 standard deviation
pk15sig = find(abs(Cor(pks))-1.5*sig>0,1,'first');
lag15sig = pks(pk15sig)-1 ; % -1 to account for starting at 0 lag
if isempty(pk15sig); lag15sig = NaN; end

% find first peak outside 2 standard deviations
pk2sig = find(abs(Cor(pks))-2*sig>0,1,'first')  ; 
lag2sig = pks(pk2sig) - 1;
if isempty(pk2sig); lag2sig = NaN; end

% find first peak with no qualification
lag0sig = pks(1)-1;
lagout = [lag0sig; lag1sig; lag15sig; lag2sig];

if doplot
    subplot(3,1,3)
    plot(tm, [nan(1,lagout) a(1:end-lagout)],tm,b)
    title( {['maximum correlation lag= ' num2str(lagout)],...
        [ ' corr no0 = ' num2str(corrout)]})
    figure;
    plot(0:maxlag,Cor,0:maxlag,abs(Cor),'--',[0 maxlag],[sig sig],[0 maxlag],1.5*[sig sig],[0 maxlag],2*[sig sig])
    legend('Cor','|Cor|','1 \sigma','2 \sigma')

end
Corrout = Cor;
