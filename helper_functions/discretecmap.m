function [DataOut,cmap,clim,cticks,clabels] = discretecmap(Data,lvls,bifurcate,forceLvls,bluered )
% DiSCRETECMAP creates a discrete colormap appropriate for the given levels
% and maps in Cdata to appropriate values for the cooresponding colormap
%
% written by: Bethany Sutherland
% last edited: 7/8/2019
%
% Inputs:
%   Data   = data to be put into discrete color bins
%   levels  = values at which to separate each bin
%   bifurcate   = (optional) data below bifurcate will be plotted in black
%                   and white, data above will be plotted in color

% Outputs: 
%   DataOut    = Cdata mapped to descrete bins according to levels 
%   cmap    = 


if nargin <3; bifurcate = {}; end
if nargin <4; forceLvls = false; end
if nargin <5; bluered = false; end

DataOut = NaN(size(Data));
lvls = sort(lvls);
origlvls = lvls;

if ~forceLvls
    % remove levels outside of range of data
    Mx = max(max(Data)); Mn = min(min(Data));
    lvls = lvls(lvls>Mn & lvls<Mx);
    % check if bifurcation point is a set level
    if ~isempty(bifurcate)
        if ~ismember(bifurcate,lvls)
            lvls = [lvls bifurcate];
            lvls = sort(lvls);
        end
    end
end
if isempty(lvls) % correct for if min and max values both fall between two levels
    lvls = [origlvls(find(origlvls<Mn,1,'last')) origlvls(find(origlvls>Mx,1,'first'))];
end

lvls = [-Inf lvls Inf];
LvlsOut = lvls;

numClrs = length(lvls)-1; % number of colors to plot

for k = 1:numClrs
    DataOut(Data>lvls(k) & Data<=lvls(k+1)) = k;
end


% design colormap
if ~isempty(bifurcate)
    % check if bifurcation point is in range of data
    if ismember(bifurcate,lvls)
        bw = max(find(lvls <= bifurcate));
        if bluered
            bw = bw-1;
        end
    else
        bw = 0;
    end
    
else
    bw = 0;
end


if bluered
    % create blue to white to red colormap
    red_top     = [1 0 0];
    white_middle= [1 1 1];
    blue_bottom = [0 0 1];
    numClrs1side = max([bw numClrs-bw]);
    for ci = 1:3
        redcmap(:,ci) = linspace(white_middle(ci),red_top(ci),numClrs1side);
        bluecmap(:,ci) = linspace(blue_bottom(ci),white_middle(ci),numClrs1side);
    end
    cmap(1:bw,:) = bluecmap(end-bw+1:end,:);
    cmap(bw+1:numClrs,:) = redcmap(1:numClrs-bw,:);
else
    % use black and white and parula colormap
    tmpcmap = flipud(gray(bw*3));
    cmap(1:bw,:) = tmpcmap(2:2+bw-1,:);
    cmap(bw+1:numClrs,:) = parula(numClrs-bw);
end

clim = [0.5 numClrs+0.5];
cticks = [1:numClrs-1]*(clim(2) - clim(1))/numClrs + clim(1);
clabels = lvls(2:end-1);

% ensure more than one tick mark
if numClrs == 2 
    cticks(2) = cticks(1);
    cticks(1) = 0.5; cticks(3) = 2.5;
    I = find(origlvls==lvls(2),1);
    clabels = origlvls(I-1:I+1);
end



end

