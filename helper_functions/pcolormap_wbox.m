function ax = pcolormap(lat,lon,box,z,varargin)
% wrapper to pcolorm with typical values used in paper
persistent coast

opts = setOptions(varargin{:});

if isempty(coast)
    coast = load('coast.mat');
end

ax = setMapAGU();
[lonExt, zExt] = extendLon(lon, z);
pcolorm(lat, lonExt, zExt.'); hold on

if opts.showCoast
    plotm(coast.lat,coast.long,'k');
end


colormap(ax, opts.ColorMap);
caxis(ax, opts.caxis);


if opts.showColorbar
    colorbar
end

if opts.showBox
    % Add box around subgrid area
    plotm(box.lat,box.lon,'c','LineWidth',1.2)
end

if opts.PadRight ~= 0
    lonlim = getm(ax,'maplonlim');
    setm(ax,'maplonlim',[lonlim(1) lonlim(end) + opts.PadRight]);
end

tightmap;
end

function [lonExtended, zExtended] = extendLon(lon, z)
% wrap the longitude to get rid of ugly white stripe
lonExtended = [lon(:); lon(1)]; 
zExtended = [z; z(1,:)];
end

function plotOptions = setOptions(varargin)
% try doing this using the input parser:
p = inputParser;

% default options:
p.addParameter('caxis','auto')
p.addParameter('showCoast',true,@islogical)
p.addParameter('showColorbar',true,@islogical)
p.addParameter('showBox',true,@islogical)
p.addParameter('ColorMap',parula)
p.addParameter('PadRight', -3);


% parse inputs
p.parse(varargin{:})
plotOptions = p.Results;

end