function [ax, varargout] = pcolormap(lat,lon,z,varargin)
% wrapper to pcolorm with typical values used in paper
persistent coast

opts = setOptions(varargin{:});

if isempty(coast)
    coast = load('coast.mat');
end

ax = setMapAGU();
[lonExt, zExt] = extendLon(lon, z);
pcolorm(lat, lonExt, zExt.'); hold on
setm(gca,'MeridianLabel','on','ParallelLabel','on','MLineLocation',60,'PLineLocation',[-60 -30 0 30 60],'MLabelLocation',60,'PLabelLocation',[-60 -30 0 30 60])

if opts.showCoast
    plotm(coast.lat,coast.long,'k');
end


colormap(ax, opts.ColorMap);
caxis(ax, opts.caxis);


if opts.showColorbar
    cbax = colorbar;
    varargout{1} = cbax;
end

if ~isempty(opts.box)
    % Add box around perturbed area
    blat = [max(opts.box.lat); max(opts.box.lat); min(opts.box.lat); min(opts.box.lat); max(opts.box.lat); ];
    blon = [min(opts.box.lon); max(opts.box.lon); max(opts.box.lon); min(opts.box.lon); min(opts.box.lon); ]; 
    plotm(blat,blon,'c','LineWidth',1.1)
end

if opts.PadRight ~= 0
    lonlim = getm(ax,'maplonlim');
    %     setm(ax,'maplonlim',[lonlim(1) lonlim(end) + opts.PadRight]);
    setm(ax,'maplonlim',[lonlim(1) lonlim(end)])
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
p.addParameter('ColorMap',parula)
p.addParameter('PadRight', -3);
p.addParameter('box',[])

% parse inputs
p.parse(varargin{:})
plotOptions = p.Results;
end