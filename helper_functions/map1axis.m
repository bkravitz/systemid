function [] = map1axis(Xx,Yy,Zz,Cc,lbl,figname,box,clim,Fpert,discretize)
% MAP4AXIS creates a 3d plot and 3 axis projections onto each plane in
% one figure
% 
% 
% written by: Bethany Sutherland
% last edited: 6/14/2019
%
% Inputs:
%   X       = x axis (longitude) values for each point to be plotted
%   Y       = y axis (latitude) values for each point to be plotted
%   Z       = z axis (frequency) values for each point to be plotted
%   C       = data for coloring of each point
%   ttl     = string for title of the figure
%   lbl     = label for colorbar
%   figname = string of name to use when saving figure. Does not save
%               figure if figname is not specified
%   box     = include box info (from getbox4pcolormap()) to plot box around
%             a specific area
%   discretize = if you want to plot with a discrete colormap include
%               desired levels here
%

if nargin < 6; figname = []; end
if nargin < 7; box = []; end
if nargin <8; clim = []; end
if nargin <9 || isempty(Fpert); Fpert = 0; end
if nargin <10
    discretize = [];
    cmap = [];
    cticks = [];
    clabels = [];
end
% convert frequency to period
ZZ = Zz./(2*pi);
% load map plotting data
if ~exist('coast','var')
    coast = load('coast.mat');
end
% adjust coast longitudes from east west to east
coast.long(coast.long<0) = coast.long(coast.long<0)+ 360;
for ci = 1:length(coast.long)
    if abs(coast.long(ci+1) - coast.long(ci))>355
        coast.long(ci+2:end+1) = coast.long(ci+1:end);
        coast.long(ci+1) = NaN;
        coast.lat(ci+2:end+1) = coast.lat(ci+1:end);
        coast.lat(ci+1) = NaN;
    end
end
    set(gca,'colorscale','log')


figure('pos',[1 1 500 800])

subplot(2,1,2)
plot(coast.long,coast.lat,'k'); hold on
if isempty(discretize)
    scatter3(Xx,Yy,Zz,[],Cc)
else
    [Cc,cmap,clim,cticks,clabels] = discretecmap(Cc,discretize);
    scatter3(Xx,Yy,Zz,[],Cc)
%     scatter3(Xx(Cc>discretize),Yy(Cc>discretize),Zz(Cc>discretize),[],Cc(Cc>discretize))
%     scatter3(Xx(Cc<discretize),Yy(Cc<discretize),Zz(Cc<discretize),[],Cc(Cc<discretize),'x')
end
if ~isempty(box)
    plot3(box.lon,box.lat,ones(size(box.lon))*max(Zz),'r','linewidth',1.5)
end
xlabel('longitude')
ylabel('latitude')
zlabel('frequency')
xlim([0 360])
ylim([-90 90])
c = colorbar;
set(gca,'colorscale','log')
% c.Ticks = [0 0.1 0.5 1 1.5 2 5 10 20]
% subplot(2,2,2)
% plot(coast.long,coast.lat,'k'); hold on
% scatter3(Xx,Yy,Zz,[],Cc)
% %     plot3([360 360],[-90 90], [Fpert Fpert],'r','linewidth',1.5)
% plot3(ones(size(box.lon))*max(Zz),box.lat,ones(size(box.lon))*Fpert,'r','linewidth',1.5)
% xlabel('longitude')
% ylabel('latitude')
% zlabel('frequency')
% if ~isempty(clim); caxis(clim); end
% view(90,0)
% hold off
% 
% subplot(2,2,3)
% plot(coast.long,coast.lat,'k'); hold on
% scatter3(Xx,Yy,Zz,[],Cc)
% plot3(box.lon,ones(size(box.lat))*min(Zz),ones(size(box.lon))*Fpert,'r','linewidth',1.5)
% %     plot3([0 360],[-90 -90], [Fpert Fpert],'r','linewidth',1.5)
% xlabel('longitude')
% ylabel('latitude')
% zlabel('frequency')
% if ~isempty(clim); caxis(clim); end
% view(0,0)
% hold off

if ~isempty(clim);
    caxis(clim); 
end
if ~isempty(cmap); colormap(cmap); end

if ~isempty(cticks); c.Ticks = cticks; end
if ~isempty(clabels); c.TickLabels = clabels; end

subplot(2,1,1)
plot3(coast.long,coast.lat,ones(size(coast.long))/365/5,'k'); hold on
scatter3(Xx,Yy,Zz,[],Cc)
% plot3([0 360 360],[90 90 -90], [Fpert Fpert Fpert],'r','linewidth',1.5)
plot3(box.lon,box.lat,ones(size(box.lon))*Fpert,'r','linewidth',1.5)
xlabel('longitude')
ylabel('latitude')
zlabel('frequency')
if ~isempty(clim);
    caxis(clim); 
end
view(3)
hold off
% xlim([0 360])
% ylim([-90 90])
set(gca,'ZScale','log',...
    'ZTick',sort([1/2 1/7 1/30 1/91 1/182 1/365 1/365/2 1/365/5 1/365/10]),...
    'ZTickLabel',{'10yr','5yr','2yr','1yr','6mo','3mo','1mo','1wk','2d'},...
    'xlim',[0 360],'ylim',[-90 90],'zlim',[1/365/5 1/2],'colorscale','log');



c.Label.String = lbl;
c.Ticks = [0.1 0.5 1:10 15 20]

% save figure
if ~isempty(figname)
    print(gcf,'-dpng','-r300',figname);
end

end

