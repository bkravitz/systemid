function [] = map1axis_2colorscales(Xx,Yy,Zz,Cc,Cc2,lbl,figname,box,clim,clim2,Fpert,discretize)
% MAP1AXIS_2COLORSCALES creates 2 2d plots for use in paper
% 
% 
% written by: Bethany Sutherland
% last edited: 5/19/2020
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

if nargin < 7; figname = []; end
if nargin < 8; box = []; end
if nargin <9; clim = []; end
if nargin <10; clim2 = []; end
if nargin <11 || isempty(Fpert); Fpert = 0; end
if nargin <12
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

% _top plot_

subplot(2,1,1)
plot(coast.long,coast.lat,'k'); hold on

scatter3(Xx,Yy,Zz,[],Cc)

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
if ~isempty(clim)
    caxis(clim); 
end
c.Label.String = 'period';
c.Ticks = sort([1/2 1/7 1/30 1/91 1/182 1/365 1/365/2 1/365/5 1/365/10]);
c.TickLabels = {'10yr','5yr','2yr','1yr','6mo','3mo','1mo','1wk','2d'};


% _ bottom plot_

subplot(2,1,2)
plot(coast.long,coast.lat,'k'); hold on

scatter3(Xx,Yy,Zz,[],Cc2)

if ~isempty(box)
    plot3(box.lon,box.lat,ones(size(box.lon))*max(Zz),'r','linewidth',1.5)
end
xlabel('longitude')
ylabel('latitude')
zlabel('frequency')
xlim([0 360])
ylim([-90 90])
c2 = colorbar;
set(gca,'colorscale','log')
if ~isempty(clim2)
    caxis(clim2); 
end
c2.Label.String = lbl;
c2.Ticks = [0.1 0.5 1:10 15 20];


% save figure
if ~isempty(figname)
    print(gcf,'-dpng','-r300',figname);
end


end

