function [RegMap,num_regions] = split_into_regions(lons,lats,regsz,PertMap,do_plot)
% SPLIT_INTO_REGIONS creates a map of different regions of a given size
%
% written by: Bethany Sutherland
% last edited: 1/25/2019
%
% Inputs:
%   lons    = longitudinal values of the model output
%   lats    = latitudinal values of the model output
%   regsz   = [lon lat] = number of grid points to use in each direction
%   PertMap = Map which indicates what grid points were perturbed. Should
%             be the same dimensions as the model
%   do_plot = Logical argument which determines if a figure of the regions
%             is created
%
% Outputs:
%
%   yRegMap = Map of what region each lat/lon grid point belongs to
%   num_regions = the total number of regions created


lonreg = [];
regnum = 1;
Mat = [];
while size(Mat,1)<length(lons)
    Mat = [Mat; regnum*ones(regsz)];
    regnum = regnum+1;
end
colblock = regnum - 1;

regnum = 0;
RegMap = [];
while size(RegMap,2)<length(lats)
    RegMap = [RegMap Mat + colblock*regnum];
    regnum = regnum +1;
end

% check that RegMap is correct size
RegMap = RegMap(1:length(lons),1:length(lats));

num_regions = max(max(RegMap));

% define colormap for region plot
Clrs = [1 1 0;      % yellow
    0 1 0;          % green
    1 0 1;          % magenta
    0 1 1;          % cyan
    1 0 0;          % red
    0 0 1           % blue
    0 0.5 0;        % army green
    1, 0.64, 0.31;  % orange 
    0.5, 0.18, 0.55;% purple   
    0.6, 1, 0;      % yellow-green
    0, 0.8, 0.22;   % turquoise
    0.4, 0.6, 0.2;  % maroon
    0.62, 0.71, 0.8;% slate gray
    0.93, 0.71, 0.71;% dusty rose
    0.36, 0.2, 0.09;% chocolate brown
    0.5, 0 1;       % bright purple
    ];          
    
% limit colors used to number of longitudinal regions
while colblock > length(Clrs); Clrs = [Clrs; 0.67*Clrs(end-15:end,:)]; end
Clrs = Clrs(1:colblock,:);

lp = 0;
Clrmap = [];
while size(Clrmap,1) < num_regions
    lp = lp + 1;
    Clrmap = [Clrmap; (1-0.1*lp)*Clrs];
end

Clrmap2 = [Clrmap; 1 1 1];
map_opts2 = {'showColorbar', false,'ColorMap',Clrmap2};
PertMap = PertMap + RegMap;
PertMap(PertMap>99) = num_regions+1;


if do_plot == true    
    figure
    pcolormap(lats, lons, PertMap ,map_opts2{:})
    title({'Regions'})
end

end

