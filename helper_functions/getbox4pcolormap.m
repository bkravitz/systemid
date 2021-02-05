function BOX = getbox4pcolormap(DataMatrix,lons,lats)

% function  BOX = getbox4pcolormap(DataMatrix,lons,lats)
% This function gets coordinate values of perturbed region from name of
% data matrix and creates a box struct to use with pcolormap to draw a box
% around perturbed region.
%
% Inputs:
%   DataMatrix = data to be put into the struct
%   lons       = array of grid longitude values
%   lats       = array of grid latitude values
%
% Outputs: 
%   BOX       = the structure to use with pcolormap
%
%
% System identification for teleconnections
% Copyright (C) 2020  Bethany Sutherland and Ben Kravitz
% Contact:  ben.kravitz.work@gmail.com
%
% Written by Bethany Sutherland
% Last updated 29 October 2020
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

% get box coordinates
li = strfind(DataMatrix,'LAT');
if DataMatrix(li-4) == 'n'; c1 = -1; else; c1 = 1; end
if DataMatrix(li+3) == 'n'; c2 = -1; else; c2 = 1; end
lattmp = [c1*str2num(DataMatrix(li-3:li-1)) c2*str2num(DataMatrix(li+4:li+6))];
li = strfind(DataMatrix,'LON');
lontmp = [str2num(DataMatrix(li-3:li-1)) str2num(DataMatrix(li+3:li+5))];

% convert defined perturbed area to actual area perturbed 
%   (grid edges halfway between lat/lon values and area perturbed are all 
%   lat/lon < or > but not equal to defined area. )

lonlow = min(lontmp); 
    I = find(lons>lonlow,1);
    lonlow = lons(I); 

lonhigh = max(lontmp);
    I = find(lons<lonhigh,1,'last');
    lonhigh = lons(I+1); % ADDING ONE BECAUSE AREA PLOTTED IS DEFINED BY BOTTOM EDGE OF BOX
    % CHECK LATER!!!!
    
latlow = min(lattmp);
    I = find(lats>latlow,1);
    latlow = lats(I); 
    
lathigh = max(lattmp);
    I = find(lats<lathigh,1,'last');
    lathigh = lats(I+1); 

BOX.lon = [lonlow lonhigh lonhigh lonlow lonlow];
BOX.lat = [lathigh lathigh latlow latlow lathigh];
fprintf('NOTE! DOUBLE CHECK HOW PCOLORM PLOTS THINGS!!!')
