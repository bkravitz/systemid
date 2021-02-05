function dat = nc_varget(varargin)

% function  dat = nc_varget(varargin)
% This function adds a \ before all underscores so that they will be
% displayed as underscores and not subscripts in matlab plot titles.
%
% Inputs:
%   textin  = input text with '_' (you can have as many as you want)
%
% Outputs: 
%   textout = output text with '\_'
%
%
% System identification for teleconnections
% Copyright (C) 2020  Bethany Sutherland and Ben Kravitz
% Contact:  ben.kravitz.work@gmail.com
%
% Written by Bethany Sutherland
% Last updated 1 November 2020
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

% squeeze out singleton dimensions
dat = squeeze(ncread(varargin{:}));

% convert to MATLAB native double
if isnumeric(dat) && ~isa(dat,'double')
   dat = double(dat);
end
