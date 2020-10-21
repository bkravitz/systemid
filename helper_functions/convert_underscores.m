function textout = convert_underscores(textin)

% function textout = convert_underscores(textin)
% This function adds a \ before all underscores so that they will be
% displayed as underscores and not subscripts in matlab plot titles.
%
% Inputs:
%   textin  = input text with '_'
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
% Last updated 21 October 2020
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

k = strfind(textin,'_');
textout = textin;
for n = 1:length(k)
    textout(k(n)+1:end+1) = textout(k(n):end);
    textout(k(n)) = '\';
    k = strfind(textout,'_');
end
