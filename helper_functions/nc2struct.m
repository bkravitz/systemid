function [out, atts] = nc2struct(ncfiles, vars, newnames)
%NC2STRUCT reads data from netcdf file into a structure array
% From Chris Jones 2017
%
% out = NC2STRUCT(ncfiles): Reads all variables into structure array "out"
% out = NC2STRUCT(ncfiles, vars): Reads variables contained in cell array "vars"
%   into structure array "out."
% out = NC2STRUCT(ncfiles, vars, newnames): Renames vars(1:n) with
%   newnames in the output structure, where "n" is the length of the
%   newnames cell array.
% [out,atts] = NC2STRUCT(...): out contains variable content; atts
%   structure contains attributes of variables.
% 

ncfiles = parseFilenames(ncfiles);

nf = numel(ncfiles);
doReadAll = false;

if (~exist('vars','var') || isempty(vars))
   % get all variables in file:
   vars = getVars(ncfiles);
   doReadAll = true;
elseif ischar(vars)
   vars = {vars};
end

nv = numel(vars);

% customize varnames:
if ~exist('newnames','var')
    newVars = customvarname(vars);
else
    newVars = customvarname(vars, newnames);
end    

for f=nf:-1:1
   for vv=1:nv
       vin = vars{vv}; % input variable name
       vout = newVars{vv}; % output variable name
       out(f).(vout) = nc_varget(ncfiles{f},vin);
   end
   out(f).ncfile = ncfiles{f};
end

if nargout>1
   if doReadAll
      for f=nf:-1:1
         nctmp = ncinfo(ncfiles{f});
         v2 = {nctmp.Variables.Name};
         for j=1:nv
            v = vars{j};
            vout = newVars{j};
            ix = find(strcmp(v,v2));
            if isempty(ix)
               atts(f).(vout) = [];
            else
               attNames = {nctmp.Variables(ix).Attributes.Name};
               nAtt = numel(attNames);
               for a=1:nAtt
                  a1 = matlab.lang.makeValidName(attNames{a}); % sanitize field name
                  atts(f).(vout).(a1) = nctmp.Variables(ix).Attributes(a).Value;
               end
               % add dimensions/size as well:
               atts(f).(vout).Dimensions = {nctmp.Variables(ix).Dimensions.Name};
               atts(f).(vout).Size = nctmp.Variables(ix).Size;
            end
         end
      end
   else
      for f=nf:-1:1
         for j=1:nv
            v = vars{j};
            vout = newVars{j};
            try 
               nctmp = ncinfo(ncfiles{f}, v);
               attNames = {nctmp.Attributes.Name};
               nAtt = numel(attNames);
               for a=1:nAtt
                  a1 = matlab.lang.makeValidName(attNames{a}); % sanitize field name
                  atts(f).(vout).(a1) = nctmp.Attributes(a).Value;
               end
               % add dimensions as well:
               atts(f).(vout).Dimensions = {nctmp.Dimensions.Name};
               atts(f).(vout).Size = nctmp.Size;
            catch ME
               atts(f).(vout) = [];
            end
         end
      end
   end
end

end

function vars = getVars(ncfiles)
% gets list of all unique variables in files ncfiles
vars = {};
nf = numel(ncfiles);
   for f=nf:-1:1
      nctmp = ncinfo(ncfiles{f});
      v2 = {nctmp.Variables.Name};
      vars = unique([v2, vars]);
   end
end

function filelist = parseFilenames(ncfiles)
% ncfiles could be single file (string), directory structure, or list of
% filenames.
if ischar(ncfiles)
    % convert single filename to cell array
    filelist = {ncfiles};
elseif isstruct(ncfiles)
    % structure indicates ncfiles obtained from "dir"
    filelist = strcat({ncfiles.folder}, filesep, {ncfiles.name});
elseif iscell(ncfiles)
    filelist = ncfiles;
else
    error('nc2struct:ncfiles', 'unknown input file format')
end


end

function vout = customvarname(vin, vnew)
% replace offending characters in output file names
% drop characters before "-"
vtmp = vin;

if exist('vnew', 'var') && ~isempty(vnew)
    for vi=1:numel(vnew)
        vtmp{vi} = vnew{vi};
    end
end

% hasdash = contains(vin,'-');
% vtmp(hasdash) = cellstr(extractBefore(vin(hasdash),'-'));

vout = matlab.lang.makeValidName(vtmp);

end
