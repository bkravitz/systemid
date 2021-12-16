function [reg_w_thisF, regsel] = regions2calculate(Fapprox, F, num_regions, significant_Fs, significant_Fs_list, regions_w_signal, reg_crit, pct)
% REGIONS2CALCULATE determines which regions meet the regional section
% criteria
%
% written by: Bethany Sutherland
% last edited: 4/29/2019
%
% Inputs:
%   Fapprox         = [low high] the approximate frequencies to use for filter
%   F               = array of the DFT frequency bins
%   num_regions     = the total number of regions globe has been seperated into
%   significant_Fs  = structure array listing all the significant
%                     frquencies for each region
%   significant_Fs_list = same as significant_Fs, but with all F in single
%                     list
%   regions_w_signal= list of regions which include some signal
%   reg_crit        = determines what criterion is used to decide which regions to
%             calculate
%               Region criteria           calculate projections in region...
%               reg_crit = 1;     % for all regions
%               reg_crit = 2;     % if all frequencies within F_range are significant
%               reg_crit = 3;     % if any frequencies within F_range are significant
%               reg_crit = 4;     % if a given percentage
%   pct             = threshold to use if using reg_crit = 4 (as a ratio)
%
% Outputs:
%
%   reg_w_thisF     = list of regions which meet the selection criteria for
%                     calculation
%   F_range         = The cutoff frequencies to be used for the bandpass filtering
%   regsel          = string indicating the reg_crit used (for use with plotting
%                     titles later


[~, indxlow] = min(abs(F-Fapprox(1)));
[~, indxhigh] = min(abs(F-Fapprox(2)));

Fsinrange = F(indxlow:indxhigh);

% determine regions which have this range
reg_w_thisF = [];

% for all regions
if reg_crit == 1 
    reg_w_thisF = [1:num_regions]';
    regsel = 'all regions';


% if all frequencies within F_range are significant
elseif reg_crit == 2
    for r = 1:numel(fieldnames(significant_Fs))
        fld = ['region' num2str(regions_w_signal(r))];
        
        if sum(ismember(Fsinrange,significant_Fs_list.(fld))) == length(Fsinrange)
            % do this section
            reg_w_thisF = [reg_w_thisF; regions_w_signal(r)];
        end
    end
    regsel = 'regs. with ALL F in range';


% if any frequencies within F_range are significant
elseif reg_crit == 3
    for r = 1:numel(fieldnames(significant_Fs))
        fld = ['region' num2str(regions_w_signal(r))];
        
        if sum(ismember(Fsinrange,significant_Fs_list.(fld))) > 0
            % do this section
            reg_w_thisF = [reg_w_thisF; regions_w_signal(r)];
        end
    end
    regsel = 'regs. with ANY F in range';

% if a percentage of frequencies within F_range are significant
elseif reg_crit == 4
    for r = 1:numel(fieldnames(significant_Fs))
        fld = ['region' num2str(regions_w_signal(r))];
        
        if sum(ismember(Fsinrange,significant_Fs_list.(fld)))/length(Fsinrange) > pct
            % do this section
            reg_w_thisF = [reg_w_thisF; regions_w_signal(r)];
        end
    end
    regsel = ['regs. with ' num2str(100*pct) '% F in range'];
end


end

