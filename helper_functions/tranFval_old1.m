function [] = tranFval(x,y,cntrl,nfft,Fsamp,Fmax,flag,fignum )
% TRANFVAL a wrapper for the code gest2.m (written and provided by D.
% MacMartin)
%
% Compares a transfer function between input x and output y
% with a transfer function between input x and a control signal to
% determine regions where transfer function could be caused by natural
% variability
%
% written by: Bethany Sutherland
% last edited: 12/19/2018
%
% Inputs:
%   x       = input signal
%   y       = output signal
%   cntrl   = noise estimate from control run
%   nfft    = number of sampling points used to calculate the DFT
%   Fsamp   = sampling frequency
%   Fmax    = maximum real frequency to use in estimate (e.g. 6/year)
%               if length(Fmax)=2, then first element is minimum freq
%   flag    = define fitting (see options below)
%   figinfo = figure plotting parameters
%               for default values use setfiginfo
%               individual parameters can be changed (e.g. setfiginfo('Phfig',105)))
%
% Outputs:
%   Sxx     = Auto-correlation of input signal
%   Syy     = auto-correlation of output signal
%   Sxy     = correlation between input and output signals
%   C       = magnitude squared coherence between input and output signals
%   H       = Transfer function estimate
%   F       = Frequency of coherence estimate bin centers
%
%
% flag=-1, then skip fitting
% flag=0 then fit a constant
% flag=1 then multiply xfer by s before averaging (i.e. fit an integrator)
% flag=2 then fit mu/(s+epsilon) (and ghat has two entries for mu, epsilon)
% flag=3, fit identical to flag=2, but reported as g/(1+tau s)
% flag=4, fit mu/(1+2*zeta*(s/om)+(s/om)^2)  (DGM 12/21/12)
% %%% now part of figinfo!! pltstr=0 then no plot is generated, otherwise pltstr is the string label
% optional fitval over-rides estimates for fit to use on plot
% figinfo provides additional arguments to use to control plotting
% of multiple transfer functions... figure # for magnitude, phase, # of
% subplot columns & rows, subplot number, and font size


%% Adjustable Parameters


%% Set Defaults for missing arguments

%% Transfer functions
figinfo = setfiginfo('colour','b','Magfig',fignum);
[ghat,gerr,P,C,F,he,hep,hmag,hph] = gest2(x,y,nfft,Fsamp,Fmax,flag,figinfo);
figinfo.colour = 'r';
[ghat_c,gerr_c,P_c,C_c,F_c,he_c,hep_c,hmag_c,hph_c] = gest2(x,cntrl,nfft,Fsamp,Fmax,flag,figinfo);


prj = sum(y.*x)/sum(x.^2); % the sensitivity is just the projection of the output onto the input
prj_c = sum(cntrl.*x)/sum(x.^2);

%% Determine Overlapping Areas

subplot(figinfo.nrow,figinfo.ncol,figinfo.subplt(1))
yl = ylim;

% determine where overlap does not occur
freq2include = false(size(he.XData));
freq2include(he.YData - he.YPositiveDelta > he_c.YData + he_c.YPositiveDelta) = true;
        % Always use YPositiveDelta for analysis because YNegativeData
        % is altered to prevent negative values to avoid plotting issues
listfreq = he.XData(freq2include);

% plot gray squares behind non overlaping areas
tmparray = [];
celindx = 1;
xF = {[]};
for findx = 2:length(freq2include)
    findx
    fcg = freq2include(findx)-freq2include(findx-1)
    
    if sum(freq2include) == length(freq2include) % all 1s
        fprintf('loop 1 \n')
        xF = {he.XData};
        break
    elseif sum(freq2include) == 0 % all 0s
        fprintf('loop 2 \n')
        xF = {[]};
        break
    end
    if fcg == 0 && freq2include(findx) == 0
        % no transition
        fprintf('loop 3 \n')
        continue
    elseif fcg == 1
        % fi is beginning of group
        fprintf('loop 4 \n')
        tmparray = [tmparray he.XData(findx)];
    elseif fcg == 0 && freq2include(findx) == 1
        % fi is also in group
        fprintf('loop 5 \n')
        tmparray = [tmparray he.XData(findx)];
        if findx == length(freq2include)
            % end of array, include last group
            fprintf('loop 6 \n')
            xF(celindx,:) = {tmparray};
        elseif fcg == -1
            % group has ended
            fprintf('loop 7 \n')
            xF(celindx,:) = {tmparray};
            celindx = celindx + 1;
            tmparray = [];
        end
    end
end


for k1 = 1:size(xF,1)
    q = xF{k1};
    
    % add overhanging area hafway to adjacent frequencies
    if length(q)>=1
        Ldf = sort(abs(he.XData - q(1)));
        LOh = q(1) - Ldf(2)/2;
        Rdf = sort(abs(he.XData - q(end)));
        ROh = q(end) + Rdf(2)/2;
        
        q = [LOh q ROh];
    end
    
    qx = [min(q) max(q)  max(q)  min(q)];
    
    if length(qx)>1
        qy = [[1 1]*yl(1) [1 1]*yl(2)];
        
        p = patch(qx, qy, [1 1 1]*0.8);
        p.FaceAlpha = 0.5;
        p.EdgeAlpha = 0.5;
    end
end



%% Plotting

plotTitle = @(p1, p2 ) {...
    ['Projection onto signal = ' num2str(p1)]; ...
    ['Projection onto control = ' num2str(p2)]; ...
    [' ']};     % add extra space between title and plots

subplot(figinfo.nrow,figinfo.ncol,figinfo.subplt(end))
legend('model','control')

subplot(figinfo.nrow,figinfo.ncol,figinfo.subplt(1))
title(plotTitle(prj, prj_c));



end

