function [ghat,gerr,H,C,F,prj, prj_c, avgH] = tranFval3sig(x, y1, y2, y3, y4, cntrl,nfft,Fsamp,Fmax,ctoff,flag,fignum )
% TRANFVAL a wrapper for the code gest2.m (which was written and provided by 
% D. MacMartin)
%
% Compares a transfer function between input x and 3 different outputs y1,
% y2, y3
%
% written by: Bethany Sutherland
% last edited: 1/30/2019
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
%
%   C       = magnitude squared coherence between input and output signals
%   H       = Transfer function estimate
%   F       = Frequency of coherence estimate bin centers
%   lF      = list of frequencies that do not overlap with control 
%   xF      = same as lF but with each group of isolated frequencies in a
%               cell
%   prj     = projection of output onto input
%   prj_c   = projection of control onto input
%   avgH    = average magnitude of the transfer function
%   avgGreyH= average magnitude of the transfer function which is distince
%             from the control run
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



%% Transfer functions
figinfo = setfiginfo('colour','r','Magfig',fignum);
[ghat_c,gerr_c,H_c,C_c,F_c,he_c,hep_c,hmag_c,hph_c] = gest2(x,cntrl,nfft,Fsamp,Fmax,flag,figinfo);
figinfo.colour = 'b';
[ghat,gerr,H,C,F,he,hep,hmag,hph] = gest2(x,y1,nfft,Fsamp,Fmax,flag,figinfo);
figinfo.colour = 'g';
[ghat2,gerr2,H2,C2,F2,he2,hep2,hmag2,hph2] = gest2(x,y2,nfft,Fsamp,Fmax,flag,figinfo);
figinfo.colour = 'c';
[ghat2,gerr2,H2,C2,F2,he2,hep2,hmag2,hph2] = gest2(x,y3,nfft,Fsamp,Fmax,flag,figinfo);
figinfo.colour = 'y';
[ghat2,gerr2,H2,C2,F2,he2,hep2,hmag2,hph2] = gest2(x,y4,nfft,Fsamp,Fmax,flag,figinfo);


%% Determine Overlapping Areas
% 
% subplot(figinfo.nrow,figinfo.ncol,figinfo.subplt(1))
% yl = ylim;
% xl = xlim;
% 
% % determine where overlap does not occur
% freq2include = false(size(he.XData));
% freq2include(he.YData - he.YPositiveDelta > he_c.YData + he_c.YPositiveDelta) = true;
% freq2include(he.YData + he.YPositiveDelta < he_c.YData - he_c.YPositiveDelta) = true;
%         % Always use YPositiveDelta for analysis because YNegativeData
%         % is altered to prevent negative values to avoid plotting issues
% lF = he.XData(freq2include);
% 
% % plot gray squares behind non overlaping areas
% tmparray = [];
% celindx = 1;
% xF = {[]};
% freqlabels = [0 freq2include 0];
% for m = 2:length(freqlabels)
%     
%     findx = m-1;
%     fcg = freqlabels(m)-freqlabels(m-1);
%     if fcg == 1 % start of group
%         tmparray = [tmparray he.XData(findx)];
%     elseif fcg == 0 && freqlabels(m) ==1 % still in group
%         tmparray = [tmparray  he.XData(findx)];
%     elseif fcg == -1 % end of group
%         xF(celindx,:) = {tmparray};
%         tmparray = [];
%         celindx = celindx + 1;
%     elseif fcg == 0 && freqlabels(m) == 0 % outside of group
%         continue
%     end
% end
% 
% 
% for k1 = 1:size(xF,1)
%     q = xF{k1};
%     
%     % add overhanging area hafway to adjacent frequencies
%     if length(q)>=1
%         Ldf = sort(abs(he.XData - q(1)));
%         LOh = q(1) - Ldf(2)/2;
%         Rdf = sort(abs(he.XData - q(end)));
%         ROh = q(end) + Rdf(2)/2;
%         
%         q = [LOh q ROh];
%     end
%     
%     qx = [min(q) max(q)  max(q)  min(q)];
%     
%     if length(qx)>1
%         qy = [[1 1]*yl(1) [1 1]*yl(2)];
%         
%         p = patch(qx, qy, [1 1 1]*0.8);
%         % make boxes semi-transparent
%         p.FaceAlpha = 0.35; 
%         p.EdgeAlpha = 0.35;
%         % ensures addition of box does not change scale of axis
%         xlim(xl) 
%         ylim(yl)
%     end
% end

%% Calculate Projections

% calculate original projections
prj = sum(y1.*x)/sum(x.^2); % the sensitivity is just the projection of the output onto the input
prj2 = sum(y2.*x)/sum(x.^2);
prj3 = sum(y3.*x)/sum(x.^2);
prj4 = sum(y4.*x)/sum(x.^2);
prj_c = sum(cntrl.*x)/sum(x.^2);

avgH = mean(abs(H));
% avgGreyH = mean(abs(H(freq2include)));

%% Plotting

plotTitle = @(p1, p2, p3, p4, p5,ct ) {...
        [num2str(ct) '% lag ' ]; ...    
    ['Projection onto control (noise) = ' num2str(p1)]; ...    
    ['Projection onto x + noise = ' num2str(p2)]; ...    
    ['Projection onto x time delay removed = ' num2str(p3)]; ...
    ['Projection onto x time delayed = ' num2str(p4)]; ...
    ['Projection onto x phase matched = ' num2str(p5)]; ...        
    [' ']};     % add extra space between title and plots

subplot(figinfo.nrow,figinfo.ncol,figinfo.subplt(end))
legend('control','x+noise','x time delay removed','x time delayed','phase matched')

subplot(figinfo.nrow,figinfo.ncol,figinfo.subplt(1))
title(plotTitle( prj_c, prj, prj2, prj3, prj4 ,ctoff));




end

