function [ghat,gerr,H,C,F,lF,xF] = tranFval_withphaseshifter(x,y,cntrl,nfft,Fsamp,Fmax,flag,fignum )
% TRANFVAL a wrapper for the code gest2.m (which was written and provided by 
% D. MacMartin)
%
% Compares a transfer function between input x and output y
% with a transfer function between input x and a control signal to
% determine regions where transfer function could be caused by natural
% variability
%
% written by: Bethany Sutherland
% last edited: 1/4/2019
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
figinfo = setfiginfo('colour','b','Magfig',fignum);
[ghat,gerr,H,C,F,he,hep,hmag,hph] = gest2(x,y,nfft,Fsamp,Fmax,flag,figinfo);
figinfo.colour = 'r';
[ghat_c,gerr_c,H_c,C_c,F_c,he_c,hep_c,hmag_c,hph_c] = gest2(x,cntrl,nfft,Fsamp,Fmax,flag,figinfo);


%% Determine Overlapping Areas

subplot(figinfo.nrow,figinfo.ncol,figinfo.subplt(1))
yl = ylim;
xl = xlim;

% determine where overlap does not occur
freq2include = false(size(he.XData));
freq2include(he.YData - he.YPositiveDelta > he_c.YData + he_c.YPositiveDelta) = true;
freq2include(he.YData + he.YPositiveDelta < he_c.YData - he_c.YPositiveDelta) = true;
        % Always use YPositiveDelta for analysis because YNegativeData
        % is altered to prevent negative values to avoid plotting issues
lF = he.XData(freq2include);

% plot gray squares behind non overlaping areas
tmparray = [];
celindx = 1;
xF = {[]};
freqlabels = [0 freq2include 0];
for m = 2:length(freqlabels)
    
    findx = m-1;
    fcg = freqlabels(m)-freqlabels(m-1);
    if fcg == 1 % start of group
        tmparray = [tmparray he.XData(findx)];
    elseif fcg == 0 && freqlabels(m) ==1 % still in group
        tmparray = [tmparray  he.XData(findx)];
    elseif fcg == -1 % end of group
        xF(celindx,:) = {tmparray};
        tmparray = [];
        celindx = celindx + 1;
    elseif fcg == 0 && freqlabels(m) == 0 % outside of group
        continue
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
        % make boxes semi-transparent
        p.FaceAlpha = 0.35; 
        p.EdgeAlpha = 0.35;
        % ensures addition of box does not change scale of axis
        xlim(xl) 
        ylim(yl)
    end
end

%% Remove Phase and Do Projections

% calculate original projections
prj = sum(y.*x)/sum(x.^2); % the sensitivity is just the projection of the output onto the input
prj_c = sum(cntrl.*x)/sum(x.^2);

% remove phase
 nfft = 512; segSz = nfft; Fs = 1;
window = hamming(nfft);
OL = 0.5; fprintf('un hard code segSz, nfft, ect. later')

[C2,f2,Sxx2,Syy2,Sxy2,P1X2,P1Y2,X_k,Y_k] = mscest(x,y,nfft,Fs,segSz,window,OL);

X = mean(X_k,2); Y = mean(Y_k,2); 
PhiX = angle(X);

Ymag = abs(Y);
PhiYnew = PhiX;

IFinal = Ymag .* sin(PhiYnew);
RFinal = Ymag .* cos(PhiYnew);
spec= RFinal + 1i*IFinal;

x2 = ifft(X,'symmetric');
y2 = ifft(spec,'symmetric'); %Invert the transform

figinfo.colour = 'c';
[ghat,gerr,H,C,F,he,hep,hmag,hph] = gest2(x2,y2,nfft,Fsamp,Fmax,flag,figinfo);
prj_phrmvd = sum(y2.*x2)/sum(x2.^2); % the sensitivity is just the projection of the output onto the input

 % compute the two sided power spectrum
    P2Y = abs(spec);
    
    % compute the single sided power spectrum
    P1Y = P2Y(f2>=0);
    P1Y(2:ceil(nfft/2)) = 2*P1Y(2:ceil(nfft/2));
    
%% Plotting

plotTitle = @(p1, p2, p3, p4 ) {...
    ['Projection onto signal = ' num2str(p1)]; ...
    ['Projection onto control = ' num2str(p2)]; ...
    ['Phase removed projection onto signal = ' num2str(p3)]; ...
    ['Avg. TransFun. Magnitude = ' num2str(mean(abs(H)))];...
    ['Avg. Grey TF. Magnitude = ' num2str(mean(abs(H(freq2include))))];...
    [' ']};     % add extra space between title and plots

subplot(figinfo.nrow,figinfo.ncol,figinfo.subplt(end))
% legend('model','control')
legend('model','control','phase adjusted')

subplot(figinfo.nrow,figinfo.ncol,figinfo.subplt(1))
title(plotTitle(prj, prj_c, prj_phrmvd));

figure('rend','painters','pos',[10 10 1500 750])
subplot(2,1,1)
plot(f2(f2>0),real(Y(f2>0)),'r',f2(f2>0),imag(Y(f2>0)),'r:',...
    f2(f2>0),real(spec(f2>0)),'b',f2(f2>0),imag(spec(f2>0)),'b:')
legend('Re(y)','Im(y)','Re(y after shift)','Im(y after shift)')
title('Fourier Transforms')
subplot(2,1,2)
plot(f2(f2>=0),mean(P1Y2,2),'r',f2(f2>=0),P1Y,'b')
legend('y','y after shift')
title('one sided power spectrum')


end

