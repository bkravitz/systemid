% SYSID_TF_STEP2 uses transfer function amplitudes calculated in sysid_TF_step1.m
% to create plots
%
% written by: Bethany Sutherland
% last edited: 5/29/2020
%
% 
% The basic steps performed by this code are:
% 
% 1) load transfer function matrix created in sysid_TF_step1.m
% 2) plot results

clear all; close all; clc

addpath Data_Output Figures ../../helper_functions

% runname = {'sysid_ocn_n005LATp005_190LON240_rndp_002_0K_000_0dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_rndp_002_0K_000_0dy_anom'};
% runname = {'sysid_ocn_n005LATp005_190LON240_rnlp_002_0K_007_0dy'};
runname = {'sysid_ocn_n005LATp005_190LON240_rnlp_002_0K_007_0dy_anom'};
% runname = {'sysid_ocn_n005LATp005_190LON240_sqwp_002_0K_366_304dy'};

% runname = {'sysid_ocn_n015LATp015_120LON170_snwv_002_0K_007_0dy'}; % step 1 done: TS
% runname = {'sysid_ocn_n015LATp015_120LON170_snwv_002_0K_045_0dy'};
% runname = {'sysid_ocn_n015LATp015_120LON170_snwv_002_0K_090_0dy'};
% runname = {'sysid_ocn_n015LATp015_140LON220_rand_002_0K_000_0dy'};
% runname = {'sysid_ocn_n015LATp015_140LON220_snwv_002_0K_045_0dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_0dy'}; 
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_ady'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_bdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_cdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwv_002_0K_091_2dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_rndp_002_0K_000_0dy'}; % step 1 done: TS, PRECT, CLDTOT

% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_0dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_ady'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_bdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_cdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwv_002_0K_091_2dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_rndp_002_0K_000_0dy'};

% runname = {'sysid_ocn_p065LATp080_340LON355_snwv_002_0K_045_0dy'};
% runname = {'sysid_ocn_p065LATp080_340LON355_snwv_002_0K_091_0dy'};
% runname = {'sysid_ocn_p065LATp080_340LON355_rnlp_002_0K_007_0dy'};
% runname = {'sysid_ocn_p065LATp080_165LON195_rnlp_002_0K_007_0dy'}; % step 1 done TS, PRECT, CLDTOT

% runname = {'sysid_ocn_n005LATp005_190LON240_snwv_002_0K_091_2dyControlTest'}; % step 1 done TS, PRECT, CLDTOT


DataMatrix = [runname{:} '_TFresults'];

Vars = {'TS','PRECT'};
 ptT = [7];


% get only necessary info from
M = matfile(DataMatrix);
Run = M.RunInfo;

%% do analysis 

for vi = 1:numel(Vars)
    v = Vars{vi};
    Var = M.(v);
    Zz = []; Yy = []; Xx = []; CcDiff = []; CcRatio = []; 
    CcLag_pk = []; CcLag_1 = []; CcLag_15 = []; CcLag_2 = [];
    
    TrackAllSigF = Var(1,1).allSigF;
    
    for latindx = 1:length(Run.lat)
        for lonindx = 1:length(Run.lon)
            Mask(lonindx,latindx,:) = Var(lonindx,latindx).Finc;
            TFModel(lonindx,latindx,:) = Var(lonindx,latindx).Run.TF;
            TFControl(lonindx,latindx,:) = Var(lonindx,latindx).Control.TF;
            ModelErr(lonindx,latindx,:) = Var(lonindx,latindx).Run.Err;
            ControlErr(lonindx,latindx,:) = Var(lonindx,latindx).Control.Err;
            if sum(Var(lonindx,latindx).Finc)>0
                Zz = [Zz; Run.F(Var(lonindx,latindx).Finc)];
                numF = size(Run.F(Var(lonindx,latindx).Finc),1);
                Xx = [Xx; Run.lon(lonindx)*ones(numF,1)];
                Yy = [Yy; Run.lat(latindx)*ones(numF,1)];
                CcDiff = [CcDiff; Var(lonindx,latindx).Run.TF(Var(lonindx,latindx).Finc)' - Var(lonindx,latindx).Control.TF(Var(lonindx,latindx).Finc)'];
                CcRatio = [CcRatio; Var(lonindx,latindx).Run.TF(Var(lonindx,latindx).Finc)'./Var(lonindx,latindx).Control.TF(Var(lonindx,latindx).Finc)'];
                CcLag_pk = [CcLag_pk; Var(lonindx,latindx).lag(1)*ones(numF,1)]; 
                CcLag_1 = [CcLag_1; Var(lonindx,latindx).lag(2)*ones(numF,1)]; 
                CcLag_15 = [CcLag_15; Var(lonindx,latindx).lag(3)*ones(numF,1)]; 
                CcLag_2 = [CcLag_2; Var(lonindx,latindx).lag(4)*ones(numF,1)]; 
%             if sum(Var(lonindx,latindx).Finc)>5
%                 fprintf('stophere')
%             end
            end
        end
    end
    
    MskTFModel = TFModel; MskTFModel(~Mask) = NaN;
    MskTFControl = TFControl; MskTFControl(~Mask) = NaN;
    
    
    for f = 1:length(Run.F)
        totalSigFarea(f) = sum(sum(Mask(:,:,f)));
        tsfa_Diff(f) = sum(nansum(MskTFModel(:,:,f)-MskTFControl(:,:,f)));
        tsfa_Model(f) = sum(nansum(MskTFModel(:,:,f)));
        tsfa_Control(f) = sum(nansum(MskTFControl(:,:,f)));
    end
    

    
    %% plot results
    
    % visualize amount of significant signals for each frequency
    figure(1000)
    plot(Run.F(2:end),totalSigFarea(2:end)./sum(totalSigFarea(2:end)),Run.F(2:end),tsfa_Diff(2:end),Run.F(2:end),tsfa_Model(2:end)./sum(tsfa_Model(2:end)),Run.F(2:end),tsfa_Control(2:end)./sum(tsfa_Control(2:end)))
    yl = ylim;
    legend('num of sig F','Differnce','Model num of sig F weighted by ratio','Control num of sig F weighted by ratio')
    
    hold on
    plot([1/ptT 1/ptT],yl)
    
    Findx = find(TrackAllSigF);
    
        % visualize amount of significant signals for each frequency
    figure(1001)
    semilogy(Run.F,totalSigFarea./sum(totalSigFarea),Run.F,tsfa_Model./sum(tsfa_Model),Run.F,tsfa_Control./sum(tsfa_Control))
    yl = ylim;
    legend('num of sig F','Model num of sig F weighted by ratio','Control num of sig F weighted by ratio')
    
    hold on
    plot([1/ptT 1/ptT],yl)
    
    %% Plot lags 
    
    figure()
    xlag = 1:length(CcLag_pk);
    semilogy(xlag,CcLag_pk,xlag,CcLag_1,xlag,CcLag_15, xlag, CcLag_2)
    legend('1st local max','1st 1\sigma','1st 1.5\sigma','1st 2\sigma')
    
    %% Create 3d plot
    
    % 0 doesnt get plotted so ensure 0 lags get plotted
    CcLag_pk(CcLag_pk == 0) = 0.001;
    CcLag_pk(CcLag_1 == 0) = 0.001;
    CcLag_pk(CcLag_15 == 0) = 0.001;
    CcLag_pk(CcLag_2 == 0) = 0.001;

        
    if ~exist('coast','var')
        coast = load('coast.mat');
        
        coast.long(coast.long<0) = coast.long(coast.long<0)+ 360;
        for ci = 1:length(coast.long)
            if abs(coast.long(ci+1) - coast.long(ci))>355
                coast.long(ci+2:end+1) = coast.long(ci+1:end);
                coast.long(ci+1) = NaN;
                coast.lat(ci+2:end+1) = coast.lat(ci+1:end);
                coast.lat(ci+1) = NaN;
            end
        end
    end
    
    [box] = getbox4pcolormap(runname{:},Run.lon,Run.lat);
    
    
    sigma = std(CcDiff); sigma2 = std(CcRatio);
    avgC = mean(CcDiff); avgC2 = mean(CcRatio);
    cLims = [avgC - 1*sigma, avgC + 1*sigma];
%     cLims2 = [avgC2 - 1*sigma2, avgC2 + 1*sigma2];
    cLims2 = [0 20]
    
% _____plot difference______________________________________________ 
    ttl = [v ' Model - Control: transfer functions ' convert_underscores(runname{:}) ];
    lbl = 'model - control';
    figname = ['Figures/' runname{:} '/' ...
        runname{:} '_' v '_TF_difference'...
        ];
    map4axis(Xx,Yy,Zz,CcDiff,ttl,lbl,figname,box,cLims)
    
 
%% _____plot ratios______________________________________________ 
    ttl = [v ' Model/Control: transfer functions ' convert_underscores(runname{:}) ];
    lbl = 'model/control';
    figname = ['Figures/' runname{:} '/' ...
        runname{:} '_' v '_TF_ratio'...
        ];
    climR = [0 20];
    map4axis(Xx,Yy,Zz,CcRatio,ttl,lbl,figname,box,climR)
 
%% _____plot ratios for paper______________________________________ 

%     figname = ['Figures/' runname{:} '/' ...
%         runname{:} '_' v '_TF_ratio_for_paper'...
%         ];
%     lbl = 'model/control';
%     climR = [0.1 20];
%     map1axis(Xx,Yy,Zz,CcRatio,lbl,figname,box,climR)
    
    figname = ['Figures/' runname{:} '/' ...
        runname{:} '_' v '_TF_ratio_for_paper'...
        ];
    lbl = 'model/control';
    climR = [1/365/10 0.5];
    climR2 = [0.1 20];
    CcF = Zz;
    % add randomness to plotting
%     Xxrnd = Xx;
%     Yyrnd = Yy;
    Xxrnd = Xx + 3*rand(size(Xx))-1;
    Yyrnd = Yy + 3*rand(size(Yy))-1;
    map1axis_2colorscales(Xxrnd,Yyrnd,Zz,CcF,CcRatio,lbl,figname,box,climR,climR2)    


%% _____plot lags______________________________________________ 

    % 1st local maximum in correlation
    ttl = [v ' lag (1st local max):  transfer functions ' convert_underscores(runname{:}) ];
    lbl = 'lag (days)';
    figname = ['Figures/' runname{:} '/' ...
        runname{:} '_' v '_TF_lag_1pk'...
        ];
    levels = [0 7 14 21 28 60 90 120 150 210 270 365 730 1095];
    map4axis(Xx,Yy,Zz,CcLag_pk,ttl,lbl,figname,box,[],[],levels)
     
    
    % 1st local maximum above 1 standard deviation from average correlation
    ttl = [v ' lag (1st above 1\sigma):  transfer functions ' convert_underscores(runname{:}) ];
    lbl = 'lag (days)';
    figname = ['Figures/' runname{:} '/' ...
        runname{:} '_' v '_TF_lag_1sig'...
        ];
    map4axis(Xx,Yy,Zz,CcLag_1,ttl,lbl,figname,box,[],[],levels)
    
    % 1st local maximum above 1.5 standard deviation from average correlation
    ttl = [v ' lag (1st above 1.5\sigma):  transfer functions ' convert_underscores(runname{:}) ];
    lbl = 'lag (days)';
    figname = ['Figures/' runname{:} '/' ...
        runname{:} '_' v '_TF_lag_15sig'...
        ];
    map4axis(Xx,Yy,Zz,CcLag_15,ttl,lbl,figname,box,[],[],levels)
    
    % 1st local maximum above 2 standard deviation from average correlation
    ttl = [v ' lag (1st above 2\sigma):  transfer functions ' convert_underscores(runname{:}) ];
    lbl = 'lag (days)';
    figname = ['Figures/' runname{:} '/' ...
        runname{:} '_' v '_TF_lag_2sig'...
        ];
    map4axis(Xx,Yy,Zz,CcLag_2,ttl,lbl,figname,box,[],[],levels)
 

% _____show frequency pattern______________________________________________ 
    ttl = [v ' transfer functions ' convert_underscores(runname{:}) ];
    lbl = 'frequency';
    figname = ['Figures/' runname{:} '/' ...
        runname{:} '_' v '_TF_frequencies'...
        ];
    map4axis(Xx,Yy,Zz,Zz,ttl,lbl,figname,box)
  
  
end
    

