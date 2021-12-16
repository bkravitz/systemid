% SIMPLE_CONTROL_RUN_COMPARE compares the average for each variable of the
% perturbed and control runs for the entire modeled time
%
% written by: Bethany Sutherland
% last edited: 6/20/2019
%

clear all; close all; clc
addpath ../helper_functions Figures/simple_control_run_compare


% Runs to include:
% runname = {'sysid_ocn_n005LATp005_190LON240_rndp_002_0K_000_0dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_sqwp_002_0K_366_304dy'};
runname = {'sysid_ocn_n005LATp005_190LON240_rnlp_002_0K_007_0dy'};

% runname = {'sysid_ocnTest31'};
% runname = {'sysid_ocnTest32_addinrAttr'};
% runname = {'sysid_ocnTest31_null'};

% runname = {'sysid_ocn_n015LATp015_120LON170_snwv_002_0K_007_0dy'}; 
% runname = {'sysid_ocn_n015LATp015_120LON170_snwv_002_0K_045_0dy'};
% runname = {'sysid_ocn_n015LATp015_120LON170_snwv_002_0K_090_0dy'};
% runname = {'sysid_ocn_n015LATp015_140LON220_rand_002_0K_000_0dy'};
% runname = {'sysid_ocn_n015LATp015_140LON220_snwv_002_0K_045_0dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_0dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_ady'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_bdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_cdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwv_002_0K_091_2dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_rndp_002_0K_000_0dy'}; 
% 
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_0dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_ady'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_bdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_cdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwv_002_0K_091_2dy'};

% runname = {'sysid_ocn_p065LATp080_340LON355_snwv_002_0K_045_0dy'};
% runname = {'sysid_ocn_p065LATp080_340LON355_snwv_002_0K_091_0dy'};
% runname = {'sysid_ocn_p065LATp080_340LON355_rnlp_002_0K_007_0dy'};
% runname = {'sysid_ocn_p065LATp080_165LON195_rnlp_002_0K_007_0dy'}; 
% 

control_runname = {'Control_20yr_b'};

% Variables to use:
Vars = {'TS','PRECT'};
% Vars = {'TS'};

% plot all on one figure or each individually
plot_together = true;
% plot_together = false;

% plot shortened figure titles
for_paper = true;
% short_titles = false;

coordVars = {'lat','lon','time'};

% Data Directory 
% Dir = 'C:\Users\suth827\OneDrive - PNNL\Documents\DATA\csmresults\'; % original on B.S. PnNl laptop
Dir = 'C:\Users\bsuther\Documents\PNNL Data\DATA\'; % on beths NCSU laptop

% length of run
yrs = 20;
styr = 1;

%% Other Book Keeping

Fs = 1; % sample/day 
L = 365*yrs;
stL = 365*styr - 364;
endL = stL + L - 1;

ttl_runname = convert_underscores(runname{:})

for vi = 1:numel(Vars)
    v = Vars{vi}
    
    % ____Load Control Run Data___________
    %_____________________________________
    fprintf('Loading Control Run Data \n')

    for vi2 = 1:numel(coordVars)
        v2 = coordVars{vi2};
        OGCtrlData.(v2) = nc_varget([Dir v '_' control_runname{:} '.nc'],v2);
    end
    
    OGCtrlData.(v) = nc_varget([Dir v '_' control_runname{:} '.nc'],v);
    if max(size(size(OGCtrlData.(v)))) == 4  % only look at lowest level data (for now)
        OGCtrlData.lev = nc_varget([Dir v '_' control_runname{:} '.nc'],'lev');
        OGCtrlData.(v) = squeeze(OGCtrlData.(v)(:,:,end,:));
    end
    
    CtrlData = OGCtrlData;
    CtrlData.(v) = OGCtrlData.(v)(:,:,stL:endL);
    CtrlData.time = OGCtrlData.time(stL:endL);
    
    
    % ____Load Model Data_________________
    %_____________________________________
    
    fprintf('Loading Model Data \n')
    
    for vi2 = 1:numel(coordVars)
        v2 = coordVars{vi2};
        OGMdlData.(v2) = nc_varget([Dir v '_' runname{:} '.nc'],v2);
    end
    
    OGMdlData.(v) = nc_varget([Dir v '_' runname{:} '.nc'],v);
    if max(size(size(OGMdlData.(v)))) == 4  % only look at lowest level data (for now)
        OGMdlData.(v) = squeeze(OGMdlData.(v)(:,:,end,:));
    end
    
    MdlData = OGMdlData;
    MdlData.(v) = OGMdlData.(v)(:,:,stL:endL);
    MdlData.time = OGMdlData.time(stL:endL);
    
    clearvars OGMdlData OGCtrlData
    
    % ____Load Perturbation Data_________________
    %_____________________________________
    
    fid = fopen([Dir runname{:} '_sysidout.bin'],'rb');  
    indata = fread(fid,'single','b');
    fclose(fid);
    Ptyrs = length(indata)/144/96; Ptyrs = (Ptyrs - 5)/365;
    PertDatatmp = reshape(indata,[144*96 Ptyrs*365+5]); % number of days + 5 header info columns
    % Columns of InData:
    % (1) lon# (2) lat# (3) real_lon (4) real_lat (5) box#(which region)
    % (6+) perturbation sequence
    platlon = 7825
    Pertseries = PertDatatmp(platlon,stL+5:endL+5); 
    
    
    %% Convert Units
    
    % convert precipitation from m/s to mm/day
    if ismember('PRECT',v)
        MdlData.(v) = MdlData.(v)*1000*86400;
        CtrlData.(v) = CtrlData.(v)*1000*86400;
    end
    
    %% plot a timeseries
    
    plat = 50 % same values as figure 2
    plon = 65 % same values as figure 2
%     plat = 49
%     plon = 82
    
    if plot_together
        % create figure
        figure('pos',[10 10 1200 800])
        ha = tight_subplot(2,2,[.07 0.05],[.1 .1],[.04 .04]);
    end
    
    tmsr = squeeze(MdlData.(v)(plon,plat,:));
    tmct = squeeze(CtrlData.(v)(plon,plat,:));
    if plot_together
        axes(ha(1))
    else
        figure('pos',[10 10 700 400])
    end
   

%     if for_paper
%         figure('pos',[10 10 600 800])
%         ha = tight_subplot(3,1,[.07 0.05],[.1 .1],[.04 .04]);
%         
%         
%         plot(1:L,tmsr,1:L,tmct)
%         yl = ylim; ydif = yl(2)-yl(1);
%         if ~strcmp(v,'PRECT')
%         ylim([yl(1)-0.2*ydif yl(2)])
%         end
%         if ismember('TS',v)
%             ylabel('K');
%         elseif ismember('PRECT',v)
%             ylabel('mm/day');
%         end
%         title('(a)')
%         legend('perturbed run','control run')
% 
%     else
%         yyaxis left
%         plot(1:L,tmsr,1:L,tmct)
%         yl = ylim; ydif = yl(2)-yl(1);
%         ylim([yl(1)-0.2*ydif yl(2)])
%         if ismember('TS',v)
%             ylabel('K');
%         elseif ismember('PRECT',v)
%             ylabel('mm/day');
%         end
%         yyaxis right
%         plot(1:L,Pertseries)
%         yr = ylim; ylim([yr(1) yr(1)+ydif])
%         ylabel('K')
%         legend('model','control','perturbation')
%         title({[ttl_runname],[v ' sample timeseries from perturbed region & pertubation']})
%     end

xticks(1:365:yrs*365)
xticklabels(1:yrs) % first of each month
xtickangle(45)
xlabel('year')
xlim([0 yrs*365])
    %% Sum temperatures over all time
    vnew = ['mean_' v];
    MdlData.(vnew) = nanmean(MdlData.(v),3);
    CtrlData.(vnew) = nanmean(CtrlData.(v),3);
    
    Difference = MdlData.(vnew) - CtrlData.(vnew);

    %% plot all data 
    if ismember('TS',v)
        levels = [-2.5:0.5:2.5 0.25 -0.25];
    elseif ismember('PRECT',v)
        levels = [-2.5:0.5:2.5 0.25 -0.25 0.1 -0.1];
    end
    if for_paper
        figure('pos',[10 10 600 800])
        ha = tight_subplot(3,1,[.07 0.05],[.1 .1],[.04 .04]);
        
        axes(ha(1))
    elseif plot_together
        axes(ha(2))
    else
        figure('pos',[10 10 700 400])
    end
    opts = [];
    [~,~,c] = discritizemap(MdlData.lat,MdlData.lon, Difference, levels, 0, false, opts, true);
    if plot_together
        if for_paper
            title('(a)')
            
        else
            title({[' average ' v],[' Perturbed - Control run ']})
        end
    else
        title([v ' Perturbed - Control run ' ttl_runname]);
    end
    if ismember('TS',v)
        c.Label.String = 'K';
    elseif ismember('PRECT',v)
        c.Label.String = 'mm/day';
    end

%% Split into DJF JJA
DJFi1 = [true(59,1); false(275,1); true(31,1)]; DJFi = DJFi1;
JJAi1 = [false(152,1); true(92,1); false(121,1)]; JJAi = JJAi1;
JANi1 = [true(31,1); false(334,1)]; JANi = JANi1;
while length(DJFi) < L
    DJFi = [DJFi; DJFi1];
    JJAi = [JJAi; JJAi1];
    JANi = [JANi; JANi1];
end


        
%% ______________Plot__DJF________________
vnew = ['DJF_mean_' v];
MdlData.(vnew) = nanmean(MdlData.(v)(:,:,DJFi),3);
CtrlData.(vnew) = nanmean(CtrlData.(v)(:,:,DJFi),3);

Difference = MdlData.(vnew) - CtrlData.(vnew);
    figure('pos',[10 10 700 400])
    opts{1} = [];
    if for_paper
            axes(ha(2))
    elseif plot_together
        axes(ha(3))
        
    else
        figure('pos',[10 10 700 400])
    end
    [~,~,c] = discritizemap(MdlData.lat,MdlData.lon, Difference, levels, 0, false, [], true);
    text(2,-0.85,'(DJF)','FontSize',14)
    if plot_together
        if for_paper
            title('(b)')
        else
            title({['DJF average ' v ],' Perturbed run - Control run '});
        end
    else
        title({['DJF average ' v ],' Perturbed run - Control run ', ttl_runname});
    end
    if ismember('TS',v)
        c.Label.String = 'K'
    elseif ismember('PRECT',v)
        c.Label.String = 'mm/day'
    end

    if ~plot_together
        % save figure
        figname = ['Figures/simple_control_run_compare/' ...
            runname{:} '_avgdeviationfromcontrol_'  v '_DJF'...
            ];
        print(gcf,'-dpng','-r300',figname);
        
    end
    
%% _______________Plot_JJA________________
    
vnew = ['JJA_mean_' v];
MdlData.(vnew) = nanmean(MdlData.(v)(:,:,JJAi),3);
CtrlData.(vnew) = nanmean(CtrlData.(v)(:,:,JJAi),3);

Difference = MdlData.(vnew) - CtrlData.(vnew);
if for_paper
            axes(ha(3))
elseif plot_together
    axes(ha(4))
else
    figure('pos',[10 10 700 400])
end
opts{1} = [];
[~,~,c] = discritizemap(MdlData.lat,MdlData.lon, Difference, levels, 0, false, [], true);
    text(2,-0.85,'(JJA)','FontSize',14)

if plot_together
    if for_paper
        title('(c)')
    else
        title({['JJA average ' v ],' Perturbed run - Control run '});
    end
else
    title({['JJA average ' v ],' Perturbed run - Control run ', ttl_runname});
end
if ismember('TS',v)
    c.Label.String = 'K'
elseif ismember('PRECT',v)
    c.Label.String = 'mm/day'
end

if ~plot_together
    % save figure
    figname = ['Figures/simple_control_run_compare/' ...
        runname{:} '_avgdeviationfromcontrol_'  v '_JJA'...
        ];
    print(gcf,'-dpng','-r300',figname);
else
    plot_together
    % save figure
    figname = ['Figures/simple_control_run_compare/' ...
        runname{:} '_avgdeviationfromcontrol_'  v ...
        ];
    if for_paper
        figname = [figname '_for_paper'];
    end
    print(gcf,'-dpng','-r300',figname);
end


%% _______________Plot_January only________________

levels = [-2.5:0.5:2.5];
vnew = ['JAN_mean_' v];
MdlData.(vnew) = nanmean(MdlData.(v)(:,:,JANi),3);
CtrlData.(vnew) = nanmean(CtrlData.(v)(:,:,JANi),3);

Difference = MdlData.(vnew) - CtrlData.(vnew);
    figure('pos',[10 10 700 400])
    opts{1} = [];
    [~,~,c] = discritizemap(MdlData.lat,MdlData.lon, Difference, levels, 0, false, [], true);
    title({['January average ' v ],' Perturbed run - Control run ', ttl_runname});
    if ismember('TS',v)
        c.Label.String = 'K'
    elseif ismember('PRECT',v)
        c.Label.String = 'mm/day'
    end
    
    % save figure
        figname = ['Figures/simple_control_run_compare/' ...
            runname{:} '_avgdeviationfromcontrol_'  v '_Jan'...
            ];
        print(gcf,'-dpng','-r300',figname);
        
end


