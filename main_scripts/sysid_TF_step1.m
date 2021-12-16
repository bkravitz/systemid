% SYSID_TF_STEP1 uses only the transfer function amplitudes to analyse
% the output of the perturbed sysid runs
%
% written by: Bethany Sutherland
% last edited: 5/29/2020
%
% 
% The basic steps performed by this code are:
% 
% 1) Loads control run, model run, and perturbation data
% 2) if use_anomaly = true removes anomalies from data (takes out average for that 
%    day of the year for all years of the run)
% 2) Compute the transfer function between the input and both the sysid run
%    and the control run for each grid point 
% 3) Saves data into a matrix for future use by sysid_TF_step2.m

clear all; close all; clc
addpath ../../helper_functions Figures ../../Gest_From_MnT Data_Output


% Runs to include:
% runname = {'sysid_ocnTest31'};
% runname = {'sysid_ocnTest32_addinrAttr'};
% runname = {'sysid_ocnTest31_null'};
% runname = {'sysid_ocn_n015LATp015_120LON170_snwv_002_0K_007_0dy'}; % done before adding lags
% runname = {'sysid_ocn_n015LATp015_120LON170_snwv_002_0K_045_0dy'};
% runname = {'sysid_ocn_n015LATp015_120LON170_snwv_002_0K_090_0dy'};
% runname = {'sysid_ocn_n015LATp015_140LON220_rand_002_0K_000_0dy'};
% runname = {'sysid_ocn_n015LATp015_140LON220_snwv_002_0K_045_0dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_0dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_ady'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_bdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_cdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwv_002_0K_091_2dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_rndp_002_0K_000_0dy'}; % done only TS, PRECT,CLDOT
% 
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_0dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_ady'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_bdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwp_002_0K_091_cdy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_snwv_002_0K_091_2dy'};
% runname = {'sysid_ocn_n005LATp005_190LON240_rndp_002_0K_000_0dy'};

% runname = {'sysid_ocn_p065LATp080_340LON355_snwv_002_0K_045_0dy'};
% runname = {'sysid_ocn_p065LATp080_340LON355_snwv_002_0K_091_0dy'};
% runname = {'sysid_ocn_p065LATp080_340LON355_rnlp_002_0K_007_0dy'};
% runname = {'sysid_ocn_p065LATp080_165LON195_rnlp_002_0K_007_0dy'}; % done only TS, prect, cldtot
runname = {'sysid_ocn_n005LATp005_190LON240_rnlp_002_0K_007_0dy'};


% number of fourier transform points
nfft = [512,floor(10^2.8),floor(10^3),floor(10^3.2)];
control_runname = {'Control_20yr_b'};

% Variables to use:
Vars = {'TS','PRECT','CLDTOT'};
% Vars = {'TS'};

coordVars = {'lat','lon','time'};

% Data Directory 
% Dir = 'C:\Users\suth827\OneDrive - PNNL\Documents\DATA\csmresults\';
Dir = 'C:\Users\bsuther\Documents\PNNL Data\DATA\';

% length of run
yrs = 20;
styr = 1;

% remove anomalies
use_anomaly = true;

makeTFplots = false;

%% Other Book Keeping

Fs = 1; % sample/day 
L = 365*yrs;
stL = 365*styr - 364;
endL = stL + L - 1;

if use_anomaly
    dataoutName = ['Data_Output/' ...
        runname{:} '_anomTFresults.mat']
else
    dataoutName = ['Data_Output/' ...
        runname{:} '_TFresults.mat']
end

%% Do analysis

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
    
    
    %____Load Perturbation Data___________
    %_____________________________________
    
    if ~exist('PertData','var')
        fprintf('Loading Perturbation Data \n')
        fid = fopen([Dir runname{:} '_sysidout.bin'],'rb');  % MAKE GENERAL LATER!!!!
        indata = fread(fid,'single','b');
        fclose(fid);
        Ptyrs = length(indata)/144/96; Ptyrs = (Ptyrs - 5)/365;
        PertDatatmp = reshape(indata,[144*96 Ptyrs*365+5]); % number of days + 5 header info columns
        % Columns of InData:
        % (1) lon# (2) lat# (3) real_lon (4) real_lat (5) box#(which region)
        % (6+) perturbation sequence
        PertData(:,1:5) = PertDatatmp(:,1:5);
        PertData(:,6:L+5) = PertDatatmp(:,stL+5:endL+5);
        
        clear indata
        
        
        pertregs = PertData(:,5);
        numperts = max(pertregs);
        
        for k = 1:numperts
            [a,~] = find(pertregs == k);
            rw = min(a);
            
            % determine input
            x(1,:) = PertData(rw,6:end);  % HARDCODED!
            
        end
        
        % (1) lon# (2) lat# (3) real_lon (4) real_lat (5) box#(which region)
        PertMap = zeros(size(MdlData.(Vars{1})));
        PertMap = squeeze(PertMap(:,:,1));
        for k = 1:size(PertData,1)
            PertMap(PertData(k,1),PertData(k,2)) = PertData(k,5);            
        end
        PertMap(PertMap ~= 0) = PertMap(PertMap ~= 0) + 1000; % to differentiate regions (ensures pertmap values are higher than the number of regions)
    end
    
    % Convert Units
    
    % convert precipitation from m/s to mm/day
    if ismember('PRECT',v)
        MdlData.(v) = MdlData.(v)*1000*86400;
        CtrlData.(v) = CtrlData.(v)*1000*86400;
    end
    
    % Calculate and remove anomalies (if use_anomaly = true)
    %_____________________________________
    
    if use_anomaly 
        
        dayNum = 1:365;
        dayNum = repmat(dayNum',yrs,1);
        for n = 1:365
            includedays = false(L,1);
            includedays(dayNum == n) = true;
            dayClim(:,:,n) = mean(CtrlData.(v)(:,:,includedays),3);
        end
        dayClim = repmat(dayClim,1,1,yrs);
        
        MdlData.(v) = MdlData.(v) - dayClim;
        CtrlData.(v) = CtrlData.(v) - dayClim; 
    end
    
    
    
    % Calculate transfer functions
    %_____________________________________
    totr = length(MdlData.lon)*length(MdlData.lat);
    r = 0;
    fprintf(['calculating TF for ' v '\n'])
%     for latindx = 1:length(MdlData.lat)
%         for lonindx = 1:length(MdlData.lon)
    for latindx = 49:50
        for lonindx = 80:81
            fprintf(['lat = ' num2str(latindx) '\n'])
            fprintf(['lon = ' num2str(lonindx) '\n'])
            r = r+1;
            
              % indicate progress
            if r == floor(totr/100) 
                fprintf('1%% done \n')
            elseif r == floor(totr/10)
                fprintf('10%% done \n')
            elseif r == floor(totr/4)
                fprintf('25%% done \n')
            elseif r == floor(2*totr/4)
                fprintf('50%% done \n')
            elseif r == floor(3*totr/4)
                fprintf('75%% done \n')
            end
%             
%             if latindx < 45 && latindx >38
%                 if lonindx > 46 && lonindx < 60
                    makeTFplots = true;
%                 end
%             end
            [H, H_c, C, F, lF, xF, he, he_c] = tranFval2(x,squeeze(MdlData.(v)(lonindx,latindx,:))',squeeze(CtrlData.(v)(lonindx,latindx,:))',...
                nfft,Fs,[],-1,1000*vi+r, makeTFplots );
            [H, H_c, C, F, lF, xF, he, he_c] = tranFval_forpaper(x,squeeze(MdlData.(v)(lonindx,latindx,:))',squeeze(CtrlData.(v)(lonindx,latindx,:))',...
                nfft,Fs,[],-1,1000*vi+r, makeTFplots );
%             if latindx < 45 && latindx >38
%                 if lonindx > 46 && lonindx < 60
%                     title(['lat = ' num2str(MdlData.lat(latindx)) ' lon = ' num2str(MdlData.lon(lonindx)) 'pert?=' num2str(PertMap(lonindx,latindx)) ])
%                     makeTFplots = false;
%                 end
%             end
            % determine optimal lags 
            [lagout,Corr] = lagcorrect(x,squeeze(MdlData.(v)(lonindx,latindx,:))',3*365, false);
            F = F(F>0);  % remove zero from F list
            
            Finc = ismember(F,lF);
            
            VarData.(v)(lonindx,latindx).Run.H = H(2:end); % remove F = 0 data
            VarData.(v)(lonindx,latindx).Run.TF = he.YData;
            VarData.(v)(lonindx,latindx).Run.Err = he.YPositiveDelta;
            VarData.(v)(lonindx,latindx).Control.H = H_c(2:end);
            VarData.(v)(lonindx,latindx).Control.TF = he_c.YData;
            VarData.(v)(lonindx,latindx).Control.Err = he_c.YPositiveDelta;
            VarData.(v)(lonindx,latindx).Finc = Finc;
            VarData.(v)(lonindx,latindx).lF = lF;
            VarData.(v)(lonindx,latindx).lag = lagout;
            VarData.(v)(lonindx,latindx).correlation = Corr;
            if lonindx == 1 && latindx == 1
                allSigFlist = false(size(F));
                allSigFmagt = zeros(size(F));
            end
            
            allSigFlist(Finc) = true;
            allSigFmagt = allSigFmagt + Finc;
        end
    end
    VarData.(v)(1,1).allSigF = allSigFlist;
    VarData.(v)(1,1).allSigFmag = allSigFmagt;
    
    %____________Save  Data_______________
    %_____________________________________
%     
    % check if data matrix already exists
    previousMatrix = exist(dataoutName);
    % dataout = matfile(dataoutName,'Writable',true);
    if ~previousMatrix % if not create add general info to new matrix
        fprintf('Creating a New Matrix \n')
        RunInfo.runname = runname{:};
        RunInfo.F = F;
        RunInfo.PertMap = PertMap;
        RunInfo.perturbation = x;
        RunInfo.control_runname = control_runname{:};
        RunInfo.lon = MdlData.lon;
        RunInfo.lat = MdlData.lat;
        save(dataoutName, 'RunInfo','-v7.3')
    end
    
    
    if ismember('TS',v)
        TS = VarData.(v);
        save(dataoutName, 'TS','-append')
    elseif ismember('PRECT',v)
        PRECT = VarData.(v);
        save(dataoutName, 'PRECT','-append')
    elseif ismember('CLDTOT',v)
        CLDTOT = VarData.(v);
        save(dataoutName, 'CLDTOT','-append')
    elseif ismember('PS',v)
        PS = VarData.(v);
        save(dataoutName, 'PS','-append')
    else
        fprintf('VARIABLE DATA NOT SAVED! \n ')
        return
    end
    
    clearvars VarData MdlData CtrlData
end
