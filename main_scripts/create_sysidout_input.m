% CREATE_SYSIDOUT_INPUT creates a binary input matrix for a specified 
% perturbation for a region in CESM between specified lon and lat values
%
% uses helper_functions/sysidout_generic_22regions.bin as a template

% clearvars; close all; clc
addpath helper_functions

% Use regions from ben's original 22 regions
use_bens_regions = false;
reg = [];       % choose regions to perturb

% OR input region(s) to perturb
lons = [
%     140 220;
%         120 170
        190 240
%         340 355
%         165 195
        ];   % % nino 3.4 lats are [190 240] East (170-120W)
lats = [
%         -15 15;
        -5 5;
%         65 80;
%         -70 -55
        ];    % nino 3.4 lats are [5 -5]

 


% input perturbation type
% perttype = 1; % Sinewave
% perttype = 1.5; % Sinewave ALL POSITIVE
% perttype = 1.6; % multiple sinewaves added together
% perttype = 2; % Random
% perttype = 2.5; % Random ALL POSITIVE
% perttype = 2.6; % Random with low pass filter
% perttype = 3; % Squarewave
perttype = 4; % one sided squareish pulse waves
% perttype = 4.5; % two sided squareish pulse waves

% Input amplitude
amp = 2;

% Input length of run
yrs = 20;
%% Sinewave Inputs

% Input period (days)
% 32 days = 10*2qrt(10)
% T = 91;
T = 330;
ph = 0;         % 91_0
% ph = pi/2;    % 91_a
% ph = pi;      % 91_b
% ph = 3*pi/2;  % 91_c

%% Load Template

fid = fopen('sysidout_generic_22regions.bin','rb');
indata = fread(fid,'single','b');
fclose(fid);
outputmatrix = reshape(indata,[144*96 5]); % number of days + 5 header info columns,
% Columns of InData:
% (1) lon# (2) lat# (3) real_lon (4) real_lat (5) box#(which region)
% (6+) perturbation sequence

% adjust size of matrix for the length of the run
outputmatrix = [outputmatrix zeros(length(outputmatrix), yrs*365)];

t = 0:yrs*365-1;


%% Define Regions

if use_bens_regions == false
    % redefine regions from lat/lon values
    outputmatrix(:,5) = 0;
    reg = 1:size(lons,1);
    for nr = 1:size(lons,1)
        for k = 1:length(outputmatrix)
            if outputmatrix(k,3) > lons(nr,1) && outputmatrix(k,3) < lons(nr,2)
                if outputmatrix(k,4) > lats(nr,1) && outputmatrix(k,4) < lats(nr,2)
                    outputmatrix(k,5) = nr;
                end
            end
        end
    end
end


%% Determine Perturbation

if perttype == 1 % SineWave Perturbation
    pert = amp*sin(2*pi*t/T+ph);
elseif perttype == 1.5 % SineWave Perturbation ALL POSITIVE
    pert = amp*sin(2*pi*t/T+ph)+amp;
elseif perttype == 1.6 % 3 SineWave Perturbations ALL POSITIVE
    pert = amp*sin(2*pi*t/7) + amp*sin(2*pi*t/45) + amp*sin(2*pi*t/91);
elseif perttype == 2 % Random perturbation
    pert = 2*amp*rand(1,length(t))-amp;
elseif perttype == 2.5 % Random perturbation ALL POSITIVE
    pert = amp*rand(1,length(t));
elseif perttype == 2.6 % Random perturbation with low pass filter
    pert1 = rand(1,length(t)) - 0.5;
    [b,a] = butter(5,2/7,'low');
    pert = filtfilt(b,a,pert1);
    mx = max(abs(pert));
    pert = amp*pert/mx;
elseif perttype == 3 % square wave perturbation
elseif perttype == 4 || perttype == 4.5 % two sided squareish pulse waves
    pd = 180; % pulse duration
    ramp = 92; % time to ramp up to full perturbation
    gap = [60 365*2-pd-2*(ramp+1)-60]; % gap between pulses
    rampup = (0:ramp)'/ramp;
    rampdn = flipud(rampup);
    if perttype == 4
        pertcycle = [rampup; ones(pd,1); rampdn; zeros(gap(2),1)];
    pert = [zeros(gap(1),1); pertcycle];
    elseif perttype == 4.5
        pert = [rampup; ones(pd,1); rampdn; zeros(gap,1); -1*rampup; -1*ones(pd,1); -1*rampdn; zeros(gap,1); ];
    end
    
    while length(pert) < length(t)
        pert = [pert; pertcycle];
    end
    pert = amp*pert(1:length(t))+(amp/50)*rand(length(t),1);

end

% plot perturbation to visualize
figure
plot(1:10*T,pert(1:10*T))

%% add perturbation to signal

for rn = 1:length(reg)
    for k = 1:size(outputmatrix,1)
        if outputmatrix(k,5) == reg(rn)
            outputmatrix(k,6:end) = pert;
        end
    end
end


%% Save output

fid=fopen('sysidout.bin','wb');
fwrite(fid,outputmatrix,'single','b');
fclose(fid);



