Workflow for System Identification - Transfer Function Method - Codes

author: Bethany Sutherland
last edited: 9/25/2019

These are the codes used to create the figures in: 
System identification techniques for detection of teleconnections within climate models
(Bethany Sutherland, Ben Kravitz, Philip J. Rasch, and Hailong Wang )

simple_control_run_compare.m -- used to create Figures 3 and 4
% SIMPLE_CONTROL_RUN_COMPARE.m compares the average for each variable of the
% perturbed and control runs for the entire modeled time


sysid_TF_step1.m -- used to create Figure 2 
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


sysid_TF_step2.m -- used to create Figures 5 and 6
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


% CREATE_SYSIDOUT_INPUT creates a binary input matrix for a specified 
% perturbation for a region in CESM between specified lon and lat values
