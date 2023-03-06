%% FFNAV_FAEKF_RunSim =====================================================
% FFNAV Fuzzy Adaptive Extended Kalman Filter Implementation ==============
% Description: This script executes a fuzzy adaptive extended kalman filter
% for relative navigation of a spacecraft formation. The exact nonlinear
% equations of relative motion are utilized in the dynamics model, and GPS
% measurements are simulated using an orbit propogator with various
% perturbations and noise sources.
%
% Inputs:
%   Parameter_filename  - Spacecraft and formation parameters filename
%   EKF_flag            - Select filter (EKF, MLE-EKF, or FAEKF)
%   Q_adapt             - Choose Q-adaptations (ON/OFF)
%   R_adapt             - Choose R-adaptations (ON/OFF)
%   Smooth_flag         - Choose residual fixed-window smoothing (ON/OFF)
%   N_window            - Fixed-window size for smoothing
%   PropOptions.        - Structure for orbit propagator settings
%   time_start          - Start time for the simulation
%   time_step           - Time step of the simulation
%   orbit_num           - Number of orbits to simulate
%   time_end            - End time of the simulation
%
%   save_flag           - Choose to save output to a .mat file
%   post_flag           - Choose to post-process the data
%   n_start             - Starting point for post-processing
%   n_end               - Ending point for post-processing
%   plot_flag           - Choose to create output plots
%   print_flag          - Choose to print plots to .eps files
%   profiler_flag       - Choose to show performance profiler
%
% Outputs:
%   Results.mat file
%   Figures & animations (optional)
%   RMS Error Analysis (optional)
%
% Other Functions Called:
%   FFNAV_EKF_MakeParams.m      - Create FAEKF parameters (if none found)
%   OrbitalMechanicsConstants   - Initialize astrodynamics constants
%
% Created by:  Cory Fraser - JAN 01, 2023
% Latest Edit: Cory Fraser - FEB 23, 2023
% Copyright(c) 2023 by Cory Fraser
% =========================================================================
%% Initialization
clear
close all
profiler_flag = false;

fprintf('\n--------------------------------------------------------------\n')
fprintf(' \t\t FFNAV EXTENDED KALMAN FILTER SIMULATION \t\t             \n')
fprintf('--------------------------------------------------------------\n')

% Define path to parameters/output data
paths.pwd           = fullfile(fileparts('C:\SoftwareRepos\FFNAV_Pose\'));
paths.EKF           = fullfile([paths.pwd, '\EKF\']);
paths.orbit_prop    = fullfile([paths.pwd, '\Orbit_Propagator\']);
paths.utils         = fullfile([paths.pwd, '\utilities\']);
paths.output        = fullfile(fileparts('C:\SoftwareRepos\FFNAV_Pose_Output\'));
input_data.paths    = paths;

% Add paths
thesePaths = fieldnames(paths);
for iPath = 1:length(thesePaths)
    addpath(paths.(thesePaths{iPath}))
end

% Run the simulation
tic;
simulate_EKF(input_data);
t_run = toc;

fprintf('EKF simulation complete \n')
fprintf('Runtime = %f \n', t_run)
fprintf('..............\\(-_-)/........................\n')