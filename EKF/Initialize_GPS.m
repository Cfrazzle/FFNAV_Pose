%% InitializeGPS ======================================================
% Description: This script initializes the GPS constellation.
%
% Inputs:
%   Parameter_filename  - Spacecraft and formation parameters filename
%   PropOptions.        - Structure for orbit propagator settings
%   time_start          - Start time for the simulation
%
% Outputs:
%   Results.mat file
%
% Other Functions Called:
%
% Created by:  Cory Fraser - NOV 18, 2018
% Latest Edit: Cory Fraser - NOV 18, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Initialization

%fprintf('\n--------------------------------------------------------------\n')
fprintf('\n Loading GPS Measurement Parameters...')


%% GPS Measurement Specifications

% GPS Errors - Standard deviations for white gaussian noise
GPS_r_STD = [1.20 1.20 1.20]; %(m)  
GPS_v_STD = [0.03 0.03 0.03]; %(m/s)

% Measurement Sampling Time
r_target_sensor_dt = 1;
r_chaser_sensor_dt = 1;
v_target_sensor_dt = 1;
v_chaser_sensor_dt = 1;

% Sensor Flags (1 = Sensor On)
r_target_sensor_flag = 1;
r_chaser_sensor_flag = 1;
v_target_sensor_flag = 1;
v_chaser_sensor_flag = 1;

% Measurement Flags (1 = Measurement is good, process it)
r_rel_meas_flag = 1;
v_rel_meas_flag = 1;
theta_meas_flag = 1;

% Noise Seeds
r_target_NoiseSeed = [1 2 3];
r_chaser_NoiseSeed = [4 5 6];
v_target_NoiseSeed = [7 8 9];
v_chaser_NoiseSeed = [10 11 12];

fprintf(' Complete \n')
% =========================================================================