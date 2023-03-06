function simulate_EKF(input_data)
%% FFNAV_FAEKF_RunSim =====================================================
% FFNAV Fuzzy Adaptive Extended Kalman Filter Implementation ==============
% Description: This script executes a fuzzy adaptive extended kalman filter
% for relative navigation of a spacecraft formation. The exact nonlinear
% equations of relative motion are utilized in the dynamics model, and GPS
% measurements are simulated using an orbit propogator with various
% perturbations and noise sources.
%
% Inputs:
%
% Outputs:
%
% Other Functions Called:
%
%
% Created by:  Cory Fraser - FEB 05, 2023
% Latest Edit: Cory Fraser - FEB 24, 2023
% Copyright(c) 2023 by Cory Fraser
% =========================================================================
%% Initialization

Day     = 29;   %Range: 1 - 31
Month   = 11;   %Range: 1-12
Year    = 2018; %Range 1901-2099

Hour    = 0;    %Range: 0-23
Minute  = 0;    %Range: 0-60
Second  = 0;    %Range: 0-60

%Calculate Julian Day Number
JulianDate0 = JulianDate(Day,Month,Year,Hour,Minute,Second);

% Load initial parameters for real-world simulator
Initialize_RWSW

% Load initial parameters for GPS constellation
Initialize_GPS

% Load initial parameters for on-board software
Initialize_RelNavEKF_OBSW

%% Select Simulation Options
time.start     = 0;                           % Start time of the simulation
RWSW_dt        = 1;                           % Time step for the simulation
time.step      = RWSW_dt;
orbit_num      = 2;                           % Number of orbits to simulate
time.stop      = ceil(orbit_num*COE_target.T);    % End time of the simulation
time.sim_time  = time.start:time.step:time.stop;

% Post Processing Options
save_flag     = true;     % Save data to a .mat file
post_flag     = true;     % Post-process the simulation data
    n_start   = 1; %floor(time_end*0.5); 
    n_end     = time.stop+1;  
    cov_flag  = 'on';    % Perform analysis of covariances
plot_flag     = 0;     % Create output plots
print_flag    = 'off';    % Print plots to .eps files
profiler_flag = false;    % Show profile simulink performance

%% Real-world Simulator

% Propagate target and chaser orbits
ode_options = odeset('relTol', 1e-6, 'AbsTol', 1e-6);
[~,RV_target] = ode45(@(t,y) accelerations_diff_eq(t,y,Constants,PropOptions_target), time.sim_time, [R0_target, V0_target],ode_options);
[~,RV_chaser] = ode45(@(t,y) accelerations_diff_eq(t,y,Constants,PropOptions_chaser), time.sim_time, [R0_chaser, V0_chaser],ode_options);

% Seperate position/velocity data for target/chaser
r_target_ECI = [RV_target(:,1), RV_target(:,2), RV_target(:,3)];
v_target_ECI = [RV_target(:,4), RV_target(:,5), RV_target(:,6)];
r_chaser_ECI = [RV_chaser(:,1), RV_chaser(:,2), RV_chaser(:,3)];
v_chaser_ECI = [RV_chaser(:,4), RV_chaser(:,5), RV_chaser(:,6)];

%% Add measurement noise to ECI position/velocities

% Initialize noise variables
GPS_r_target_noise = zeros(length(time.sim_time),3);
GPS_v_target_noise = zeros(length(time.sim_time),3);
GPS_r_chaser_noise = zeros(length(time.sim_time),3);
GPS_v_chaser_noise = zeros(length(time.sim_time),3);

% Old methods of generating noise
%GPS_r_target_noise = GPS_r_STD(1).*randn(length(time.sim_time),1); % Used in thesis
%GPS_v_target_noise = GPS_v_STD(1).*randn(length(time.sim_time),1); % Used in thesis
%GPS_r_chaser_noise = GPS_r_STD(1).*randn(length(time.sim_time),1); % Used in thesis
%GPS_v_chaser_noise = GPS_v_STD(1).*randn(length(time.sim_time),1); % Used in thesis
%GPS_r_target_noise = GPS_r_STD.*randn(length(time.sim_time),3); % for 3x1 noise
%GPS_v_target_noise = GPS_v_STD.*randn(length(time.sim_time),3); % for 3x1 noise
%GPS_r_chaser_noise = GPS_r_STD.*randn(length(time.sim_time),3); % for 3x1 noise
%GPS_v_chaser_noise = GPS_v_STD.*randn(length(time.sim_time),3); % for 3x1 noise

% Generate awgn channel for each x,y,z coordinate of position/velocity
for iSeed = 1:length(r_target_NoiseSeed)
    rng(r_target_NoiseSeed(iSeed), 'v4')
    GPS_r_target_noise(:,iSeed) = GPS_r_STD(iSeed)*randn(length(time.sim_time),1);

    rng(r_chaser_NoiseSeed(iSeed), 'v4')
    GPS_r_chaser_noise(:,iSeed) = GPS_r_STD(iSeed)*randn(length(time.sim_time),1);

    rng(v_target_NoiseSeed(iSeed), 'v4')
    GPS_v_target_noise(:,iSeed) = GPS_v_STD(iSeed)*randn(length(time.sim_time),1);

    rng(v_chaser_NoiseSeed(iSeed), 'v4')
    GPS_v_chaser_noise(:,iSeed) = GPS_v_STD(iSeed)*randn(length(time.sim_time),1);
end

% Add noise to the truth data
GPS_r_target = r_target_ECI + GPS_r_target_noise;
GPS_v_target = v_target_ECI + GPS_v_target_noise;
GPS_r_chaser = r_chaser_ECI + GPS_r_chaser_noise;
GPS_v_chaser = v_chaser_ECI + GPS_v_chaser_noise;

%% Calculate relative position/velocity in LVLH Frame, add noise

% Initialize Measurement Variables
r_clean             = zeros(length(time.sim_time),3);
v_clean             = zeros(length(time.sim_time),3);
h_target            = zeros(length(time.sim_time),1);
theta_clean         = zeros(length(time.sim_time),1);
theta_dot_clean     = zeros(length(time.sim_time),1);
rt_clean            = zeros(length(time.sim_time),1);
rt_dot_clean        = zeros(length(time.sim_time),1);
r_meas              = zeros(length(time.sim_time),3);
v_meas              = zeros(length(time.sim_time),3);
h_target_meas       = zeros(length(time.sim_time),1);
theta_target        = zeros(length(time.sim_time),1);
theta_dot_meas      = zeros(length(time.sim_time),1);
rt_meas             = zeros(length(time.sim_time),1);
rt_dot_meas         = zeros(length(time.sim_time),1);

for i_step = 1:length(time.sim_time)

    % Truth Model Calculations - Calculate relative position/velocity in LVLH
    r_rel = r_chaser_ECI(i_step,:) - r_target_ECI(i_step,:);
    v_rel = v_chaser_ECI(i_step,:) - v_target_ECI(i_step,:);
    [r_clean(i_step,:), v_clean(i_step,:)] = ECI_To_LVLH(r_target_ECI(i_step,:)', r_rel', v_target_ECI(i_step,:)', v_rel');

    % Calculating Target Orbit Parameters (theta, rt)
    [COE] = COEfromRV(r_target_ECI(i_step,:), v_target_ECI(i_step,:), Constants.mu_Earth);
    arg_of_latitude = wrapTo360(COE.TA + COE.AoP);
    theta_clean(i_step) = arg_of_latitude*pi/180;

    rt_clean(i_step) = norm(r_target_ECI(i_step,:));
    h_target(i_step) = norm(cross(r_target_ECI(i_step,:),v_target_ECI(i_step,:)));

    % Perturbations due to J2
%     dtheta_dt_J2 = (3*Constants.J2*Constants.mu_Earth*Constants.R_Earth^2) / (2*COE.ECC*h_target(i_step)*rt_clean(i_step)^3) * ...
%                                 ( (h_target(i_step)^2 / (Constants.mu_Earth*rt_clean(i_step)) ) * cosd(COE.TA)*(3*sind(COE.INC)*sind(COE.INC)*sind(COE.TA+COE.AoP)*sind(COE.TA+COE.AoP)-1) ...
%                                 + (2+COE.ECC*cosd(COE.TA)) * sind(2*(COE.TA+COE.AoP)) *sind(COE.INC)*sind(COE.INC) * sind(COE.TA) );
% 
%     domega_dt_J2 = (3*Constants.J2*Constants.mu_Earth*Constants.R_Earth^2) / (2*COE.ECC*h_target(i_step)*rt_clean(i_step)^3) * ...
%                                 ( (h_target(i_step)^2 / (Constants.mu_Earth*rt_clean(i_step)) ) * cosd(COE.TA)*(1-3*sind(COE.INC)*sind(COE.INC)*sind(COE.TA+COE.AoP)*sind(COE.TA+COE.AoP)) ...
%                                 - (2+COE.ECC*cosd(COE.TA)) * sind(2*(COE.TA+COE.AoP)) *sind(COE.INC)*sind(COE.INC) * sind(COE.TA) ...
%                                 + 2*COE.ECC*cosd(COE.INC)*cosd(COE.INC)*cosd(COE.TA+COE.AoP)*cosd(COE.TA+COE.AoP) );
%     theta_dot_clean(i_step) = h_target(i_step)/(rt_clean(i_step)^2) + dtheta_dt_J2 + domega_dt_J2; % [rad/s]
%     dh_dt_J2 = (3*Constants.J2*Constants.mu_Earth*Constants.R_Earth^2) / (2*rt_clean(i_step)^3) * ...
%                     (sind(COE.INC)*sind(COE.INC)*sind(COE.TA+COE.AoP)*sind(COE.TA+COE.AoP) );

    theta_dot_clean(i_step) = h_target(i_step)/(rt_clean(i_step)^2); % [rad/s]
    rt_dot_clean(i_step)  = (Constants.mu_Earth/h_target(i_step))*COE.ECC*sin(theta_clean(i_step));

    % Measurement Calculations - Calculate relative position/velocity in LVLH
    r_rel_m = GPS_r_chaser(i_step,:) - GPS_r_target(i_step,:);
    v_rel_m = GPS_v_chaser(i_step,:) - GPS_v_target(i_step,:);
    [r_meas(i_step,:), v_meas(i_step,:)] = ECI_To_LVLH(GPS_r_target(i_step,:)', r_rel_m', GPS_v_target(i_step,:)', v_rel_m');

    % Calculating Target Orbit Parameters (theta, rt)
    [COE_m]             = COEfromRV(GPS_r_target(i_step,:), GPS_v_target(i_step,:), Constants.mu_Earth);
    arg_of_latitude     = wrapTo360(COE_m.TA + COE_m.AoP);
    theta_target(i_step)  = arg_of_latitude*pi/180;

    rt_meas(i_step)         = norm(GPS_r_target(i_step,:));
    h_target_meas(i_step)   = norm(cross(GPS_r_target(i_step,:),GPS_v_target(i_step,:)));
    theta_dot_meas(i_step)  = h_target_meas(i_step)/(rt_meas(i_step)^2); % [rad/s]
    rt_dot_meas(i_step)     = (Constants.mu_Earth/h_target_meas(i_step))*COE_m.ECC*sin(theta_target(i_step));
end

% Unwrap angular values to monotonically increasing theta
theta_clean = unwrap(theta_clean);
theta_target = unwrap(theta_target);
GPS_output = [r_meas theta_target v_meas];

% Ensure initial measurement and initial state estimate of theta are
% aligned to the same 2pi revolution
if state0(4,1) - theta_clean(1) > pi
    state0(4) = state0(4,1) - 2*pi;
end

%% Integrate NERM Dynamics from initial state estimate (for comparison)

% Integration of the dynamics equations
[~,NERM_output] = ode45(@(t,y) FFNAV_relative_motion_diff_eqns(t,y,Constants), time.sim_time, state0, ode_options);
% Note: Expect errors in the intial state to lead to inaccurate state estimates over time

% Extract parameters of interest
int_theta       = NERM_output(:,4);
int_theta_dot   = NERM_output(:,9);
int_rt          = NERM_output(:,5);
int_rt_dot      = NERM_output(:,10);

%% Run the EKF

% Assemble the measurement vector Z0
meas_vec = [r_meas theta_target v_meas];
meas0    = [meas_vec(1,:)]';

% Vectorize the initial covariance matrix P0
P0 = [ P0(1,:) P0(2,:) P0(3,:) P0(4,:) P0(5,:) P0(6,:) P0(7,:) P0(8,:) P0(9,:) P0(10,:)]; 

% Vectorize the Initial EKF State (1x117)   
EKF0 = [ state0' P0 meas0' ];

% Specify when measurements are available for processing
flag_is_measurement_available = ones(length(time.sim_time),1);
%flag_is_measurement_available(60:5:end) = 1;

% EKF Simulation Loop
EKF_output = zeros(length(time.sim_time),243);
EKF_output(1,1:117) = EKF0;
for i_step = 2:length(time.sim_time)

    % Define the previous state
    state_pre = EKF_output(i_step-1,:)';

    % Propagate previous state estimates to a priori estimates using dynamics model
    state_priori = FFNAV_EKF_Propagation(state_pre, time.step, Constants.mu_Earth, Q0);

    % Correct a priori state estimates to a posteriori estimates using measurements
    if flag_is_measurement_available(i_step)
        EKF_output(i_step,:) = FFNAV_EKF_Correction(state_priori, meas_vec(i_step,:)', R0);
    else
        EKF_output(i_step,1:110) = state_priori;
    end
end

% =========================================================================
%% Clean-up and save data
if profiler_flag
    profile off
end

clearvars -except   EKF_output   GPS_output     NERM_output             ...
                    EKF_flag     EKF_name       time        t_run       ...
                    r_clean      v_clean        R_Earth     mu          ...
                    r_target_ECI v_target_ECI                           ...
                    r_chaser_ECI v_chaser_ECI                           ...
                    theta_target rt_clean       rt_dot_clean            ...
                    theta_clean  theta_dot_clean                        ...
                    Pv_obs       Pv_smoothed    Pr_theo                 ...
                    P_error      Pdot_error     P0                      ...
                    Q_output     R_output       Q0          R0          ...
                    lambda       NIS            deltaDOD    DOD_true    ...
                    save_flag    post_flag      plot_flag   print_flag  ...
                    pathname     FF_Config      Sim_filename...
                    File_specifier Constants...
                    profiler_flag  cov_flag Q_adapt R_adapt n_start n_end paths input_data;
                    
if save_flag
    file_out = strcat(FF_Config, ['_', EKF_name,...
                            '_',File_specifier]);
    clear save_output
    file_out = fullfile(input_data.paths.output, '\',file_out);
    save(file_out);        
end

% =========================================================================
%% Post Processing  
if post_flag
    tic;
    FFNAV_EKF_PostProcess(file_out)
    t_process = toc;
    fprintf('Post-processing Time = %f \n', t_process)
end

% Plotting
if plot_flag
    FFNAV_EKF_Plotter(file_out)
end

% Computation Profile Viewer
if profiler_flag
    profile viewer
end
fprintf('--------------------------------------------------------------\n')

end