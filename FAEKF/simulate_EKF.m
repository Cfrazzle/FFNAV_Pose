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

% Load initial parameters for real-world simnulator
InitializeRWSW

% Load initial parameters for GPS constellation
InitializeGPS

% Load initial parameters for on-board software
InitializeRelNavEKF_OBSW

%% Select Simulation Options
time.start     = 0;                           % Start time of the simulation
RWSW_dt        = 1;                           % Time step for the simulation
time.step      = RWSW_dt;
orbit_num      = 2;                           % Number of orbits to simulate
time.stop      = ceil(orbit_num*COE_target.T);    % End time of the simulation
time.sim_time  = time.start:time.step:time.stop;

%Initial Position/Velocity vectors for the spacecraft
[R0_target, V0_target] = RVfromCOE(Constants.mu_Earth, COE_target);
[R0_chaser, V0_chaser] = RVfromCOE(Constants.mu_Earth, COE_chaser);

%% Real-world Simulator

% Propagate target and chaser orbits
ode_options = odeset('relTol', 1e-10, 'AbsTol', 1e-10);
[~,RV_target] = ode45(@(t,y) accelerations_diff_eq(t,y,Constants,PropOptions_target), time.sim_time, [R0_target, V0_target],ode_options);
[~,RV_chaser] = ode45(@(t,y) accelerations_diff_eq(t,y,Constants,PropOptions_chaser), time.sim_time, [R0_chaser, V0_chaser],ode_options);

r_target_ECI = [RV_target(:,1),RV_target(:,2),RV_target(:,3)];
v_target_ECI = [RV_target(:,4),RV_target(:,5),RV_target(:,6)];
r_chaser_ECI = [RV_chaser(:,1),RV_chaser(:,2),RV_chaser(:,3)];
v_chaser_ECI = [RV_chaser(:,4),RV_chaser(:,5),RV_chaser(:,6)];

%% Add measurement noise to ECI position/velocities
GPS_r_target_noise = zeros(length(time.sim_time),3);
GPS_v_target_noise = zeros(length(time.sim_time),3);
GPS_r_chaser_noise = zeros(length(time.sim_time),3);
GPS_v_chaser_noise = zeros(length(time.sim_time),3);

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
%GPS_r_noise = GPS_r_STD(1).*randn(length(time.sim_time),1); % Used in thesis
%GPS_v_noise = GPS_v_STD(1).*randn(length(time.sim_time),1); % Used in thesis
%GPS_r_target_noise = GPS_r_STD.*randn(length(time.sim_time),3); % for 3x1 noise
%GPS_v_target_noise = GPS_v_STD.*randn(length(time.sim_time),3); % for 3x1 noise
GPS_r_target = r_target_ECI + GPS_r_target_noise;
GPS_v_target = v_target_ECI + GPS_v_target_noise;

%GPS_r_noise = GPS_r_STD(1).*randn(length(time.sim_time),1);
%GPS_v_noise = GPS_v_STD(1).*randn(length(time.sim_time),1);
%GPS_r_chaser_noise = GPS_r_STD.*randn(length(time.sim_time),3); % for 3x1 noise
%GPS_v_chaser_noise = GPS_v_STD.*randn(length(time.sim_time),3); % for 3x1 noise
GPS_r_chaser = r_chaser_ECI + GPS_r_chaser_noise;
GPS_v_chaser = v_chaser_ECI + GPS_v_chaser_noise;

% Initialize RWSW Simulator Variables
r_clean = zeros(length(time.sim_time),3);
v_clean = zeros(length(time.sim_time),3);
h_target = zeros(length(time.sim_time),1);
theta_clean = zeros(length(time.sim_time),1);
theta_dot_clean = zeros(length(time.sim_time),1);
rt_clean = zeros(length(time.sim_time),1);
rt_dot_clean = zeros(length(time.sim_time),1);
r_meas = zeros(length(time.sim_time),3);
v_meas = zeros(length(time.sim_time),3);
h_target_meas = zeros(length(time.sim_time),1);
theta_target = zeros(length(time.sim_time),1);
theta_dot_meas = zeros(length(time.sim_time),1);
rt_meas = zeros(length(time.sim_time),1);
rt_dot_meas = zeros(length(time.sim_time),1);

%% RWSW Simulator - calculate relative position and velocities in LVLH Frame
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

    % Due to J2
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

%% Run the EKF

% Assemble the measurement vector
meas_vec = [r_meas theta_target v_meas];
meas0 = [meas_vec(1,:)]';

P0= [ P0(1,:) P0(2,:) P0(3,:) P0(4,:) P0(5,:) ...
      P0(6,:) P0(7,:) P0(8,:) P0(9,:) P0(10,:)]; 

% Initial EKF Vector (1x117)   
if state0(4,1) - theta_clean(1) > pi
    state0(4) = state0(4,1) - 2*pi;
end
EKF0 = [ state0' P0 meas0' ];

% EKF Simulation Loop
EKF_output = zeros(length(time.sim_time),243);
EKF_output(1,1:117) = EKF0;
for i_step = 2:length(time.sim_time)

    % Data from Previous Time Step
    state_pre = EKF_output(i_step-1,:)';

    % Propagate state estimates using dynamics model
    state_priori = FFNAV_FAEKF_Propagation(state_pre, time.step, Constants.mu_Earth, Q0);

    % Correct state estimates using measurements
    EKF_output(i_step,:) = FFNAV_FAEKF_Correction(state_priori, meas_vec(i_step,:)', R0);
end

% NERM Integration for target orbit parameters (theta, rt)
[~,NERM_output] = ode45(@(t,y) FFNAV_relative_motion_diff_eqns(t,y,Constants), time.sim_time, state0, ode_options);
%int_theta       = wrapTo2Pi(NERM_output(:,4));
int_theta       = NERM_output(:,4);
int_theta_dot   = NERM_output(:,9);
int_rt          = NERM_output(:,5);
int_rt_dot      = NERM_output(:,10);

%% Post-processing
%{
% Plot ECI Positions - Truth
figure; hold on
plot3(r_target_ECI(:,1),r_target_ECI(:,2),r_target_ECI(:,3))
plot3(r_chaser_ECI(:,1),r_chaser_ECI(:,2),r_chaser_ECI(:,3))
title('Truth ECI Positions')

% Plot ECI Positions - Measured
figure; hold on
plot3(GPS_r_target(:,1),GPS_r_target(:,2),GPS_r_target(:,3))
plot3(GPS_r_chaser(:,1),GPS_r_chaser(:,2),GPS_r_chaser(:,3))
title('Measured ECI Positions')

figure; hold on
plot(time.sim_time, GPS_r_target(:,1)-r_target_ECI(:,1))
plot(time.sim_time, GPS_r_target(:,2)-r_target_ECI(:,2))
plot(time.sim_time, GPS_r_target(:,3)-r_target_ECI(:,3))
legend('x', 'y', 'z')
title('Measured ECI Positions - Errors')

% Plot LVLH Positions
figure; hold on
plot3(r_clean(:,1),r_clean(:,2),r_clean(:,3))
xlabel('x position (m)')
axis([-500 500 -1000 0 -200 200])
ylabel('y position (m)')
zlabel('z position (m)')
scatter3(0,0,0, 'k', 'filled')
view(3)

% Target orbit Parameters
figure;
subplot(2,2,1)
plot(time.sim_time,theta_clean*180/pi); hold on
if exist('int_theta', 'var')
    plot(time.sim_time,int_theta*180/pi)
end
ylabel('\theta (deg)')
subplot(2,2,3)
plot(time.sim_time, rt_clean); hold on
if exist('int_rt', 'var')
    plot(time.sim_time,int_rt)
end
ylabel('r_t (km)')
subplot(2,2,2)
plot(time.sim_time,theta_dot_clean*180/pi); hold on
if exist('int_theta_dot', 'var')
    plot(time.sim_time,int_theta_dot*180/pi)
end
plot(time.sim_time(2:end),wrapTo2Pi(diff(theta_clean))/time.step*180/pi)
ylabel('\theta-dot (deg)')
subplot(2,2,4)
plot(time.sim_time, rt_dot_clean); hold on
if exist('int_rt_dot', 'var')
    plot(time.sim_time,int_rt_dot)
end
plot(time.sim_time(2:end),diff(rt_clean)/time.step)
legend('Calculated', 'Integrated Eqns', 'Calc via diff()')
ylabel('r_t-dot (km/s)')

% Plot clean measurements
figure; hold on
subplot(5,2,1)
plot(time.sim_time,r_clean(:,1))    
ylabel('x (m)')
subplot(5,2,3)
plot(time.sim_time,r_clean(:,2))
ylabel('y (m)')
subplot(5,2,5)
plot(time.sim_time,r_clean(:,3))
ylabel('z (m)')
subplot(5,2,2)
plot(time.sim_time,v_clean(:,1))    
ylabel('x-dot (m/s)')
subplot(5,2,4)
plot(time.sim_time,v_clean(:,2))
ylabel('y-dot (m/s)')
subplot(5,2,6)
plot(time.sim_time,v_clean(:,3))
ylabel('z-dot (m/s)')
subplot(5,2,7)
plot(time.sim_time,theta_clean*180/pi); hold on
if exist('int_theta', 'var')
    plot(time.sim_time,int_theta*180/pi)
    legend('Calculated', 'Integrated Eqns')
end
ylabel('\theta (deg)')
subplot(5,2,9)
plot(time.sim_time, rt_clean); hold on
if exist('int_rt', 'var')
    plot(time.sim_time,int_rt)
end
ylabel('r_t (km)')
subplot(5,2,8)
plot(time.sim_time,theta_dot_clean*180/pi); hold on
if exist('int_theta_dot', 'var')
    plot(time.sim_time,int_theta_dot*180/pi)
end
ylabel('\theta-dot (deg)')
subplot(5,2,10)
plot(time.sim_time, rt_dot_clean); hold on
if exist('int_rt_dot', 'var')
    plot(time.sim_time,int_rt_dot)
end
ylabel('r_t-dot (km/s)')

% LVLH Position Errors
figure; hold on
subplot(3,1,1)
plot(time.sim_time, r_clean(:,1))
ylabel('x')
subplot(3,1,2)
plot(time.sim_time, r_meas(:,1))
ylabel('x measured')
subplot(3,1,3)
plot(time.sim_time, r_clean(:,1)-r_meas(:,1))
ylabel('x error')

% LVLH Velocity Errors
figure; hold on
subplot(3,1,1)
plot(time.sim_time, v_clean(:,1))
ylabel('x-dot')
subplot(3,1,2)
plot(time.sim_time, v_meas(:,1))
ylabel('x-dot measured')
subplot(3,1,3)
plot(time.sim_time, v_clean(:,1)-v_meas(:,1))
ylabel('x-dot error')
%}

% =========================================================================
%% Clean-up and save data
profiler_flag = false;
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
                    
% Post Processing Options
save_flag     = true;     % Save data to a .mat file
post_flag     = true;     % Post-process the simulation data
    n_start   = 1; %floor(time_end*0.5); 
    n_end     = time.stop+1;  
    cov_flag  = 'on';    % Perform analysis of covariances
plot_flag     = false;     % Create output plots
print_flag    = 'off';    % Print plots to .eps files
profiler_flag = false;    % Show profile simulink performance

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
    FFNAV_FAEKF_PostProcess(file_out)
    t_process = toc;
    fprintf('Post-processing Time = %f \n', t_process)
end

% Plotting
if plot_flag
    FFNAV_FAEKF_Plotter(file_out)
end

% Computation Profile Viewer
if profiler_flag
    profile viewer
end

%load gong
%sound(y,Fs) 

fprintf('--------------------------------------------------------------\n')
%}
% =========================================================================
end