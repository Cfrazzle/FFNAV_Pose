%% Run Attitude Simulation Script
% =========================================================================
% Description: This scripts performs a reorientation manuever for a 
% spacecraft regulation case involving external torques.
%
% Inputs ==================================================================
%
% Outputs =================================================================
%
% Notes ===================================================================
%   quaternions - defined here as: 
%           q = q1_ihat + q2_jhat + q3_khat + q0
%      in some nomenclature, it is common to see the scalar part (q0)
%      referred to as q4
%      
%
% References ==============================================================
%   [1] Fundamentals of Spacecraft Attitude Determination and Control,
%       by Markley and Crassidis (Example 7.1, p)
% 
% Created by: Cory Fraser - JAN 12, 2023
% Last edits: Cory Fraser - FEB 04, 2023
% =========================================================================
%% Premable

import attitude_sim.*
set(0,'DefaultFigureWindowStyle','docked')

clc
close all
clear

%% Time and Inertia Matrix
J_inertia = diag([10000 9000 12000]);
inv_J_inertia = inv(J_inertia);
mu = 398600.4418;

% Initial and Desired Quaternion q = [vec scalar]
Qb_desired0 = [ 0 0 0 1]';

% Orbit Data for STK Scenario Comparison
a           = 7000; 
ecc         = 0.05; 
inc         = 55*pi/180;
big_omega   = 0*pi/180;
w           = 0*pi/180;
big_m       = 0*pi/180;

%% Initial Position/Velocity calculation

% Solve Kepler's Equation
% Initial Guess for Eccentric Anomaly
big_e=big_m;
delta_e=10;eps=1e-15;
max_iter=100;count=1;

while abs(delta_e) > eps
  delta_e=(big_m-(big_e-ecc*sin(big_e)))/(1-ecc*cos(big_e));
  big_e=big_e+delta_e;
  count = count + 1;
  if count == max_iter
      disp(' Maximum Number of Iterations Achieved')
      break
  end
end

% Get Initial Position and Velocity
rmag=a*(1-ecc*cos(big_e));
rp=a*(cos(big_e)-ecc);rq=a*sqrt(1-ecc^2)*sin(big_e);
vp=-sqrt(mu*a)/rmag*sin(big_e);vq=sqrt(mu*a*(1-ecc^2))/rmag*cos(big_e);
c11=cos(big_omega)*cos(w)-sin(big_omega)*sin(w)*cos(inc);
c12=-cos(big_omega)*sin(w)-sin(big_omega)*cos(w)*cos(inc);
c21=sin(big_omega)*cos(w)+cos(big_omega)*sin(w)*cos(inc);
c22=-sin(big_omega)*sin(w)+cos(big_omega)*cos(w)*cos(inc);
c31=sin(w)*sin(inc);
c32=cos(w)*sin(inc);
r1=c11*rp+c12*rq;r2=c21*rp+c22*rq;r3=c31*rp+c32*rq;
v1=c11*vp+c12*vq;v2=c21*vp+c22*vq;v3=c31*vp+c32*vq;

% Initial State Vector for Satellite Orbit
r0=[r1;r2;r3];
v0=[v1;v2;v3];

% Orbit Period
orb_per=2*pi/sqrt(mu)*(a^(3/2));

% STK Scenario Initial Conditions (in Inertial Frame)
pos0 = [6650	0	0]'; % S/C position in inertial frame
vel0 = [0	4.550342	6.498562]'; % Spacecraft velocity in inertial Frame
w0  = [0	-0.068352	0]'*pi/180; % Spacecraft angular rates in body Frame
q0  = [-0.212631	-0.67438	0.212631	0.67438]; % Spacecaft attitude in inertial frame

% Initial Angular Rates and quaternion for Crassidas example case 7.1
%q0 = [ 0.685 0.695 0.153 0.153]';
%q0 = q0/norm(q0);
%w0 = [0.53 0.53 0.053]'*pi/180;

%% Converting between Inertial and Spacecraft Frames

% Calculate orientation of the spacecraft frame, seen in inertial frame
C_LI = rotmat_ECI_to_VVLH(pos0, vel0);
C_IL = C_LI';
Q_IL = my_quaternion.rotmat2quat(C_IL); % [scalar vector]

% Convert inertial s/c quaternion into body-frame
q_att_b0 = my_quaternion.quat_mult(Q_IL,q0);

%% Initial state vector for Attitude Dynamics

time_step = 1;
time_start = 0;
time_stop = 1*orb_per; %300;
time_vec = (time_start:time_step:time_stop)';
n_steps = length(time_vec);
state_vec = zeros(n_steps,7+6);
Qb_desired = zeros(n_steps,4);

% Define initial state vector
state_vec(1,:) = [q_att_b0 w0' pos0' vel0'];

% Controller Gains
flag_use_controller = true;
kp = zeros(1,1);
kd = zeros(1,1);
if flag_use_controller
    kp = 50;
    kd = 500;
end

%% Simulation Loop - Numerical integration of the dynamics equations
for i = 1:n_steps-1

    % Speicify the desired attitude for this time step
    Qb_desired(i,:) = Qb_desired0;

    f1 = time_step * attitude_sim.euler_quaternion_dot(state_vec(i,:), J_inertia,inv_J_inertia,Qb_desired(i,:)',kp,kd);
    f2 = time_step * attitude_sim.euler_quaternion_dot(state_vec(i,:) + 0.5*f1', J_inertia,inv_J_inertia,Qb_desired(i,:)',kp,kd);
    f3 = time_step * attitude_sim.euler_quaternion_dot(state_vec(i,:) + 0.5*f2', J_inertia,inv_J_inertia,Qb_desired(i,:)',kp,kd);
    f4 = time_step * attitude_sim.euler_quaternion_dot(state_vec(i,:) + f3', J_inertia,inv_J_inertia,Qb_desired(i,:)',kp,kd);
    state_vec(i+1,:) = state_vec(i,:) + 1/6*(f1'+2*f2'+2*f3'+f4');
end
Qb_desired(i+1,:) = Qb_desired0;

%% Post-processing
Qb_out = state_vec(:,1:4);
wb_out = state_vec(:,5:7);
r_out = state_vec(:,8:10);
v_out = state_vec(:,11:13);

% Calculate Error Quaternion over time
q_error = my_quaternion.quat_err(Qb_out,Qb_desired);
for iStep = 2:size(q_error,1)

    if sign(q_error(iStep,4)) ~= sign(q_error(iStep-1,4))
        q_error(iStep,4) = -q_error(iStep,4);
    end
end

% Torque command over time
torque = -kp*kron(sign(q_error(:,4)),[1 1 1]).*q_error(:,1:3)-kd*wb_out;
%torque = sign(torque).*min(abs(torque),5);

% Organize Output Data into a single struct
OutputData.time_vec = time_vec;
OutputData.x        = state_vec;
OutputData.qerr     = q_error;
OutputData.torque   = torque;

% Create plots
attitude_sim.create_plots(OutputData)

%% Convert to ECI Frame

Qb_out_I = zeros(size(state_vec,1), 4);
omega_I_out = zeros(size(state_vec,1), 3);
for iStep = 1:size(state_vec,1)

    % Calculate rotation matrix from VVLH frame to inertial frame
    C_LI = attitude_sim.rotmat_ECI_to_VVLH(r_out(iStep,:)', v_out(iStep,:)');
    Q_LI = my_quaternion.rotmat2quat(C_LI); %[q = scalar vector]

    % Compute angular rate of the VVLH frame w.r.t. ECI frame
    omega_LI =  cross(r_out(iStep,:)',v_out(iStep,:)')/norm(r_out(iStep,:)')^2;
    
    % Add body rates to rates of VVLH frame rate converted from ECI to VVLH frame
    omega_I_out(iStep,:) = (C_LI*omega_LI)' + wb_out(iStep,:);

    % Convert body attitude to inertial frame
    Q_mult = my_quaternion.quat_mult(Q_LI, Qb_out(iStep,:)); %Expects/outputs q = [vector scalar]
    Qb_out_I(iStep,:) = my_quaternion.slerp(Q_LI, Q_mult,1); % Expects: q= [scalar vector]
    
    % Continuity of quaternion using Slerp
     if iStep > 1
         Qb_out_I(iStep,:) = my_quaternion.slerp(Qb_out_I(iStep-1,:), Qb_out_I(iStep,:),1);
     end
end

%% Compare with STK Scenario
path_name = 'C:\Users\coryt\OneDrive\Obruta Contract\STK_Ref_Scenarios\';
file_name.two_body = 'TwoBody_Satellite1_Inertial_Position_Velocity_Quaternion_Angular_Rates.csv';
file_name.two_body_J2 = 'TwoBody_wJ2_Satellite1_Inertial_Position_Velocity_Quaternion_Angular_Rates.csv';
file_name.HPOP = 'HPOP_defaults_Satellite1_Inertial_Position_Velocity_Quaternion_Angular_Rates.csv';

compare_case = 'two_body'; % two_body, two_body_J2, HPOP
STK_data = attitude_sim.import_STK_pose_data(fullfile(path_name,file_name.(compare_case)));
STK_data = STK_data(1:n_steps,:);

% Calculate the difference between calculated data and STK data
delta_r = r_out - [STK_data.xkm STK_data.ykm STK_data.zkm];
delta_v = v_out - [STK_data.vxkmsec STK_data.vykmsec STK_data.vzkmsec];
delta_q = Qb_out_I - [STK_data.q1vec STK_data.q2vec STK_data.q3vec STK_data.q4scalar];
delta_w = omega_I_out*180/pi - [STK_data.wxdegsec STK_data.wydegsec STK_data.wzdegsec];

% Position
figure
subplot(3,2,1)
plot(time_vec, r_out(:,1), time_vec, STK_data.xkm)
xlabel('Time (s)')
ylabel('Position - x (km)')
legend('Simulator','STK')
subplot(3,2,3)
plot(time_vec, r_out(:,2),time_vec, STK_data.ykm)
xlabel('Time (s)')
ylabel('Position - y (km)')
subplot(3,2,5)
plot(time_vec, r_out(:,3),time_vec, STK_data.zkm)
xlabel('Time (s)')
ylabel('Position - z (km)')
subplot(3,2,2)
plot(time_vec,delta_r(:,1))
xlabel('Time (s)')
ylabel('Position Error - x (km)')
subplot(3,2,4)
plot(time_vec,delta_r(:,2))
xlabel('Time (s)')
ylabel('Position Error - y (km)')
subplot(3,2,6)
plot(time_vec,delta_r(:,3))
xlabel('Time (s)')
ylabel('Position Error - z (km)')

% Velocity
figure
subplot(3,2,1)
plot(time_vec, v_out(:,1),time_vec, STK_data.vxkmsec)
xlabel('Time (s)')
ylabel('Velocity - x (km/s)')
legend('Simulator','STK')
subplot(3,2,3)
plot(time_vec, v_out(:,2),time_vec, STK_data.vykmsec)
xlabel('Time (s)')
ylabel('Velocity - y (km/s)')
subplot(3,2,5)
plot(time_vec, v_out(:,3),time_vec, STK_data.vzkmsec)
xlabel('Time (s)')
ylabel('Velocity - z (km/s)')
subplot(3,2,2)
plot(time_vec,delta_v(:,1))
xlabel('Time (s)')
ylabel('Velocity Error - x (km/s)')
subplot(3,2,4)
plot(time_vec,delta_v(:,2))
xlabel('Time (s)')
ylabel('Velocity Error - y (km/s)')
subplot(3,2,6)
plot(time_vec,delta_v(:,3))
xlabel('Time (s)')
ylabel('Velocity Error - z (km/s)')

% Angular Rates
figure; hold on
subplot(3,2,1)
plot(time_vec, omega_I_out(:,1)*180/pi,time_vec, STK_data.wxdegsec)
xlabel('Time (s)')
ylabel('Angular Rates - \omega_x (deg/s)')
legend('Simulator','STK')
subplot(3,2,3)
plot(time_vec, omega_I_out(:,2)*180/pi,time_vec, STK_data.wydegsec)
xlabel('Time (s)')
ylabel('Angular Rates - \omega_y (deg/s)')
subplot(3,2,5)
plot(time_vec, omega_I_out(:,3)*180/pi,time_vec, STK_data.wzdegsec)
xlabel('Time (s)')
ylabel('Angular Rates - \omega_z (deg/s)')
subplot(3,2,2)
plot(time_vec,delta_w(:,1))
xlabel('Time (s)')
ylabel('Angular Rate Error - \omega_x (deg/s)')
subplot(3,2,4)
plot(time_vec,delta_w(:,2))
xlabel('Time (s)')
ylabel('Angular Rate Error - \omega_y (deg/s)')
subplot(3,2,6)
plot(time_vec,delta_w(:,3))
xlabel('Time (s)')
ylabel('Angular Rate Error - \omega_z (deg/s)')

% Attitude
figure; hold on
subplot(4,2,1)
plot(time_vec, Qb_out_I(:,1), time_vec, STK_data.q1vec)
xlabel('Time (s)')
ylabel('Attitude - q_1')
legend('Simulator','STK')
subplot(4,2,3)
plot(time_vec, Qb_out_I(:,2), time_vec, STK_data.q2vec)
xlabel('Time (s)')
ylabel('Attitude - q_2')
subplot(4,2,5)
plot(time_vec, Qb_out_I(:,3), time_vec, STK_data.q3vec)
xlabel('Time (s)')
ylabel('Attitude - q_3')
subplot(4,2,7)
plot(time_vec, Qb_out_I(:,4), time_vec, STK_data.q4scalar)
xlabel('Time (s)')
ylabel('Attitude - q_4 Scalar')
subplot(4,2,2)
plot(time_vec, delta_q(:,1))
xlabel('Time (s)')
ylabel('Attitude Error - q_1')
subplot(4,2,4)
plot(time_vec, delta_q(:,2))
xlabel('Time (s)')
ylabel('Attitude Error - q_2')
subplot(4,2,6)
plot(time_vec, delta_q(:,3))
xlabel('Time (s)')
ylabel('Attitude Error - q_3')
subplot(4,2,8)
plot(time_vec, delta_q(:,4))
xlabel('Time (s)')
ylabel('Attitude Error - q_4 Scalar')
