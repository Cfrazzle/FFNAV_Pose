function [state_dot] = euler_quaternion_dot(state, J, inv_J, qd, kp, kd)
% Rigid Body Dynamics using Quaternions ===================================
% Description: This function defines Eurler's Equations in quaternion form. 
% When passed a state vector, this function calculates the respective state
% derivatives and returns them to the calling function in a vector.
%
% Inputs:
%   state - Current state vector (quaternion and angular rates)
%   J - Inertia matrix
%   inv_J - Inverse of the inertia matrix
%   qd - Desired attitude quaternion
%   kp - Control Gain - Proportional Component
%   kd - Control Gain - Derivative Component
%
% Outputs:
%   f       - Diffential equation output 
%
% Created by: Cory Fraser - JAN 12, 2023
% Last Edits: Cory Fraser - JAN 12, 2023
% Copyright(c) 2023 by Cory Fraser
% =========================================================================
%% Extract quaternion and angular rates from the state vector

quat = state(1:4)';
omega = state(5:7)';

%% Calculate the error quaternion between the current and desired attitudes

% Function form - different formulation, works row-wise across multiple quaternions
quat_err = my_quaternion.quat_err(quat',qd')';

%% Calculate the Torque Command
%u = -kp*qerr(1:3) - kd*omega;
u = -kp * sign(quat_err(4)) * quat_err(1:3) - kd*omega;
% Try plus or minus
%u = -kp * sign(quat_err(4))*quat_err(1:3) - kd*(1+quat_err(1:3)'*quat_err(1:3))*omega;
%u = -kp * sign(quat_err(4))*quat_err(1:3) - kd*(1-quat_err(1:3)'*quat_err(1:3))*omega;

% Implement torque limit
%u = sign(u).*min(abs(u),5);

%% Attitude Dynamics

% Calculate the cross-produce matrix for angular rates
omega_skew = attitude_sim.skew_sym_matrix(omega);

% Assemble 
om = [-omega_skew omega
      -omega' 0];

% Quaternion Rates - Crassidas, 7.1a (pg 289)
q_dot = 0.5 * om * quat;

% Angular Rates - Crassidas, 7.1b (pg 289)
% J * omega_dot = -omega_skew*J*omega + u;
omega_dot = inv_J*(-omega_skew*J*omega + u);

%% Orbital Mechanics

r_sat = state(8:10)';
v_sat = state(11:13)';
mu_Earth = 398600.4418;
r_ddot = -(mu_Earth/(norm(r_sat))^3)*r_sat;

%% Assemble the State Derivate Vector

state_dot = [q_dot 
     omega_dot 
     v_sat
     r_ddot];

end