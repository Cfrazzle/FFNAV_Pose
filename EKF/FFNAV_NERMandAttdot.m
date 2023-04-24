function [state_dot] = FFNAV_NERMandAttdot(state, mu)
% FFNAV EKF Dot ===========================================================
% Description: This function defines the nonlinear equations of relative
% motion for formation flying spacecraft. When passed a state vector, this
% function calculates the respective accelerations for the x, y, z, rt and
% theta functions, and returns them to the calling function in a vector.
%
% Inputs:
%   state - The current state vector
%   mu    - Earth's gravitational parameter
%
% Outputs:
%   f       - Diffential equation output accelerations
%
% Created by: Cory Fraser - FALL... 2015
% Last Edits: Cory Fraser - JUL 04, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% State Vector

% Equations of relative orbital motion dynamics
x          = state(1);
y          = state(2);
z          = state(3);
theta      = state(4);
rt         = state(5);
x_dot      = state(6);
y_dot      = state(7);
z_dot      = state(8);
theta_dot  = state(9);
rt_dot     = state(10);

% Equations of relative attitude dynamcis
delta_omega1_L0 = state(11);
delta_omega2_L0 = state(12);
delta_omega3_L0 = state(13);
q1 = state(14);
q2 = state(15);
q3 = state(16);
q4 = state(17);

delta_omega_L0 = [delta_omega1_L0 delta_omega2_L0 delta_omega3_L0];
quaternion_att = [q1 q2 q3 q4];

% Required inputs
I0 = diag([1 1 1]);
I1 = diag([1 1 1]);
N0 = diag([0 0 0]);
N1 = diag([0 0 0]);

% Calculate chaser parameters
rc = sqrt((rt + x)^2 + y^2 + z^2); % Radius of the chaser
omega0_L0 = 0; %TBD % Angular rate of the chaser, expressed in chaser frame

%% Derivatives of the state variables
 
% Relative orbital motion dynamics
%x_dot     = x_dot_priori;
%y_dot     = y_dot_priori;
%z_dot     = z_dot_priori;
%theta_dot = theta_dot_priori;
%rt_dot    = rt_dot_priori;
x_ddot     = x*theta_dot^2 + 2*theta_dot*(y_dot-y*rt_dot/rt)+mu/rt^2-mu*(rt+x)/(rc^3);
y_ddot     = y*theta_dot^2 - 2*theta_dot*(x_dot-x*rt_dot/rt)-mu*y/(rc^3);
z_ddot     = -mu*z/rc^3;
theta_ddot = -2*rt_dot*theta_dot/rt;
rt_ddot    = rt*theta_dot^2-mu/rt^2;

% Construct quaternion multiplication matrix Q_beta from quaternion
Q_quat = [ -q1 -q2 -q3
            q4 -q3 q2
            q3 q4 -q1
            -q2 q1 q4 ];

% Construct rotation matrix D_beta from quaternion (q1-q4)
D_beta = my_quaternion.rotmat(quaternion_att,'frame');

% Angular rate equations
delta_omega_L0_dot = inv(I0)* ...
        (I0*D_beta*inv(I1)*(N1 - cross(D_beta'*(delta_omega_L0 + omega0_L0),I1*D_beta'*(delta_omega_L0 + omega0_L0)) )...
        - cross(I0*omega0_L0,delta_omega_L0) ...
        - (N0 - cross(omega0_L0,I0*omega0_L0)) );

w1_dot = delta_omega_L0_dot(1);
w2_dot = delta_omega_L0_dot(2);
w3_dot = delta_omega_L0_dot(3);

% Attitude equations
quaternion_dot = 0.5*Q_quat*delta_omega_L0; % Verify L0 or L1 will work, equivalent?
q1_dot = quaternion_dot(1); 
q2_dot = quaternion_dot(2);
q3_dot = quaternion_dot(3); 
q4_dot = quaternion_dot(4); 

%% Assembling the derivative of the state vector
state_dot = [ x_dot; y_dot; z_dot; theta_dot; rt_dot; x_ddot; y_ddot; z_ddot; theta_ddot; rt_ddot;...
      w1_dot; w2_dot; w3_dot; q1_dot; q2_dot; q3_dot; q4_dot ];

end
% =========================================================================