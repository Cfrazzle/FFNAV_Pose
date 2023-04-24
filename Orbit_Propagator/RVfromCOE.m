function [R_ECI,V_ECI] = RVfromCOE(mu, COE)
% ========================================================================
% Description: This funciton calculates the initial position and velocity
% vectors for a spacecraft, given the classical orbital elements defined 
% in the initialization of the simulation.
%
% Inputs:
%   mu   -
%   COE. -Structure for classical orbital elements
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

% Calculate semi-latus rectum and radius
p   = COE.SMA*(1-COE.ECC^2);                                                %Defining the semilatus rectum (Eq. 2.31)
r   = p/(1+COE.ECC*cos(COE.TA));                               %Defining the radius (Eq. 2.32)

% Calculate perifocal position/velocity vectors
R_P = [r*cos(COE.TA) r*sin(COE.TA) 0]';                           %(r in perifocal frame - Eq. 2.86) 
V_P = [-sqrt(mu/p)*sin(COE.TA) sqrt(mu/p)*(COE.ECC+cos(COE.TA)) 0]';    %(v in perifocal frame - Eq. 2.87)

%Rotation Matrix from P to ECI Frame
C_IP = [    cos(COE.RAAN)*cos(COE.AOP)-sin(COE.RAAN)*cos(COE.INC)*sin(COE.AOP) -cos(COE.RAAN)*sin(COE.AOP)-sin(COE.RAAN)*cos(COE.INC)*cos(COE.AOP) sin(COE.RAAN)*sin(COE.INC);
            sin(COE.RAAN)*cos(COE.AOP)+cos(COE.RAAN)*cos(COE.INC)*sin(COE.AOP) -sin(COE.RAAN)*sin(COE.AOP)+cos(COE.RAAN)*cos(COE.INC)*cos(COE.AOP) -cos(COE.RAAN)*sin(COE.INC);
                        sin(COE.INC)*sin(COE.AOP)                            sin(COE.INC)*cos(COE.AOP)                       cos(COE.INC)  ];

%Definition of ECI Frame
X_ECI   = [1 0 0];
Y_ECI   = [0 1 0];
Z_ECI   = [0 0 1];
F_ECI   = [X_ECI; Y_ECI; Z_ECI]; %Reference Frame ECI

%r and v vectors in ECI Frame
R_ECI   = F_ECI*C_IP*R_P;
V_ECI   = F_ECI*C_IP*V_P;

%r_mag = norm(R_ECI);
%r_hat = R_ECI/r_mag;

end
% ========================================================================