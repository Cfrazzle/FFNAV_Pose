function [] = FFNAV_EKF_MakeParams(paths)    
%% FFNAV EKF Parameters ==================================================
% ========================================================================
% Description: This script constructs the initial conditions for the 
% extended Kalman filtering algorithm.
%
% References:
%   - 
%   - 
%
% Created by:  Cory Fraser - JUL 10, 2017
% Latest Edit: Cory Fraser - OCT 21, 2018
% Copyright(c) 2018 by Cory Fraser
%
%% Initialization and parameter definitions
fprintf('\n--------------------------------------------------------------\n')
fprintf(' \t\t FFNAV EXTENDED KALMAN FILTER PARAMETERS \t\t             \n')
fprintf('--------------------------------------------------------------\n')

% Define path to parameters/output data
%data_path = fileparts('C:\Users\Cory\Desktop\FFNAV Data\');
%addpath('C:\Users\Cory\Desktop\FFNAV Data\')

% Check for Astrodynamic Constants
if isempty(dir([paths.output, '\AstroConstants.mat']))
    fprintf('Astrodynamic constants were not found in this directory...\n')
    fprintf('Generating constants... \n')
    
    %Create orbital mechanics constants file
    OrbitalMechanicsConstants
else
    load('AstroConstants.mat')
end
fprintf('Astrodynamic constants  : %s \n', ConstantsModel)

mu = Constants.mu_Earth;
loaded = load('AstroConstants.mat', 'deg2rad');
deg2rads = loaded.deg2rad;

%Starting time for the simulation
t0 = 0;

%% Select Formation Configuration
FF_Config = 'PRISMA';
    
    % PRISMA     = PRISMA Mission (LEO, low eccentricty)
    % PROBA-3    = PROBA-3 Mission (LEO, high eccentricity)
    % PEOinLEO   = Modified PRISMA, with inclination change (LEO,low e,delta i)
    % BUSSE      = Busse Article   (LEO, 
    % Kuiack33   = Figure 3.3 in Thesis (In-plane Elliptical Orbit, low e)
    % Kuiack34   = Figure 3.4 in Thesis (In-plane Elliptical Orbit, e=0.2)
    % Kuiack35   = Figure 3.5 in Thesis (In-plane Elliptical Orbit, e=0.4)
    % Kuiack511  = Figure 5.11 in Thesis (Highly Elliptical Orbit, e=0.8)
    % Kuiack_J2invariant = Figure 5.15 in Thesis ( Elliptical Orbit, e=0.8)
    % LF_Hold    = Leader-follower Along-track holding position
    
switch FF_Config
    
    % =====================================================================
    case 'PRISMA'
        
        %Initial COE for Target Spacecraft
        %MeanAnomaly_target   = -1.09332379*deg2rad;  %True Anomaly [rad]
        theta0_target   = 358.90349028*deg2rads;
        a_target        = 7087.29755686634e3;   %Semi-major Axis [m]
        e_target        = 0.00145443;           %Eccentricity
        i_target        = 98.18528613*deg2rads;  %Inclination [rad]
        RAAN_target     = 189.8914*deg2rads;     %RAAN [rad]
        AoP_target      = 1.097451382*deg2rads;  %Argument of Perigee [rad]
        n_target        = sqrt(mu/a_target^3);  %Mean orbital motion [rad/s]
        T_target        = 2*pi/n_target;        %Orbital Period [s]
        
        COE_target = [theta0_target a_target e_target ...
                        i_target RAAN_target AoP_target];
        
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        theta0_chaser   = 0.00*deg2rads;         %True Anomaly [rad]
        a_chaser        = 7087.29767733179e3;   %Semi-major Axis [m]
        e_chaser        = 0.00145908;           %Eccentricity - delta e
        i_chaser        = 98.18466676*deg2rads;  %Inclination [rad]
        RAAN_chaser     = 189.8908602*deg2rads;  %Right-ascension of the Ascending Node [rad]
        AoP_chaser      = 0.00*deg2rads;            %Argument of Perigee [rad]
        n_chaser        = sqrt(mu/a_chaser^3);  %Mean orbital motion [rad/s]
        T_chaser        = 2*pi/n_chaser;        %Orbital Period [s]
        
        COE_chaser = [theta0_chaser a_chaser e_chaser ...
                        i_chaser RAAN_chaser AoP_chaser];
        
        %Additional parameters for initial conditions
        p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot
        
        %Initial Conditions for In-Plane Elliptical Formation in LVLH Frame
        x0_chaser       = -34.718875426019130;    %Initial x-position of chaser(m)
        y0_chaser       = -1.068178082293263e+02; %Initial y-position of chaser(m)
        z0_chaser       = 65.994822637671480;     %Initial z-position of chaser(m)
        theta0          = theta0_target;          %Initial true anomaly of target (rad)
        rt0             = rt0_target;             %Initial radius of target (m)
        
        x_dot0_chaser   = 0.208732017363154;  %Initial x-velocity of chaser(m/s)
        y_dot0_chaser   = 0.073700404718392;  %Initial y-velocity of chaser(m)
        z_dot0_chaser   = 0.081187656642932;  %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;  %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;     %Initial change in radius (m/s)
        fprintf('Formation Configuration : %s \n', FF_Config)   
    
        % Target (Tango) Spacecraft Properties for Orbit Propagator
        PropOptions_target.M_Sat   = 42.5;    %Spacecraft mass
        PropOptions_target.A_Drag  = 0.38;    %Area exposed to drag forces 
        PropOptions_target.Cd      = 2.25;    %Coefficient of Drag
        PropOptions_target.A_Solar = 0.55;    %Area exposed to solar radiation
        PropOptions_target.Cr      = 1.2;     %Coefficient of Reflection

        % Chaser (Mango) Spacecraft Properties for Orbit Propagator
        PropOptions_chaser.M_Sat   = 154.4;   %Spacecraft mass
        PropOptions_chaser.A_Drag  = 2.75;%1.3;     %Area exposed to drag forces
        PropOptions_chaser.Cd      = 2.5;     %Coefficient of Drag
        PropOptions_chaser.A_Solar = 2.5;     %Area exposed to solar radiation
        PropOptions_chaser.Cr      = 1.32;    %Coefficient of Reflection
    
    % =====================================================================
    case 'PEOinLEO'
        
        %Initial COE for Target Spacecraft
        theta0_target   = 0;                    %True Anomaly [rad]
        a_target        = 7500e3;               %Semi-major Axis [m]
        e_target        = 0.1;                  %Eccentricity
        i_target        = 98.1877*deg2rads;       %Inclination [rad]
        RAAN_target     = 189.8914*deg2rads;      %RAAN [rad]
        AoP_target      = 1.0938*deg2rads;        %Argument of Perigee [rad]
        n_target        = sqrt(mu/a_target^3);  %Mean orbital motion [rad/s]
        T_target        = 2*pi/n_target;        %Orbital Period [s]
        
        COE_target = [theta0_target a_target e_target ...
                        i_target RAAN_target AoP_target];
        
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        theta0_chaser   = theta0_target;        %True Anomaly [rad]
        a_chaser        = a_target;             %Semi-major Axis [m]
        e_chaser        = e_target + 0.00005;   %Eccentricity - delta e
        i_chaser        = i_target-0.01*deg2rads;   %Inclination [rad]
        RAAN_chaser     = RAAN_target;          %Right-ascension of the Ascending Node [rad]
        AoP_chaser      = AoP_target;           %Argument of Perigee [rad]
        n_chaser        = sqrt(mu/a_chaser^3);  %Mean orbital motion [rad/s]
        T_chaser        = 2*pi/n_chaser;        %Orbital Period [s]
        
        COE_chaser = [theta0_chaser a_chaser e_chaser ...
                        i_chaser RAAN_chaser AoP_chaser];
        
        %Additional parameters for initial conditions
        p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot
        
        %Initial Conditions for In-Plane Elliptical Formation in LVLH Frame
        x0_chaser       = -3.750000374605406e+02;          %Initial x-position of chaser(m)
        y0_chaser       =   -0.00196;             %Initial y-position of chaser(m)
        z0_chaser       =  -22.48775;             %Initial z-position of chaser(m)
        theta0          = theta0_target; %Initial true anomaly of target (rad)
        rt0             = rt0_target;    %Initial radius of target (m)
        x_dot0_chaser   = -4.68572441684328e-06;%Initial x-velocity of chaser(m/s)
        y_dot0_chaser   = 0.854695440915271;    %Initial y-velocity of chaser(m)
        z_dot0_chaser   = -1.40647980533729;    %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;    %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;       %Initial change in radius (m/s)
        
        fprintf('Formation Configuration : %s \n', FF_Config)   

        %Using PRISMA Spacecraft Properties
        % Target (Tango) Spacecraft Properties for Orbit Propagator
        PropOptions_target.M_Sat   = 42.5;    %Spacecraft mass
        PropOptions_target.A_Drag  = 0.38;    %Area exposed to drag forces 
        PropOptions_target.Cd      = 2.25;    %Coefficient of Drag
        PropOptions_target.A_Solar = 0.55;    %Area exposed to solar radiation
        PropOptions_target.Cr      = 1.2;     %Coefficient of Reflection

        % Chaser (Mango) Spacecraft Properties for Orbit Propagator
        PropOptions_chaser.M_Sat   = 154.4;   %Spacecraft mass
        PropOptions_chaser.A_Drag  = 1.3;     %Area exposed to drag forces
        PropOptions_chaser.Cd      = 2.5;     %Coefficient of Drag
        PropOptions_chaser.A_Solar = 2.5;     %Area exposed to solar radiation
        PropOptions_chaser.Cr      = 1.32;    %Coefficient of Reflection
    
    % =====================================================================
    case 'Busse'
        
        %Definitions for RWOP - Initial COE for Target Spacecraft
        theta0_target   = 135*deg2rads;           %True Anomaly (rad)
        a_target        = 7017.9957e3;          %Semi-major Axis (m)
        e_target        = 0.005;                %Eccentricity
        i_target        = 28.5*deg2rads;          %Inclination (rad)
        RAAN_target     = 0;                    %Right-ascension of the Ascending Node (rad)
        AoP_target      = -88.76*deg2rads;        %Argument of Perigee (rad)
        n_target        = sqrt(mu/a_target^3);  %Mean orbital motion of the target spacecraft
        T_target        = 2*pi/n_target;        %Orbital Period [s]
        
        COE_target = [theta0_target a_target e_target ...
                        i_target RAAN_target AoP_target];
                    
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        theta0_chaser   = 136.13*deg2rads;        %True Anomaly (rad)
        a_chaser        = 7017.9886e3;          %Semi-major Axis (m)
        e_chaser        = 0.005103;             %Eccentricity - delta e
        i_chaser        = i_target;             %Inclination (rad)
        RAAN_chaser     = -1.13*deg2rads;         %Right-ascension of the Ascending Node (rad)
        AoP_chaser      = AoP_target;           %Argument of Perigee (rad)
        
        COE_chaser = [theta0_chaser a_chaser e_chaser ...
                        i_chaser RAAN_chaser AoP_chaser];
        
        %Additional parameters for initial conditions
        p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot
               
        % Known Initial Conditions for In-Plane Elliptical Formation
        x0_chaser       = 0.832187577110917e3;  %Initial x-position of chaser(m)
        y0_chaser       = 16.988800024893138e3; %Initial y-position of chaser(m)
        x_dot0_chaser   = -0.017771885771207e3; %Initial x-velocity of chaser(m/s)
        theta0          = 2.356194490192345e3;  %Initial true anomaly of target (rad)
        rt0             = rt0_target;           %Initial radius of target (m)
        y_dot0_chaser   = -0.001206063704430e3; %Initial y-velocity of chaser(m)
        z0_chaser       = 45.312543744836360e3; %Initial z-position of chaser(m)
        z_dot0_chaser   = -0.051399538490091e3; %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;    %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;       %Initial change in radius (m/s)
        fprintf('Formation Configuration : %s \n', FF_Config)   

    % ===================================================================== 
    case 'PROBA-3'
        %References: 
            %Kuiack, Spacecraft Formation Guidance and Control on J2-
                %Perturbed Eccentric Orbits, MASc Thesis, 2018
            %J. S. Llorente, A. Agenjo, et al. PROBA-3: Precise Formation 
                %Flying Demonstration Mission. 
                %Acta Astronautica, 82(1):388-46, 2013.
            %Peters, Mission Analysis for PROBA-3 Nominal Operations
        
        %Definitions for RWOP - Initial COE for Target Spacecraft
        theta0_target   = 0;                    %True Anomaly (rad)
        a_target        = 36943e3;              %Semi-major Axis (m)
        e_target        = 0.8111;               %Eccentricity
        i_target        = 59*deg2rads;            %Inclination (rad)
        RAAN_target     = 84*deg2rads;            %Right-ascension of the Ascending Node (rad)
        AoP_target      = 188*deg2rads;           %Argument of Perigee (rad)
        n_target        = sqrt(mu/a_target^3);  %Mean orbital motion of the target spacecraft
        T_target        = 2*pi/n_target;        %Orbital Period [s]
        
       COE_target = [theta0_target a_target e_target ...
                        i_target RAAN_target AoP_target];        
        
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        theta0_chaser   = theta0_target;            %True Anomaly (rad)
        a_chaser        = a_target;                 %Semi-major Axis (m)
        e_chaser        = e_target + 0.000005;      %Eccentricity + delta e
        i_chaser        = i_target;                 %Inclination (rad)
        RAAN_chaser     = RAAN_target;              %RAAN (rad)
        AoP_chaser      = AoP_target;               %Argument of Perigee (rad)
        
        COE_chaser = [theta0_chaser a_chaser e_chaser ...
                        i_chaser RAAN_chaser AoP_chaser];
        
        %Additional parameters for initial conditions
        p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot
                
        % Known Initial Conditions for In-Plane Elliptical Formation
        x0_chaser       = -184.715000001348;     %Initial x-position of chaser(m)
        y0_chaser       = 1.29855237673837e-10;  %Initial y-position of chaser(m)
        z0_chaser       = -9.72395497456091e-11; %Initial z-position of chaser(m)
        theta0          = theta0_target;         %Initial true anomaly of target (rad)
        rt0             = rt0_target;            %Initial radius of target (m)

        x_dot0_chaser   = 1.79981030079546e-13;  %Initial x-velocity of chaser(m/s)
        y_dot0_chaser   = 0.417862035614425;     %Initial y-velocity of chaser(m)
        z_dot0_chaser   = -2.43138842392909e-14; %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;     %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;        %Initial change in radius (m/s)
        fprintf('Formation Configuration : %s \n', FF_Config)   

        % Target (Occulter) Spacecraft Properties for Orbit Propagator
        PropOptions_target.M_Sat   = 211;     %Spacecraft mass
        PropOptions_target.A_Drag  = 1.77;    %Area exposed to drag forces
        PropOptions_target.Cd      = 2.5;        %Coefficient of Drag
        PropOptions_target.A_Solar = 1.77;      %Area exposed to solar radiation
        PropOptions_target.Cr      = 1.9;        %Coefficient of Reflection

        % Chaser (Cronograph) Spacecraft Properties for Orbit Propagator
        PropOptions_chaser.M_Sat   = 339;     %Spacecraft mass
        PropOptions_chaser.A_Drag  = 3.34;    %Area exposed to drag forces
        PropOptions_chaser.Cd      = 2.5;        %Coefficient of Drag
        PropOptions_chaser.A_Solar = 3.34;      %Area exposed to solar radiation
        PropOptions_chaser.Cr      = 1.29;        %Coefficient of Reflection
    
    % =====================================================================       
    case 'Kuiack33'
        
        %Definitions for RWOP - Initial COE for Target Spacecraft
        theta0_target   = 0;                    %True Anomaly (rad)
        a_target        = 7106.14e3;              %Semi-major Axis (m)
        e_target        = 0.05;               %Eccentricity
        i_target        = 98.3*deg2rads;            %Inclination (rad)
        RAAN_target     = 270*deg2rads;            %Right-ascension of the Ascending Node (rad)
        AoP_target      = 0*deg2rads;           %Argument of Perigee (rad)
        n_target        = sqrt(mu/a_target^3);  %Mean orbital motion of the target spacecraft
        T_target        = 2*pi/n_target;        %Orbital Period [s]
        
       COE_target = [theta0_target a_target e_target ...
                        i_target RAAN_target AoP_target];        
        
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        theta0_chaser   = theta0_target;            %True Anomaly (rad)
        a_chaser        = a_target;                 %Semi-major Axis (m)
        e_chaser        = e_target + 0.001; %Eccentricity + delta e
        i_chaser        = i_target;                 %Inclination (rad)
        RAAN_chaser     = RAAN_target;              %RAAN (rad)
        AoP_chaser      = AoP_target;               %Argument of Perigee (rad)
        
        COE_chaser = [theta0_chaser a_chaser e_chaser ...
                        i_chaser RAAN_chaser AoP_chaser];
        
        %Additional parameters for initial conditions
        p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot
        
        %%Initializing parameters for EKF propagation
        
        % Known Initial Conditions for In-Plane Elliptical Formation
        delta_e         = e_chaser-e_target;    %Difference in eccentricity
        x0_chaser       = -a_target*delta_e;    %Initial x-position of chaser(m)        
        y0_chaser       = 0e3;                             %Initial y-position of chaser(m)
        z0_chaser       = 0e3;                             %Initial z-position of chaser(m)
        theta0          = theta0_target;        %Initial true anomaly of target (rad)
        rt0             = rt0_target;           %Initial radius of target (m)
        
        x_dot0_chaser   = y0_chaser*n_target/2;        %Initial x-velocity of chaser(m/s)
        %y_dot0_chaser   = -2*n_target*x0_chaser;       %Initial y-velocity of chaser(m)
        y_dot0_chaser   = 16.1861188067368
        z_dot0_chaser   = 0;                    %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;    %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;       %Initial change in radius (m/s)
        fprintf('Formation Configuration : %s \n', FF_Config)  

    % =====================================================================       
    case 'Kuiack34'
        
        %Definitions for RWOP - Initial COE for Target Spacecraft
        theta0_target   = 0;                    %True Anomaly (rad)
        a_target        = 9000e3;              %Semi-major Axis (m)
        e_target        = 0.2;               %Eccentricity
        i_target        = 98.3*deg2rads;            %Inclination (rad)
        RAAN_target     = 270*deg2rads;            %Right-ascension of the Ascending Node (rad)
        AoP_target      = 0*deg2rads;           %Argument of Perigee (rad)
        n_target        = sqrt(mu/a_target^3);  %Mean orbital motion of the target spacecraft
        T_target        = 2*pi/n_target;        %Orbital Period [s]
        
       COE_target = [theta0_target a_target e_target ...
                        i_target RAAN_target AoP_target];        
        
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        theta0_chaser   = theta0_target;            %True Anomaly (rad)
        a_chaser        = a_target;                 %Semi-major Axis (m)
        e_chaser        = e_target + 0.001; %Eccentricity + delta e
        i_chaser        = i_target;                 %Inclination (rad)
        RAAN_chaser     = RAAN_target;              %RAAN (rad)
        AoP_chaser      = AoP_target;               %Argument of Perigee (rad)
        
        COE_chaser = [theta0_chaser a_chaser e_chaser ...
                        i_chaser RAAN_chaser AoP_chaser];
        
        %Additional parameters for initial conditions
        p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot
        
        %%Initializing parameters for EKF propagation
        
        % Known Initial Conditions for In-Plane Elliptical Formation
        delta_e         = e_chaser-e_target;    %Difference in eccentricity
        x0_chaser       = -a_target*delta_e;    %Initial x-position of chaser(m)
        y0_chaser       = 0e3;                  %Initial y-position of chaser(m)
        y0_chaser       = 0e3;                  %Initial y-position of chaser(m)
        z0_chaser       = 0e3;                  %Initial z-position of chaser(m)
        theta0          = theta0_target;        %Initial true anomaly of target (rad)
        rt0             = rt0_target;           %Initial radius of target (m)
        
        x_dot0_chaser   = y0_chaser*n_target/2; %Initial x-velocity of chaser(m/s)
        y_dot0_chaser   = 18.6848143678902;     %Initial y-velocity of chaser(m)
        z_dot0_chaser   = 0;                    %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;    %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;       %Initial change in radius (m/s)
        fprintf('Formation Configuration : %s \n', FF_Config) 
        
    % =====================================================================       
    case 'Kuiack35'
        
        %Definitions for RWOP - Initial COE for Target Spacecraft
        theta0_target   = 0;                    %True Anomaly (rad)
        a_target        = 12000e3;              %Semi-major Axis (m)
        e_target        = 0.4;               %Eccentricity
        i_target        = 98.3*deg2rads;            %Inclination (rad)
        RAAN_target     = 270*deg2rads;            %Right-ascension of the Ascending Node (rad)
        AoP_target      = 0*deg2rads;           %Argument of Perigee (rad)
        n_target        = sqrt(mu/a_target^3);  %Mean orbital motion of the target spacecraft
        T_target        = 2*pi/n_target;        %Orbital Period [s]
        
       COE_target = [theta0_target a_target e_target ...
                        i_target RAAN_target AoP_target];        
        
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        theta0_chaser   = theta0_target;            %True Anomaly (rad)
        a_chaser        = a_target;                 %Semi-major Axis (m)
        e_chaser        = e_target + 0.001; %Eccentricity + delta e
        i_chaser        = i_target;                 %Inclination (rad)
        RAAN_chaser     = RAAN_target;              %RAAN (rad)
        AoP_chaser      = AoP_target;               %Argument of Perigee (rad)
        
        COE_chaser = [theta0_chaser a_chaser e_chaser ...
                        i_chaser RAAN_chaser AoP_chaser];
        
        %Additional parameters for initial conditions
        p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot
        
        %%Initializing parameters for EKF propagation
        
        % Known Initial Conditions for In-Plane Elliptical Formation
        delta_e         = e_chaser-e_target;    %Difference in eccentricity
        x0_chaser       = -a_target*delta_e;    %Initial x-position of chaser(m)        
        y0_chaser       = 0e3;                  %Initial y-position of chaser(m)
        z0_chaser       = 0e3;                  %Initial z-position of chaser(m)
        theta0          = theta0_target;        %Initial true anomaly of target (rad)
        rt0             = rt0_target;           %Initial radius of target (m)
        
        x_dot0_chaser   = y0_chaser*n_target/2;        %Initial x-velocity of chaser(m/s)
        %y_dot0_chaser   = -2*n_target*x0_chaser;       %Initial y-velocity of chaser(m)
        y_dot0_chaser   = 25.1647559803533;
        z_dot0_chaser   = 0;                    %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;    %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;       %Initial change in radius (m/s)
        fprintf('Formation Configuration : %s \n', FF_Config) 
        
        % =====================================================================       
        case 'Kuiack511'
        
        %Definitions for RWOP - Initial COE for Target Spacecraft
        theta0_target   = 0;                    %True Anomaly (rad)
        a_target        = 35000e3;              %Semi-major Axis (m)
        e_target        = 0.8;               %Eccentricity
        i_target        = 59*deg2rads;            %Inclination (rad)
        RAAN_target     = 84*deg2rads;            %Right-ascension of the Ascending Node (rad)
        AoP_target      = 20*deg2rads;           %Argument of Perigee (rad)
        n_target        = sqrt(mu/a_target^3);  %Mean orbital motion of the target spacecraft
        T_target        = 2*pi/n_target;        %Orbital Period [s]
        
       COE_target = [theta0_target a_target e_target ...
                        i_target RAAN_target AoP_target];        
        
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        theta0_chaser   = theta0_target;        %True Anomaly (rad)
        a_chaser        = a_target;             %Semi-major Axis (m)
        e_chaser        = e_target + 0.00001;   %Eccentricity + delta e
        i_chaser        = i_target;             %Inclination (rad)
        RAAN_chaser     = RAAN_target;          %RAAN (rad)
        AoP_chaser      = AoP_target;           %Argument of Perigee (rad)
        
        COE_chaser = [theta0_chaser a_chaser e_chaser ...
                        i_chaser RAAN_chaser AoP_chaser];
        
        %Additional parameters for initial conditions
        p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot
        
        %%Initializing parameters for EKF propagation
        
        % Known Initial Conditions for In-Plane Elliptical Formation
        delta_e         = e_chaser-e_target;    %Difference in eccentricity
        x0_chaser       = -a_target*delta_e;    %Initial x-position of chaser(m)        
        y0_chaser       = 0e3;                  %Initial y-position of chaser(m)
        z0_chaser       = 0e3;                  %Initial z-position of chaser(m)
        theta0          = theta0_target;        %Initial true anomaly of target (rad)
        rt0             = rt0_target;           %Initial radius of target (m)
        
        x_dot0_chaser   = y0_chaser*n_target/2;        %Initial x-velocity of chaser(m/s)
        %y_dot0_chaser   = -2*n_target*x0_chaser;       %Initial y-velocity of chaser(m)
        y_dot0_chaser   = 0.787439601286276;
        z_dot0_chaser   = 0;                    %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;    %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;       %Initial change in radius (m/s)
        fprintf('Formation Configuration : %s \n', FF_Config) 
        
        % =====================================================================       
        case 'Kuiack_J2invariant'
        
        %Definitions for RWOP - Initial COE for Target Spacecraft
        theta0_target   = 0;                    %True Anomaly (rad)
        a_target        = 35000e3;              %Semi-major Axis (m)
        e_target        = 0.8;               %Eccentricity
        i_target        = 59*deg2rads;            %Inclination (rad)
        RAAN_target     = 84*deg2rads;            %Right-ascension of the Ascending Node (rad)
        AoP_target      = 20*deg2rads;           %Argument of Perigee (rad)
        n_target        = sqrt(mu/a_target^3);  %Mean orbital motion of the target spacecraft
        T_target        = 2*pi/n_target;        %Orbital Period [s]
        
       COE_target = [theta0_target a_target e_target ...
                        i_target RAAN_target AoP_target];        
        
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        theta0_chaser   = theta0_target;        %True Anomaly (rad)
        a_chaser        = 35000.0126e3;             %Semi-major Axis (m)
        e_chaser        = 0.800010039;   %Eccentricity + delta e
        i_chaser        = 59.00305*deg2rads;             %Inclination (rad)
        RAAN_chaser     = RAAN_target;          %RAAN (rad)
        AoP_chaser      = AoP_target;           %Argument of Perigee (rad)
        
        COE_chaser = [theta0_chaser a_chaser e_chaser ...
                        i_chaser RAAN_chaser AoP_chaser];
        
        %Additional parameters for initial conditions
        p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot
        
        %%Initializing parameters for EKF propagation
        
        % Known Initial Conditions for In-Plane Elliptical Formation
        delta_e         = e_chaser-e_target;    %Difference in eccentricity
        x0_chaser       = -348.846286613834;    %Initial x-position of chaser(m)        
        y0_chaser       = -0.00318741255775024; %Initial y-position of chaser(m)
        z0_chaser       = 127.439860611236;     %Initial z-position of chaser(m)
        theta0          = theta0_target;        %Initial true anomaly of target (rad)
        rt0             = rt0_target;           %Initial radius of target (m)
        
        x_dot0_chaser   = -9.22026039551493e-06;%Initial x-velocity of chaser(m/s)
        y_dot0_chaser   = 0.785032789127976;    %Initial y-velocity of chaser(m)
        z_dot0_chaser   = 0.785032789127976;    %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;    %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;       %Initial change in radius (m/s)
        fprintf('Formation Configuration : %s \n', FF_Config) 

        % Target (Tango) Spacecraft Properties for Orbit Propagator
        PropOptions_target.M_Sat   = 42.5;    %Spacecraft mass
        PropOptions_target.A_Drag  = 0.38;    %Area exposed to drag forces 
        PropOptions_target.Cd      = 2.25;    %Coefficient of Drag
        PropOptions_target.A_Solar = 0.55;    %Area exposed to solar radiation
        PropOptions_target.Cr      = 1.2;     %Coefficient of Reflection

        % Chaser (Mango) Spacecraft Properties for Orbit Propagator
        PropOptions_chaser.M_Sat   = 154.4;   %Spacecraft mass
        PropOptions_chaser.A_Drag  = 1.3;     %Area exposed to drag forces
        PropOptions_chaser.Cd      = 2.5;     %Coefficient of Drag
        PropOptions_chaser.A_Solar = 2.5;     %Area exposed to solar radiation
        PropOptions_chaser.Cr      = 1.32;    %Coefficient of Reflection
        
    case 'LF_Hold'
        
        %Initial COE for Target Spacecraft
        theta0_target   = 0;                    %True Anomaly [rad]
        a_target        = 7500e3;               %Semi-major Axis [m]
        e_target        = 0.1;                  %Eccentricity
        i_target        = 98.1877*deg2rads;       %Inclination [rad]
        RAAN_target     = 189.8914*deg2rads;      %RAAN [rad]
        AoP_target      = 1.0938*deg2rads;        %Argument of Perigee [rad]
        n_target        = sqrt(mu/a_target^3);  %Mean orbital motion [rad/s]
        T_target        = 2*pi/n_target;        %Orbital Period [s]
        
        COE_target = [theta0_target a_target e_target ...
                        i_target RAAN_target AoP_target];
        
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        theta0_chaser   = theta0_target-1*deg2rads;        %True Anomaly [rad]
        a_chaser        = a_target;             %Semi-major Axis [m]
        e_chaser        = e_target;   %Eccentricity - delta e
        i_chaser        = i_target;             %Inclination [rad]
        RAAN_chaser     = RAAN_target;          %Right-ascension of the Ascending Node [rad]
        AoP_chaser      = AoP_target;           %Argument of Perigee [rad]
        n_chaser        = sqrt(mu/a_chaser^3);  %Mean orbital motion [rad/s]
        T_chaser        = 2*pi/n_chaser;        %Orbital Period [s]
        
        COE_chaser = [theta0_chaser a_chaser e_chaser ...
                        i_chaser RAAN_chaser AoP_chaser];
        
        %Additional parameters for initial conditions
        p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot
        
        %Initial Conditions for In-Plane Elliptical Formation in LVLH Frame
        x0_chaser       = -375;          %Initial x-position of chaser(m)
        y0_chaser       = 0;             %Initial y-position of chaser(m)
        z0_chaser       = 0;             %Initial z-position of chaser(m)
        theta0          = theta0_target; %Initial true anomaly of target (rad)
        rt0             = rt0_target;    %Initial radius of target (m)
        
        x_dot0_chaser   = y0_chaser*n_chaser/2; %Initial x-velocity of chaser(m/s)
        y_dot0_chaser   = -2*n_chaser*x0_chaser;%Initial y-velocity of chaser(m)
        z_dot0_chaser   = 0;                    %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;    %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;       %Initial change in radius (m/s)
        
        %Values extracted from RWOP (Needed for precise orbit propagation)
        %x0_chaser       = -374.999999999265;          %Initial x-position of chaser(m)
        %y0_chaser       = -1.27897692436818e-12;      %Initial y-position of chaser(m)
        %z0_chaser       = -9.41520195141266e-11;      %Initial z-position of chaser(m)
        %x_dot0_chaser   = -2.51534904016637e-14;      %Initial x-velocity of chaser(m/s)
        %y_dot0_chaser   = 0.854818112337925;          %Initial y-velocity of chaser(m)
        %z_dot0_chaser   = 1.35710886972618e-13;       %Initial z-velocity of chaser(m)
         fprintf('Formation Configuration : %s \n', FF_Config)   

end

%% GPS Measurement Specifications

%GPS Errors - Standard deviations for white gaussian noise
GPS_r_STD = [1.20 1.20 1.20]; %(m)  
GPS_v_STD = [0.03 0.03 0.03]; %(m/s)
                
%Initial measurements - Assuming measured positions/velocities are correct
x0_m            = x0_chaser;   
y0_m            = y0_chaser;
z0_m            = z0_chaser;
theta0_m        = theta0;
%rt0_m        = rt0;
x_dot0_m        = x_dot0_chaser;
y_dot0_m        = y_dot0_chaser;
z_dot0_m        = z_dot0_chaser;

%Initial measurement vector (7x1 Matrix)
meas0 =     [   x0_m 
                y0_m 
                z0_m 
                theta0_m
                %rt0_m
                x_dot0_m 
                y_dot0_m 
                z_dot0_m   ];

%% Initial EKF Parameters

%Initial state vector (10x1 Matrix)
state0 =    [   x0_chaser - 20             % x          [m]
                y0_chaser + 20             % y          [m]
                z0_chaser - 20             % z          [m]
                theta0                     % theta      [rad]
                rt0 + 50                   % rt         [m]
                x_dot0_chaser + 2/1000    % x_dot      [m/s]
                y_dot0_chaser - 2/1000    % y_dot      [m/s]
                z_dot0_chaser + 2/1000    % z_dot      [m/s]
                theta_dot0                 % theta_dot  [rad/s]
                rt_dot0    ];              % rt_dot     [m/s]

%Initial state error convariance (10x10 Matrix)
%P0 = eye(10)*10^1;     %Original Settings

sig2_r          = 100;                  %Position covariance [m^2]
sig2_theta      = ( sqrt(1)*(pi/180) )^2 ;    %True anomaly covariance [rad^2]
sig2_rt         = 10000;                 %Target radius covariance [m^2]
sig2_v          = 1;                    %Velocity covariance [m^2/s^2]
sig2_thetadot   = ( sqrt(1e-2)*(pi/180) )^2;   %Anomaly rate covariance [rad^2/s^2]
sig2_rtdot      = 100;                   %Target radius rate covariance [m^2/s^2]

P0              = [ sig2_r        0 0 0 0 0 0 0 0 0 
                    0 sig2_r        0 0 0 0 0 0 0 0 
                    0 0 sig2_r        0 0 0 0 0 0 0 
                    0 0 0 sig2_theta    0 0 0 0 0 0 
                    0 0 0 0 sig2_rt       0 0 0 0 0 
                    0 0 0 0 0 sig2_v        0 0 0 0 
                    0 0 0 0 0 0 sig2_v        0 0 0 
                    0 0 0 0 0 0 0 sig2_v        0 0 
                    0 0 0 0 0 0 0 0 sig2_thetadot 0 
                    0 0 0 0 0 0 0 0 0 sig2_rtdot    ];

% =========================================================================
%Initial Q Matrix (Covariance of the process noise, 10x10)

%Q Tuning Set 1 - Used in Papers
%{
    q2_r           = 2*10^-1;             %Position covariance [m^2]
    q2_theta       = 5*10^-3; %16.414 deg^2 %True anomaly covariance [rad^2]
    q2_rt          = 5*10^-1;             %Target radius covariance [m^2]
    q2_v           = 5*10^-3;             %Velocity covariance [m^2/s^2]
    q2_thetadot    = 5*10^-5; %0.164 deg^2  %Anomaly rate covariance [rad^2/s^2]
    q2_rtdot       = 5*10^-5;             %Target radius rate covariance [m^2/s^2]
%}
%Q Tuning Set 2
%
    q2_r           = 2*10^-1;             %Position covariance [m^2]
    q2_theta       = 3.046*10^-7;%1e-3 deg^2   %True anomaly covariance [rad^2]
    q2_rt          = 5*10^-3;             %Target radius covariance [m^2]
    q2_v           = 5*10^-3;             %Velocity covariance [m^2/s^2]
    q2_thetadot    = 3.046*10^-10;%0.1 deg %Anomaly rate covariance [rad^2/s^2]
    q2_rtdot       = 5*10^-5;             %Target radius rate covariance [m^2/s^2]
%}

    Q0              = [ q2_r        0 0 0 0 0 0 0 0 0 
                        0 q2_r        0 0 0 0 0 0 0 0 
                        0 0 q2_r        0 0 0 0 0 0 0 
                        0 0 0 q2_theta    0 0 0 0 0 0 
                        0 0 0 0 q2_rt       0 0 0 0 0 
                        0 0 0 0 0 q2_v        0 0 0 0 
                        0 0 0 0 0 0 q2_v        0 0 0 
                        0 0 0 0 0 0 0 q2_v        0 0 
                        0 0 0 0 0 0 0 0 q2_thetadot 0 
                        0 0 0 0 0 0 0 0 0 q2_rtdot    ];
                    
% =========================================================================                    
%Initial R Matrix (Covariance of the measurement noise, 7x7)

%R-Tuning Set 1
    r2_rm           = 2e1;                  %Position covariance [m^2]
    r2_thetam       = (sqrt(1)*(pi/180) )^2;   %1e-2 deg^2;             
    r2_vm           = 5e-1;                 

%R-Tuning Set 2
    %r2_rm           = 2e1; 
    %r2_thetam       = 1*10^-1;    
    %r2_vm           = 5e-2;
    
    R0              = [ r2_rm     0 0 0 0 0 0 %0
                        0 r2_rm     0 0 0 0 0 %0
                        0 0 r2_rm     0 0 0 0 %0
                        0 0 0 r2_thetam 0 0 0 %0
                        0 0 0 0 r2_vm     0 0 %0
                        0 0 0 0 0 r2_vm     0 %0  
                        0 0 0 0 0 0 r2_vm     ];   %0 ];  
 
% Covariances when units were km                   
%{                   
% =========================================================================
%Initial Q Matrix (Covariance of the process noise, 10x10)
    sig2_r           = 2*10^-6;             %Variance of the position noise
    sig2_theta       = 5*10^-9;
    sig2_rt          = 5*10^-10; 
    sig2_v           = 5*10^-9; 
    sig2_thetadot    = 5*10^-15; 
    sig2_rtdot       = 5*10^-10; 

    Q0              = [ sig2_r        0 0 0 0 0 0 0 0 0 
                        0 sig2_r        0 0 0 0 0 0 0 0 
                        0 0 sig2_r        0 0 0 0 0 0 0 
                        0 0 0 sig2_theta    0 0 0 0 0 0 
                        0 0 0 0 sig2_rt       0 0 0 0 0 
                        0 0 0 0 0 sig2_v        0 0 0 0 
                        0 0 0 0 0 0 sig2_v        0 0 0 
                        0 0 0 0 0 0 0 sig2_v        0 0 
                        0 0 0 0 0 0 0 0 sig2_thetadot 0 
                        0 0 0 0 0 0 0 0 0 sig2_rtdot    ];
                    
                   %Q0 = 10^-2*eye(10,10);
                    
% =========================================================================                    
%Initial R Matrix (Covariance of the measurement noise, 7x7)
    sig2_rm           = 2*10^-4; %5*10^-3; 
    sig2_thetam       = 1*10^-6;
    sig2_vm           = 5*10^-6; %1.5*10^-4; 
    
    R0              = [ sig2_rm     0 0 0 0 0 0 
                        0 sig2_rm     0 0 0 0 0 
                        0 0 sig2_rm     0 0 0 0 
                        0 0 0 sig2_thetam 0 0 0 
                        0 0 0 0 sig2_vm     0 0 
                        0 0 0 0 0 sig2_vm     0   
                        0 0 0 0 0 0 sig2_vm      ];  
                    
                   % R0 = 10^-2*eye(7,7);
                      
% =========================================================================
%}
%========================================== 

%% Assembling Initial EKF Vector 

%Arranging the covariance matrix into a vector (1 x 100)
P0            = [ P0(1,:) P0(2,:) P0(3,:) P0(4,:) P0(5,:) ...
                    P0(6,:) P0(7,:) P0(8,:) P0(9,:) P0(10,:)]; 

% Initial EKF Vector (1x118)                 
EKF0            = [ t0 state0' P0 meas0' ];

% Naming and saving parameter file
clear Constants ConstantsModel

filename        = ['FFNAV_',FF_Config,'_Params.mat'];
file_out = fullfile(paths.output, filename);
save(file_out)
fprintf('Parameter file created  : %s \n', filename)
% =========================================================================