%% InitializeRWSW =====================================================
% Description: This script defines and loads all data needed for the real-
% world software simulators.
%
% Inputs:
%   Parameter_filename  - Spacecraft and formation parameters filename
%   PropOptions.        - Structure for orbit propagator settings
%   time_start          - Start time for the simulation
%
% Outputs:
%   Results.mat file
% Other Functions Called:
%
% Created by:  Cory Fraser - NOV 18, 2018
% Latest Edit: Cory Fraser - NOV 18, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Physical Constants
ConstantsModel = 'WGS84 Model';

%fprintf('\n--------------------------------------------------------------\n')
fprintf('\n Loading Real-World Software Parameters...')

%Gravitational Parameters (mu, GM, gm)
Constants.mu_Earth   = 3.98600435436e14;        % [m^3/s^2]; WGS84
Constants.mu_Sun     = 132712440041.939377e9; % [m^3/s^2]; DE436
Constants.mu_Moon    = 4902.800117e9;         % [m^3/s^2]; DE436
Constants.mu_Mercury = 22031.780000e9;        % [m^3/s^2]; DE436
Constants.mu_Venus   = 324858.592000e9;       % [m^3/s^2]; DE436
Constants.mu_Mars    = 42828.375214e9; 	  	 % [m^3/s^2]; DE436
Constants.mu_Jupiter = 126712764.133446e9;    % [m^3/s^2]; DE436
Constants.mu_Saturn  = 37940585.200000e9;     % [m^3/s^2]; DE436
Constants.mu_Uranus  = 5794556.465752e9;      % [m^3/s^2]; DE436
Constants.mu_Neptune = 6836527.100580e9;      % [m^3/s^2]; DE436
Constants.mu_Pluto   = 975.501176e9;   		 % [m^3/s^2]; DE436

% Earth Parameters
Constants.J2        = 1082.64e-6;        %Earth's J2 Coefficient
Constants.R_Earth   = 6378.137e3;        %Equatorial Radius - WGS84 [m]
Constants.R_Earth_p = 6356.752e3;        %Polar Radius - WGS84 [m]
Constants.w_Earth   = 7.2921151467e-5;   %Mean angular velocity - WGS84 [rad/sec]
Constants.f_Earth   = 1/298.257223563;   %Flattening; WGS-84

% Other Data
Constants.R_Sun     = 696000e3;            % Sun's radius [m]; DE436
Constants.R_Moon    = 1738e3;              % Moon's radius [m]; DE436
Constants.c_light   = 299792457.999999984; % Speed of light  [m/s]; DE436
Constants.AU        = 149597870699.999988; % Astronomical unit [m]; DE436

Constants.MJD_J2000 = 51544.5;             % Modified Julian Date of J2000
Constants.T_B1950   = -0.500002108;        % Epoch B1950

% Solar radiation pressure at 1 AU
Constants.P_Solar = 1367/Constants.c_light; % [N/m^2] (1367 W/m^2); IERS

%Conversion Factors
deg2rad          = pi/180;               % Degrees to radians
rad2deg          = 180/pi;               % Radians to degrees 
mu               = Constants.mu_Earth;

%Planetary Ephemeride Coefficients for DE436 
%   load DE436Coeff.mat
%   PC = DE436Coeff;

SphericalHarmonicCoefficients
Constants.CS = CS_JGM3;

%Caclulate ECI position of Sun and Moon (assuming fixed)
Constants.Rvec_Sun  = PositionSun(JulianDate0, Constants.AU );
Constants.Rvec_Moon = PositionMoon(JulianDate0, Constants.R_Earth);

%==========================================================================
%% Select Formation Configuration
FF_Config = 'PRISMA';
% PRISMA     = PRISMA Mission (LEO, low eccentricty)
% PROBA-3    = PROBA-3 Mission (LEO, high eccentricity)
% PEOinLEO   = Modified PRISMA, with inclination change (LEO,low e,delta i)
    
switch FF_Config
    
    case 'PRISMA'
        
        %Initial COE for Target Spacecraft
        COE_target.TA   = 358.90349028*deg2rad; %True Anomaly [rad]
        COE_target.SMA  = 7087.29755686634e3;   %Semi-major Axis [m]
        COE_target.ECC  = 0.00145443;           %Eccentricity
        COE_target.INC  = 98.18528613*deg2rad;  %Inclination [rad]
        COE_target.RAAN = 189.8914*deg2rad;     %RAAN [rad]
        COE_target.AOP  = 1.097451382*deg2rad;  %Argument of Perigee [rad]
        COE_target.n    = sqrt(mu/COE_target.SMA^3);  %Mean orbital motion [rad/s]
        COE_target.T    = 2*pi/COE_target.n;        %Orbital Period [s]
        
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        COE_chaser.TA   = 0*deg2rad;         %True Anomaly [rad]
        COE_chaser.SMA        = 7087.29767733179e3;   %Semi-major Axis [m]
        COE_chaser.ECC        = 0.00145908;           %Eccentricity - delta e
        COE_chaser.INC        = 98.18466676*deg2rad;  %Inclination [rad]
        COE_chaser.RAAN     = 189.8908602*deg2rad;  %Right-ascension of the Ascending Node [rad]
        COE_chaser.AOP      = 0.00*deg2rad;            %Argument of Perigee [rad]
        COE_chaser.n        = sqrt(mu/COE_chaser.SMA^3);  %Mean orbital motion [rad/s]
        COE_chaser.T        = 2*pi/COE_chaser.n;        %Orbital Period [s]
                
        %Additional parameters for initial conditions
        p_target        = COE_target.SMA*(1-COE_target.ECC^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+COE_target.ECC*cos(COE_target.TA));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*COE_target.ECC*sin(COE_target.TA);    %Initial r_dot
        
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

        x0_chaser       = -34.718875426019130;    %Initial x-position of chaser(m)
        y0_chaser       = -1.068178082293263e+02; %Initial y-position of chaser(m)
        z0_chaser       = 65.994822637671480;     %Initial z-position of chaser(m)
        theta0          = COE_target.TA;          %Initial true anomaly of target (rad)
        rt0             = rt0_target;             %Initial radius of target (m)
        
        x_dot0_chaser   = 0.208732017363154;  %Initial x-velocity of chaser(m/s)
        y_dot0_chaser   = 0.073700404718392;  %Initial y-velocity of chaser(m)
        z_dot0_chaser   = 0.081187656642932;  %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;  %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;     %Initial change in radius (m/s)
    
    case 'PEOinLEO'
        
        %Initial COE for Target Spacecraft
        COE_target.TA   = 0;                    %True Anomaly [rad]
        COE_target.SMA        = 7500e3;               %Semi-major Axis [m]
        COE_target.ECC        = 0.1;                  %Eccentricity
        COE_target.INC        = 98.1877*deg2rad;       %Inclination [rad]
        COE_target.RAAN     = 189.8914*deg2rad;      %RAAN [rad]
        COE_target.AOP      = 1.0938*deg2rad;        %Argument of Perigee [rad]
        COE_target.n        = sqrt(mu/COE_target.SMA^3);  %Mean orbital motion [rad/s]
        COE_target.T        = 2*pi/COE_target.n;        %Orbital Period [s]
            
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        COE_chaser.TA   = COE_target.TA;        %True Anomaly [rad]
        COE_chaser.SMA        = COE_target.SMA;             %Semi-major Axis [m]
        COE_chaser.ECC        = COE_target.ECC + 0.00005;   %Eccentricity - delta e
        COE_chaser.INC        = COE_target.INC-0.01*deg2rad;   %Inclination [rad]
        COE_chaser.RAAN     = COE_target.RAAN;          %Right-ascension of the Ascending Node [rad]
        COE_chaser.AOP      = COE_target.AOP;           %Argument of Perigee [rad]
        COE_chaser.n        = sqrt(mu/COE_chaser.SMA^3);  %Mean orbital motion [rad/s]
        COE_chaser.T        = 2*pi/COE_chaser.n;        %Orbital Period [s]
        
        %Additional parameters for initial conditions
        p_target        = COE_target.SMA*(1-COE_target.ECC^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+COE_target.ECC*cos(COE_target.TA));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*COE_target.ECC*sin(COE_target.TA);    %Initial r_dot
        
        %Initial Conditions for In-Plane Elliptical Formation in LVLH Frame
        x0_chaser       = -3.750000374605406e+02;          %Initial x-position of chaser(m)
        y0_chaser       =   -0.00196;             %Initial y-position of chaser(m)
        z0_chaser       =  -22.48775;             %Initial z-position of chaser(m)
        theta0          = COE_target.TA; %Initial true anomaly of target (rad)
        rt0             = rt0_target;    %Initial radius of target (m)
        x_dot0_chaser   = -4.68572441684328e-06;%Initial x-velocity of chaser(m/s)
        y_dot0_chaser   = 0.854695440915271;    %Initial y-velocity of chaser(m)
        z_dot0_chaser   = -1.40647980533729;    %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;    %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;       %Initial change in radius (m/s)
        
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
     
    case 'PROBA-3'
    
        %Definitions for RWOP - Initial COE for Target Spacecraft
        COE_target.TA   = 0;                    %True Anomaly (rad)
        COE_target.SMA        = 36943e3;              %Semi-major Axis (m)
        COE_target.ECC        = 0.8111;               %Eccentricity
        COE_target.INC        = 59*deg2rad;            %Inclination (rad)
        COE_target.RAAN     = 84*deg2rad;            %Right-ascension of the Ascending Node (rad)
        COE_target.AOP      = 188*deg2rad;           %Argument of Perigee (rad)
        COE_target.n        = sqrt(mu/COE_target.SMA^3);  %Mean orbital motion of the target spacecraft
        COE_target.T        = 2*pi/COE_target.n;        %Orbital Period [s]
             
        %Definitions for RWOP - Initial COE for Chaser Spacecraft
        COE_chaser.TA   = COE_target.TA;            %True Anomaly (rad)
        COE_chaser.SMA        = COE_target.SMA;                 %Semi-major Axis (m)
        COE_chaser.ECC        = COE_target.ECC + 0.000005;      %Eccentricity + delta e
        COE_chaser.INC        = COE_target.INC;                 %Inclination (rad)
        COE_chaser.RAAN     = COE_target.RAAN;              %RAAN (rad)
        COE_chaser.AOP      = COE_target.AOP;               %Argument of Perigee (rad)
        
        %Additional parameters for initial conditions
        p_target        = COE_target.SMA*(1-COE_target.ECC^2);                          %Semi-latus rectum
        h_target        = sqrt(mu*p_target);                                %Angular momentum
        rt0_target      = p_target/(1+COE_target.ECC*cos(COE_target.TA));         %Initial radius
        theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot
        rt_dot0_target  = sqrt(mu/p_target)*COE_target.ECC*sin(COE_target.TA);    %Initial r_dot
                
        % Known Initial Conditions for In-Plane Elliptical Formation
        x0_chaser       = -184.715000001348;     %Initial x-position of chaser(m)
        y0_chaser       = 1.29855237673837e-10;  %Initial y-position of chaser(m)
        z0_chaser       = -9.72395497456091e-11; %Initial z-position of chaser(m)
        theta0          = COE_target.TA;         %Initial true anomaly of target (rad)
        rt0             = rt0_target;            %Initial radius of target (m)

        x_dot0_chaser   = 1.79981030079546e-13;  %Initial x-velocity of chaser(m/s)
        y_dot0_chaser   = 0.417862035614425;     %Initial y-velocity of chaser(m)
        z_dot0_chaser   = -2.43138842392909e-14; %Initial z-velocity of chaser(m)
        theta_dot0      = theta_dot0_target;     %Initial true anomaly rate (rad/s)
        rt_dot0         = rt_dot0_target;        %Initial change in radius (m/s)
        
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
end

%Initial Position/Velcoity vectors for the spacecraft
[R0_target, V0_target] = RVfromCOE(mu, COE_target);
[R0_chaser, V0_chaser] = RVfromCOE(mu, COE_chaser);

%==========================================================================
%% Select Orbit Propagator Options (1 = ON, 0 = OFF)
PropOptions.J2              = 1*0;
PropOptions.GravitySpherical= 0;
PropOptions.Sun             = 1*0;        %Assuming fixed position of Sun
PropOptions.Moon            = 1*0;        %Assuming fixed position of Moon
PropOptions.Planets         = 0;       %Assuming fixed positions of planets
PropOptions.Drag            = 1*0;        %Using Harris-Priester Model  
PropOptions.SolarRad        = 1*0;        
    PropOptions.ShadowModel = 0;        %0=Cylindrical Shadow, 1=Geometric 
PropOptions.Relativity      = 0;        
PropOptions.GravEl_Earth    = 0;        %NOT COMPLETE
    
% Create PropOptions for target spacecraft (same for both spacecraft)
PropOptions_target.J2              = PropOptions.J2;
PropOptions_target.GravitySpherical= PropOptions.GravitySpherical;
PropOptions_target.Sun             = PropOptions.Sun;
PropOptions_target.Moon            = PropOptions.Moon;
PropOptions_target.Planets         = PropOptions.Planets;
PropOptions_target.Drag            = PropOptions.Drag;        
PropOptions_target.SolarRad        = PropOptions.SolarRad;        
    PropOptions_target.ShadowModel = PropOptions.ShadowModel; 
PropOptions_target.Relativity      = PropOptions.Relativity;        
PropOptions_target.GravEl_Earth    = PropOptions.GravEl_Earth;

% Create PropOptions for chaser spacecraft (same for both spacecraft)
PropOptions_chaser.J2              = PropOptions.J2;
PropOptions_chaser.GravitySpherical=PropOptions.GravitySpherical;
PropOptions_chaser.Sun             = PropOptions.Sun;
PropOptions_chaser.Moon            = PropOptions.Moon;
PropOptions_chaser.Planets         = PropOptions.Planets;
PropOptions_chaser.Drag            = PropOptions.Drag;        
PropOptions_chaser.SolarRad        = PropOptions.SolarRad;        
    PropOptions_chaser.ShadowModel = PropOptions.ShadowModel; 
PropOptions_chaser.Relativity      = PropOptions.Relativity;        
PropOptions_chaser.GravEl_Earth    = PropOptions.GravEl_Earth;

%==========================================================================
%% Finish

fprintf(' Complete \n')
fprintf(' \t Formation Configuration : %s \n', FF_Config)   
fprintf(' \t Astrodynamic constants  : %s \n', ConstantsModel)
fprintf(' \t Orbital Perturbations : ')
if (PropOptions.J2 || PropOptions.GravEl_Earth || PropOptions.Sun || ...
        PropOptions.Moon ||  PropOptions.Planets || PropOptions.Relativity...
        || PropOptions.Drag || PropOptions.SolarRad )
    fprintf('On \n')
else
    fprintf('None (Two-Body Motion) \n')
end
if (PropOptions.J2)
    fprintf('     - Earth''s Oblateness (J2) \n')
end
if (PropOptions.Sun)
    fprintf('     - Solar Third-Body Gravity\n')
end
if (PropOptions.Moon)
    fprintf('     - Lunar Third-Body Gravity \n')
end
if (PropOptions.Planets)
    fprintf('     - Planetary Third-Body Gravity \n')
end    
if (PropOptions.SolarRad) 
    fprintf('     - Solar Radiation Pressure \n')    
end
if (PropOptions.Drag)
    fprintf('     - Atmospheric Drag \n')
end
if (PropOptions.Relativity)
    fprintf('     - Relativistic Effects \n')
end
if (PropOptions.GravEl_Earth)
    fprintf('     - Harmonic Gravity (Elastic Earth) \n')
end

% =========================================================================