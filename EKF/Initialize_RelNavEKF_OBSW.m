%% InitializeOBSW =====================================================
% Description: This script defines and loads all data needed for the space-
% craft on-board software.
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

%% Select Kalman Filter Options

%EKF Propagation Time
EKF_dt = 1;

%EKF Measurement Sampling Time
EKF_r_rel_meas_dt = EKF_dt;
EKF_v_rel_meas_dt = EKF_dt;
EKF_theta_meas_dt = EKF_dt;

% 0 = EKF           Extended Kalman Filter        
% 1 = MLE-EKF       MLE Adaptive Extended Kalman Filter
% 2 = FAEKF         Fuzzy Adaptive Extended Kalman Filter 
EKF_flag = 0;

% Select Adaptation Options - MLE-AEKF & FAEKF (1 = ON, 0 = OFF)
Q_adapt = 1; 
R_adapt = 1; 

% Select Residuals Smoothing Options - MLE-AEKF & FAEKF (1 = ON, 0 = OFF)
Smooth_flag = '1';
N_window    = 30;  %Running-Average Window Size 
%N=5 seems to work better for MLE-AEKF


%fprintf('\n--------------------------------------------------------------\n')
fprintf('\n Loading On-board Software Parameters...')

%% Initial EKF Parameters

%Initial state vector (10x1 Matrix)
state0 =    [   x0_chaser - 20             % x          [m]
                y0_chaser + 20             % y          [m]
                z0_chaser - 20             % z          [m]
                theta0                     % theta      [rad]
                rt0 + 50                   % rt         [m]
                x_dot0_chaser + 2/1000     % x_dot      [m/s]
                y_dot0_chaser - 2/1000     % y_dot      [m/s]
                z_dot0_chaser + 2/1000     % z_dot      [m/s]
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

%Q Tuning Set 2
    q2_r           = 2*10^-1;             %Position covariance [m^2]
    q2_theta       = 3.046*10^-7;%1e-3 deg^2   %True anomaly covariance [rad^2]
    q2_rt          = 5*10^-3;             %Target radius covariance [m^2]
    q2_v           = 5*10^-3;             %Velocity covariance [m^2/s^2]
    q2_thetadot    = 3.046*10^-10;%0.1 deg %Anomaly rate covariance [rad^2/s^2]
    q2_rtdot       = 5*10^-5;             %Target radius rate covariance [m^2/s^2]

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
    
    R0              = [ r2_rm     0 0 0 0 0 0 
                        0 r2_rm     0 0 0 0 0 
                        0 0 r2_rm     0 0 0 0 
                        0 0 0 r2_thetam 0 0 0 
                        0 0 0 0 r2_vm     0 0 
                        0 0 0 0 0 r2_vm     0   
                        0 0 0 0 0 0 r2_vm     ];  

    % With rt measured (8 measurements)
    %{
    R0              = [ r2_rm     0 0 0 0 0 0 0
                        0 r2_rm     0 0 0 0 0 0
                        0 0 r2_rm     0 0 0 0 0
                        0 0 0 r2_thetam 0 0 0 0
                        0 0 0 0 10 0 0 0% rt extra
                        0 0 0 0 0 r2_vm     0 0 
                        0 0 0 0 0 0 r2_vm     0   
                        0 0 0 0 0 0 0 r2_vm     ];  
 %}
 
fprintf(' Complete \n')

% =========================================================================
switch EKF_flag
    
    %Extended Kalman Filter
    case 0
        EKF_name = 'EKF';
        fprintf(' \t Navigation Filtering Algorithm   : EKF \n')
        File_specifier = 'NOadapt_5orbits.mat';

    %MLE Adaptive Extended Kalman Filter
    case 1 
        EKF_name = 'MLE_EKF';
        fprintf(' \t Navigation Filtering Algorithm   : MLE Adaptive EKF \n')   
        fprintf('RTS Smoother          : On, N = %i \n', N_window)

        %Set Q and R Adaptation Switches
        if (Q_adapt) 
            fprintf('Q-adaptations         : On \n')
        else
            fprintf('Q-adaptations         : Off \n')
        end
        if (R_adapt)
            fprintf('R-adaptations         : On \n')
        else
            fprintf('R-adaptations         : Off \n')
        end
             
    %Fuzzy Adaptive Extended Kalman Filter
    case 2 
        EKF_name = 'FAEKF';
        fprintf(' \t Navigation Filtering Algorithm   : Fuzzy Adaptive EKF \n')
        
        % Initialize Fuzzy Inference System
        fprintf('Fuzzy system loaded   :')
        MSFplot_flag = 0; % 1 = Plot Membership Functions
        FUZZY_props = FFNAV_FAEKF_FuzzyInitialize(MSFplot_flag);
        
        %Set Smoothing Settings
        if Smooth_flag == '1'
            fprintf('Residual smoothing    : On, N = %i \n', N_window)
        else
            fprintf('Residual smoothing    : Off \n')
        end
        set_param('EKF_Sim_FAEKF/Kalman Filter/FAEKF/Adaptations/Switch: Smoothing On//Off', 'sw', Smooth_flag)
        
        %Set Q and R Adaptation Switches
        if (Q_adapt)
            fprintf('Q-adaptations         : On \n')
        else
            fprintf('Q-adaptations         : Off \n')
        end
        if (R_adapt)
            fprintf('R-adaptations         : On \n')
        else
            fprintf('R-adaptations         : Off \n')
        end
end

%Set output filename based on adaptation scheme
if (EKF_flag~=0) && (Q_adapt) && (R_adapt)
    File_specifier = 'QRadapt.mat';
elseif (EKF_flag~=0) && (Q_adapt)
    File_specifier = 'Qadapt.mat';
elseif (EKF_flag~=0) && (R_adapt)
    File_specifier = 'Radapt.mat';
else
    File_specifier = 'NOadapt.mat';
end

% =========================================================================