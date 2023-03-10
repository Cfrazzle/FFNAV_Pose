function [EKF_out, Storage_out, L, Q_update, R_update] = ...
    FFNAV_EKF_MLE(Storage_in, time_step, mu, data_pre, sensors, ...
    Q, R, Q_adapt_flag,R_adapt_flag)
% FFNAV Maximum Likelikelihood Estimation EKF =============================
% Description: This function completes both the propagation and correction
% phases of the EKF, along with adapting the process and measurement noise
% covariances Q and R using an MLE method.
%
% Inputs:
%   Storage_in   - Matrix of past data needed for the RTS Smoother
%   time_step    - Time step of the simulation
%   mu           - Earth's gravitational parameter
%   data_pre     - EKF state and covariance data from previous time step
%   sensors      - Vector of measurements
%   Q            - Process noise covariance matrix
%   Q_adapt_flag - Flag to adapt Q (0 = no adaptation) 
%   R            - Measurement noise covariance matrix
%   R_adapt_flag - Flag to adapt R (1 = adapt) 
%
% Outputs:
%   EKF_out     - Vector of EKF output data, to be passed to next EKF step
%   Storage_out - Vector of data to be stored in memory for RTS Smoother
%   L           - Value of the pseudo-likelihood cost function
%   Q_update    - Updated process noise covariance matrix
%   R_update    - Updated measurement noise covariance matrix
%
% Other Functions Called:
%   FFNAV_EKFdot  - 
%   FFNAV_EKF_STM -
%   findPSDS - Finds nearest positive semidefinite symmetric (PSDS) matrix
%   chol     - Performs a Cholesky decomposition on a PDS matrix
%
% References:
%   Bavedekar
%   _________
%
% Notes:
%   1) Data_pre comes in as a single column vector
%   2) Storage_in comes in as a matrix, with each row corresponding to the
%           full set of output data from the previous filtering iterations
%
% Created by:  Cory Fraser - AUG 22, 2017
% Latest Edit: Cory Fraser - JUL 04, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Initialize Parameters
    
T = time_step;              %Time step for integration
coder.extrinsic('legend');
coder.extrinsic('fprintf');
coder.extrinsic('chol');

% Previous State Estimates
x         = data_pre(1);
y         = data_pre(2);
z         = data_pre(3);
theta     = data_pre(4);
rt        = data_pre(5);
x_dot     = data_pre(6);
y_dot     = data_pre(7);
z_dot     = data_pre(8);
theta_dot = data_pre(9);
rt_dot    = data_pre(10);

state_pre = zeros(10,1);
state_pre =  [   x; y; z; theta; rt; ...
    x_dot; y_dot; z_dot; ...
    theta_dot; rt_dot ]; %(10 x 1)

% Unpropagated State Error Covariance (10x10)
P_pre = [ data_pre(11:20)'
    data_pre(21:30)'
    data_pre(31:40)'
    data_pre(41:50)'
    data_pre(51:60)'
    data_pre(61:70)'
    data_pre(71:80)'
    data_pre(81:90)'
    data_pre(91:100)'
    data_pre(101:110)' ];

% =========================================================================
%% EKF Step 2: Propagation Phase
% Propagation is completed using a discrete-time linear model of the system

%Dynamics Propagation using RKM Method
k1 = FFNAV_EKFdot(state_pre, mu);
k2 = FFNAV_EKFdot(state_pre + (1/3)*k1*T, mu);
k3 = FFNAV_EKFdot(state_pre + (1/6)*k1*T + (1/6)*k2*T, mu);
k4 = FFNAV_EKFdot(state_pre + (1/8)*k1*T + (3/8)*k3*T, mu);
k5 = FFNAV_EKFdot(state_pre + (1/2)*k1*T - (3/2)*k3*T + 2*k4*T, mu);

state_priori = state_pre + (1/6)*(k1 + 4*k4 + k5)*T; %(10 x 1)

%Calculating the a priori state error covariance (10 x 10)
Phi = FFNAV_EKF_STM(state_pre, T, mu);
P_priori = Phi*P_pre*Phi' + Q;

% =========================================================================
%% EKF Step 3: Correction Phase

% Measured position and velocity parameters
x_m             = sensors(1);
y_m             = sensors(2);
z_m             = sensors(3);
theta_m         = sensors(4);
x_dot_m         = sensors(5);
y_dot_m         = sensors(6);
z_dot_m         = sensors(7);

%Assembling the measured outputs (8 x 1 matrix)
Z_m = [     x_m
            y_m
            z_m
            theta_m
            x_dot_m
            y_dot_m
            z_dot_m     ];

%Defining the measurement model (7 x 10 matrix)
H = [   1 0 0 0 0 0 0 0 0 0
        0 1 0 0 0 0 0 0 0 0
        0 0 1 0 0 0 0 0 0 0
        0 0 0 1 0 0 0 0 0 0 %Providing theta
        %0 0 0 0 1 0 0 0 0 0 %Providing rt
        0 0 0 0 0 1 0 0 0 0
        0 0 0 0 0 0 1 0 0 0
        0 0 0 0 0 0 0 1 0 0];

%Calulating the Estimated Measurements (7x1 matrix)
Z_est = H*state_priori;

%Calculating the residuals (7x1 matrix)
residuals   = (Z_m - Z_est);

%Calculating the theoretical Residual Covariance (7x7 Matrix)
Pr_theo = H*P_priori*H' + R;

%Pr_inv  = inv(Pr_theo);        %This inversion gives numerical issues
Pr_inv = Pr_theo\eye(7);       %Not sufficient for MLE-AEKF

%Inverting the innovations covariance using Cholesky Decomposition
%{
Pr_inv = zeros(7);
PSDS_flag = 0;    %Flag for Positive Semi-Definite Symmetric
[M_upper, PSDS_flag] = chol(Pr_theo);

if PSDS_flag~=0
    fprintf('Pr_theo is not PDS - finding nearest PSDS Matrix \n')
    [M_upper, PSDS_flag] = chol(findPSDS(Pr_theo));
end
M_inv = M_upper \ eye(7);
Pr_inv = M_inv * M_inv';
%}

%Calculating the Kalman Gain (10x7 matrix)
K = (P_priori*H')*Pr_inv;

%Correcting the State Estimate (10x1 matrix)
correction  = K*residuals;
state_post  = state_priori + correction;

%Correct the state error covariance matrix, using Joseph's Form  (10 x 10 matrix)
P_post = (eye(10)-K*H)*P_priori*(eye(10)-K*H)' + K*R*K';

%Calculate Negative-Log Likelihood Psuedo-cost function
%L = log(abs(det(Pr_theo))) + residuals'*inv(Pr_theo)*residuals;
L = log(abs(det(Pr_theo))) + residuals'*Pr_inv*residuals;

Q_update = Q;
R_update = R;
% ========================================================================

%% Maximum Likelihood Estimation - Adaptation Algorithm   
if (Q_adapt_flag) || (R_adapt_flag)
    
    N     = size(Storage_in,1)+1;   %Size based on fixed window size
    N_use = nnz(Storage_in(:,1))+1; %Useful size (# of non-zero rows)
    
    if N_use > 2
        
        %STEP 1 - Organizing the stored data
        state_post_mem      = zeros(10,1,N);
        state_priori_mem    = zeros(10,1,N);
        P_post_mem          = zeros(10,10,N);
        P_priori_mem        = zeros(10,10,N);
        Phi_mem             = zeros(10,10,N);
        Zm_mem              = zeros(7,1,N);
        Zest_mem            = zeros(7,1,N);
        
        state_post_mem(:,1,N_use)    = state_post;
        state_priori_mem(:,1,N_use)  = state_priori;
        Phi_mem(:,:,N_use)           = Phi;
        P_post_mem(:,:,N_use)        = P_post;
        P_priori_mem(:,:,N_use)      = P_priori;
        Zm_mem(:,1,N_use)            = Z_m;
        Zest_mem(:,1,N_use)          = Z_est;
        
        % Arranging data - Flipping, so N(ewest) is at last index
        for j = 1:N_use-1
            state_post_mem(:,1,N_use-j)      = [Storage_in(j,1:10)]';
            state_priori_mem(:,1,N_use-j)    = [Storage_in(j,11:20)]';
            
            P_post_mem(:,:,N_use-j) =...
              [ Storage_in(j,21:30) %Note that the rows in Storage_in
                Storage_in(j,31:40) %correspond to each filter step
                Storage_in(j,41:50) % (different than output vector)
                Storage_in(j,51:60)
                Storage_in(j,61:70)
                Storage_in(j,71:80)
                Storage_in(j,81:90)
                Storage_in(j,91:100)
                Storage_in(j,101:110)
                Storage_in(j,111:120) ];
            
            P_priori_mem(:,:,N_use-j) =...
              [ Storage_in(j,121:130)
                Storage_in(j,131:140)
                Storage_in(j,141:150)
                Storage_in(j,151:160)
                Storage_in(j,161:170)
                Storage_in(j,171:180)
                Storage_in(j,181:190)
                Storage_in(j,191:200)
                Storage_in(j,201:210)
                Storage_in(j,211:220) ];
            
            Phi_mem(:,:,N_use-j) = ...
              [ Storage_in(j,221:230)
                Storage_in(j,231:240)
                Storage_in(j,241:250)
                Storage_in(j,251:260)
                Storage_in(j,261:270)
                Storage_in(j,271:280)
                Storage_in(j,281:290)
                Storage_in(j,291:300)
                Storage_in(j,301:310)
                Storage_in(j,311:320) ];
            
            Zm_mem(:,1,N_use-j)   = [Storage_in(j,321:327)]';
            Zest_mem(:,1,N_use-j) = [Storage_in(j,328:334)]';
        end

    %STEP 2 - Extended Kalman Smoother (for k = N-1, N-2,..., 1)
    state_smooth    = zeros(10,1,N);
    P_smooth        = zeros(10,10,N);
    G_smooth        = zeros(10,10,N);
    Zest_smooth     = zeros(7,1,N);
    
    k = N_use;
    state_smooth(:,:,k) = state_post_mem(:,:,k);
    P_smooth(:,:,k)     = P_post_mem(:,:,k);
    Zest_smooth(:,:,k)  = H*state_smooth(:,:,k);
    
    for k = N_use:-1:2
        P_priori_mem_inv = P_priori_mem(:,:,k)\eye(10,10);        %Seems sufficient for EKF, FAEKF 
        G_smooth(:,:,k-1)     = P_post_mem(:,:,k-1)*Phi_mem(:,:,k-1)'*P_priori_mem_inv;
        %G_smooth(:,:,k-1)     = P_post_mem(:,:,k-1)*Phi_mem(:,:,k-1)'*pinv(P_priori_mem(:,:,k));
        state_smooth(:,:,k-1) = state_post_mem(:,:,k-1) + G_smooth(:,:,k-1)*(state_smooth(:,k) - state_priori_mem(:,:,k));
        Zest_smooth(:,:,k-1)  = H*state_smooth(:,:,k-1);
        P_smooth(:,:,k-1)     = P_post_mem(:,:,k-1) + G_smooth(:,:,k-1)*(P_smooth(:,:,k) - P_priori_mem(:,:,k))*G_smooth(:,:,k-1)';
    end
    state_smooth(:,:,1) = state_post_mem(:,:,1);
    Zest_smooth(:,:,1)  = H*state_smooth(:,:,1);
    P_smooth(:,:,1)     = P_post_mem(:,:,1);
    
    %===============================================
    %Plots for Debugging the Smoother
    %{
    if N_useful == 31

      TF_size       = 24;               %Title font size
        AF_size       = 20;               %Axes label font size
        AT_size       = 20;               %Tick label font size
    %close all
    figure(1)
    clf
    %title('EKF Output States')
        subplot(5,2,1);
        plot(squeeze(state_post_mem(1,1,:)))
        set(gca,'FontSize',AT_size);
        title('EKF Output States')
        ylabel('X(m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(state_smooth(1,1,:)))
        plot(squeeze(F_smooth(1,1,:)))

        subplot(5,2,3);
        plot(squeeze(state_post_mem(2,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Y(m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(state_smooth(2,1,:)))
        plot(squeeze(F_smooth(2,1,:)))

        subplot(5,2,5);
        plot(squeeze(state_post_mem(3,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Z (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(state_smooth(3,1,:)))
        plot(squeeze(F_smooth(3,1,:)))

        subplot(5,2,7);
        plot(squeeze(state_post_mem(4,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('\theta (rad)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(state_smooth(4,1,:)))
        plot(squeeze(F_smooth(4,1,:)))

        subplot(5,2,9);
        plot(squeeze(state_post_mem(5,1,:)))
        set(gca,'FontSize',AT_size);
        xlabel('Time (hrs)')
        ylabel('r_t (km)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        hold on
        plot(squeeze(state_smooth(5,1,:)))
        plot(squeeze(F_smooth(5,1,:)))

        subplot(5,2,2);
        plot(squeeze(state_post_mem(6,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('X_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(state_smooth(6,1,:)))
        plot(squeeze(F_smooth(6,1,:)))

        subplot(5,2,4);
        plot(squeeze(state_post_mem(7,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Y_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(state_smooth(7,1,:)))
        plot(squeeze(F_smooth(7,1,:)))

        subplot(5,2,6);
        plot(squeeze(state_post_mem(8,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Z_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(state_smooth(8,1,:)))
        plot(squeeze(F_smooth(8,1,:)))

        subplot(5,2,8);
        plot(squeeze(state_post_mem(9,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('\theta_d_o_t (rad/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(state_smooth(9,1,:)))
        plot(squeeze(F_smooth(9,1,:)))

        subplot(5,2,10);
        plot(squeeze(state_post_mem(10,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('r_t_,_d_o_t')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        %set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(state_smooth(10,1,:)))
        plot(squeeze(F_smooth(10,1,:)))
        xlabel('Time (hrs)')
        legend('EKF State', 'Smooth State', 'F_k Estimate')

    figure(2)
    clf
    subplot(5,2,1);
        plot(squeeze(P_post_mem(1,1,:)))
        set(gca,'FontSize',AT_size);
        title('EKF State Covariances')
        ylabel('X(m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(P_smooth(1,1,:)))
        %plot(squeeze(P_lag(1,1,:)))

        subplot(5,2,3);
        plot(squeeze(P_post_mem(2,2,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Y(m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(P_smooth(2,2,:)))
        %plot(squeeze(P_lag(2,2,:)))

        subplot(5,2,5);
        plot(squeeze(P_post_mem(3,3,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Z (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(P_smooth(3,3,:)))
        %plot(squeeze(P_lag(3,3,:)))

        subplot(5,2,7);
        plot(squeeze(P_post_mem(4,4,:)))
        set(gca,'FontSize',AT_size);
        ylabel('\theta (rad)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(P_smooth(4,4,:)))
       % plot(squeeze(P_lag(4,4,:)))

        subplot(5,2,9);
        plot(squeeze(P_post_mem(5,5,:)))
        set(gca,'FontSize',AT_size);
        xlabel('Time (hrs)')
        ylabel('r_t (km)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        hold on
        plot(squeeze(P_smooth(5,5,:)))
       % plot(squeeze(P_lag(5,5,:)))

        subplot(5,2,2);
        plot(squeeze(P_post_mem(6,6,:)))
        set(gca,'FontSize',AT_size);
        ylabel('X_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(P_smooth(6,6,:)))
       % plot(squeeze(P_lag(6,6,:)))

        subplot(5,2,4);
        plot(squeeze(P_post_mem(7,7,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Y_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(P_smooth(7,7,:)))
        %plot(squeeze(P_lag(7,7,:)))

        subplot(5,2,6);
        plot(squeeze(P_post_mem(8,8,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Z_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(P_smooth(8,8,:)))
        %plot(squeeze(P_lag(8,8,:)))

        subplot(5,2,8);
        plot(squeeze(P_post_mem(9,9,:)))
        set(gca,'FontSize',AT_size);
        ylabel('\theta_d_o_t (rad/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(P_smooth(9,9,:)))
        %plot(squeeze(P_lag(9,9,:)))

        subplot(5,2,10);
        plot(squeeze(P_post_mem(10,10,:)))
        set(gca,'FontSize',AT_size);
        ylabel('r_t_,_d_o_t')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        %set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(P_smooth(10,10,:)))
        %plot(squeeze(P_lag(10,10,:)))
        legend('EKF', 'Smoothed')%, 'Lag-One')
        xlabel('Index (1-30)')

    figure(3)
    clf
    %title('EKF Estimated Measurements')
        subplot(4,2,1);
        plot(squeeze(Zest_mem(1,1,:)))
        set(gca,'FontSize',AT_size);
        title('EKF Estimates Measurements')
        ylabel('X(m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(Zest_smooth(1,1,:)))

        subplot(4,2,3);
        plot(squeeze(Zest_mem(2,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Y(m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(Zest_smooth(2,1,:)))

        subplot(4,2,5);
        plot(squeeze(Zest_mem(3,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Z (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(Zest_smooth(3,1,:)))

        subplot(4,2,7);
        plot(squeeze(Zest_mem(4,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('\theta (rad)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(Zest_smooth(4,1,:)))

        subplot(4,2,2);
        plot(squeeze(Zest_mem(5,1,:)))
        set(gca,'FontSize',AT_size);
        xlabel('Index')
        ylabel('X_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        hold on
        plot(squeeze(Zest_smooth(5,1,:)))

        subplot(4,2,4);
        plot(squeeze(Zest_mem(6,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Y_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(Zest_smooth(6,1,:)))

        subplot(4,2,6);
        plot(squeeze(Zest_mem(7,1,:)))
        set(gca,'FontSize',AT_size);
        ylabel('Z_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        plot(squeeze(Zest_smooth(7,1,:)))
        xlabel('Index')
    end
    %}
    %===============================================
    
    %STEP 3 - Maximization Step, to obtain estimates of theta_j = [Q, R]
    N_start = 3;  %Loop from k >= 2 because of 1-Lag
    BQ = zeros(7,7,N);
    
    %Q-Adaptations
    if (Q_adapt_flag)
        for  k = N_start:1:N_use
            BQ(:,:,k) = (Zm_mem(:,:,k)-Zest_smooth(:,:,k))*(Zm_mem(:,:,k)-Zest_smooth(:,:,k))';
        end
        
        BQ_sum = sum(BQ(:,:,N_start:N_use),3);
        Q_update = 1/(N_use-N_start+1)*K*BQ_sum*K';
    end
    
    %R-Adaptations
    if (R_adapt_flag)
        BR = zeros(7,7,N);
        if ~Q_adapt_flag             %If BQ was not calculated in Q-adaptations
            for  k = N_start:1:N_use
                BQ(:,:,k) = (Zm_mem(:,:,k)-Zest_smooth(:,:,k))*(Zm_mem(:,:,k)-Zest_smooth(:,:,k))';
                %BR(:,:,k) = (Zm_mem(:,:,k)-Zest_smooth(:,:,k))*(Zm_mem(:,:,k)-Zest_smooth(:,:,k))' + H*P_smooth(:,:,k)*H';
                BR(:,:,k) = BQ(:,:,k) + H*P_smooth(:,:,k)*H';
            end
        else                    %BQ used from previous Q-adaptations
            for  k = N_start:1:N_use
                %BR(:,:,k) = (Zm_mem(:,:,k)-Zest_smooth(:,:,k))*(Zm_mem(:,:,k)-Zest_smooth(:,:,k))' + H*P_smooth(:,:,k)*H';
                BR(:,:,k) = BQ(:,:,k) + H*P_smooth(:,:,k)*H';
            end
        end
        BR_sum = sum(BR(:,:,N_start:N_use),3);
        R_update = 1/(N_use-N_start+1)*BR_sum; %  + H*P_post*H'; If only using latest term    
    end
    
    % =====================================================================
    %Diagonalize
    Q_update = diag(diag(Q_update));
    R_update = diag(diag(R_update));
    
    % =====================================================================
    %Positive Semidefinite Symmetric (PSDS) Check & Correction
    %Q_update = findPSDS(Q_update);
    %R_update = findPSDS(R_update);
    
    end
end

% =========================================================================
%% Converting data into output vector format
P_post_out = [ P_post(1,:) P_post(2,:) P_post(3,:) P_post(4,:) P_post(5,:) P_post(6,:)...
           P_post(7,:) P_post(8,:) P_post(9,:) P_post(10,:) ]; %(1 x 100)
       
Pr_theo_out = [ Pr_theo(1,:) Pr_theo(2,:) Pr_theo(3,:) Pr_theo(4,:) Pr_theo(5,:) Pr_theo(6,:)...
           Pr_theo(7,:)]; %(1 x 49)
       
K_out = [ K(1,:) K(2,:) K(3,:) K(4,:) K(5,:) K(6,:) K(7,:) K(8,:) K(9,:) K(10,:)]; %(1 x 70)
       
EKF_out = [ state_post' P_post_out Z_est' Pr_theo_out residuals' K_out]'; %(1x243)

% =========================================================================
%% Preparing data vector for Memory Storage
P_priori_out = [ P_priori(1,:) P_priori(2,:) P_priori(3,:) P_priori(4,:) P_priori(5,:) P_priori(6,:)...
           P_priori(7,:) P_priori(8,:) P_priori(9,:) P_priori(10,:) ]; %(1 x 100)

Phi_out = [ Phi(1,:) Phi(2,:) Phi(3,:) Phi(4,:) Phi(5,:) Phi(6,:)...
           Phi(7,:) Phi(8,:) Phi(9,:) Phi(10,:) ]; %(1 x 100)

       
Storage_out = [state_post' state_priori' P_post_out P_priori_out Phi_out ...
               Z_m' Z_est']'; %(1 x 334)

end
% =========================================================================