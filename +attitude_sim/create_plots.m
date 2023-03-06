function [] = create_plots(OutputData)
% Create Plots of Attitude Data ===========================================
% Description: This function takes the output data of the attitude
% simulator and creates several figures.
%
% Inputs:
%   OutputData - struct containing the output data
%
% Outputs:
%   None
%
% Created by: Cory Fraser - JAN 13, 2023
% Last Edits: Cory Fraser - JAN 13, 2023
% Copyright(c) 2023 by Cory Fraser
% =========================================================================
%% Plot Results
time_vec    = OutputData.time_vec;
qerr        = OutputData.qerr;
x           = OutputData.x;
torque      = OutputData.torque;

% Plot the Orbit
figure
plot3(x(:,8), x(:,9), x(:,10))

id_fig = 1;
figure%(id_fig); id_fig = id_fig + 1;
subplot(411)
plot(time_vec,qerr(:,1))
axis([time_vec(1) time_vec(end) -1.2 1.2])
set(gca,'fontsize',12)
ylabel('dq1')
subplot(412)
plot(time_vec,qerr(:,2))
axis([time_vec(1) time_vec(end) -1.2 1.2])
set(gca,'fontsize',12)
ylabel('dq2')
subplot(413)
plot(time_vec,qerr(:,3))
axis([time_vec(1) time_vec(end) -1.2 1.2])
set(gca,'fontsize',12)
ylabel('dq3')
subplot(414)
plot(time_vec,qerr(:,4))
axis([time_vec(1) time_vec(end) -1.2 1.2])
set(gca,'fontsize',12)
ylabel('dq4')
xlabel('Time (Sec)')

figure%(id_fig); id_fig = id_fig + 1;
subplot(311)
plot(time_vec,x(:,5))
set(gca,'fontsize',12)
ylabel('w1 (rad/sec)')
subplot(312)
plot(time_vec,x(:,6))
set(gca,'fontsize',12)
ylabel('w2 (rad/sec)')
subplot(313)
plot(time_vec,x(:,6))
set(gca,'fontsize',12)
ylabel('w3 (rad/sec)')
xlabel('Time (Sec)')

figure%(id_fig); id_fig = id_fig + 1;
subplot(311)
plot(time_vec,torque(:,1))
set(gca,'fontsize',12)
ylabel('L1 (N-m)')
subplot(312)
plot(time_vec,torque(:,2))
set(gca,'fontsize',12)
ylabel('L2 (N-m)')
subplot(313)
plot(time_vec,torque(:,3))
set(gca,'fontsize',12)
ylabel('L3 (N-m)')
xlabel('Time (Sec)')

%% Pose Plot Animation
%{
figure
qs = Q_actual_row(1); Q_actual_row(end);
ps = r_out(1,:); %[0 0 0];
%p_vec = %r_out/1000;
patch = poseplot(qs,ps, 'ENU');
%zlim([1.1*min(p_vec(:,3)) 1.1*max(p_vec(:,3))])
%xlim([1.1*min(p_vec(:,1)) 1.1*max(p_vec(:,1))])
%ylim([1.1*min(p_vec(:,2)) 1.1*max(p_vec(:,2))])
%xlim([1.1*min(p_vec(:,1)) 1.1*max(p_vec(:,1))])
title('Satellite Attitude in Body Frame')
xlabel("East-x (m)")
ylabel("North-y (m)")
zlabel("Up-z (m)");
for iStep = 1:length(Q_actual_row)
    %q = slerp(qs,qf,coeff);
    position = [0 0 0]; %p_vec(iStep,:) %ps + (pf - ps)*coeff;
    set(patch,Orientation=Q_actual_row(iStep),Position=position); 
    drawnow
    %pause(0.1)
end
%}

end