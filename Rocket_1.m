%% MATLAB Exericse I: Kalman Filtering
% For the course SC42025, Filtering and Identification
% Mike Pesselse [4300564]

%% Initialisation of data and parameters
clc; clear vars; close all;

load('rocket');

dt = 0.1;   % [s] time step 
m = 100;    % [kg] rocket mass
g = 9.81;   % [m/s^2] gravitational acceleration 
y_0 = 0;    % [m] altitude of the rocket
T = (length(y)-1)*dt; % [s] Final time  
t = 0:dt:T;  % [s] time steps 

%% 1) Simulation of the discrete rocket model
% Answer: The predicted trajectory does not coincide with the measured
% trajectory. This is because the simulation assumes that the drag on 
% the rocket is constant, which is not the case. The real drag
% force is of course highly nonlinear, depending on the velocity and local
% atmospheric conditons. 

% Discrete state-space model matrices
A_1 = [1 dt; 0 1];
B_1 = [dt^2/(2*m) -dt^2/2 -dt^2/(2*m); dt/m -dt -dt/m];
C_1 = [1 0];
D_1 = 0; 

% Create state-space model 
sys_discrete_1 = ss(A_1, B_1, C_1, D_1, dt); 

% Simulate state-space model
u_1 = u; 
[sim_discrete_1, t, x_1] = lsim(sys_discrete_1, u_1, t); 

% Plot the result of the simulation 
figure
subplot(2,1,1)
plot(t, y, 'r-', t, sim_discrete_1, 'b-');
title('Question 1: Preditced trajectory vs Measured trajectory ') 
xlabel('Time [s]')
ylabel('Position [m]')
legend('Measurement trajectory', 'Predicted Trajectory')

subplot(2,1,2)
plot(t, ydottrue, 'r-', t, x_1(:,2), 'b-')
title('Measured velocity vs Predicted velocity') 
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('True velocity', 'Predicted velocity')

%% 2) Simulation using a discrete asymptotic oberserver
% Answer: The estimation of the position trajectory of Q2 is much better 
% than the estimation of the position trajectory of Q1. The estimation
% for the velocity is much better, however very noisy. So both estimations
% have improved. 

% Pole placement 
poles = [0.8 0.7];
K_2 = place(A_1',C_1',poles)';

% Discrete state-space observer model matrices
A_2 = A_1 - K_2 * C_1;
B_2 = [B_1 K_2];             
C_2 = C_1;
D_2 = 0;

% Create state-space observer model
sys_observer_2 = ss(A_2, B_2, C_2, D_2, dt);

% Simulate state-space observer model 
    % The input of the observer consits of the signals u, d 
    % (again assumed constant), g and y.                          
u_2 = [u y];                
[sim_observer_2, t, x_observer_2] = lsim(sys_observer_2, u_2, t); 

% Plot the result of the simulation 
figure
subplot(2,1,1)
plot(t, ytrue, 'r-', t, sim_observer_2, 'b-')
title('Measured position vs Predicted position (observer)') 
xlabel('Time [s]')
ylabel('Position [m]')
legend('True position', 'Predicted position')

subplot(2,1,2)
plot(t, ydottrue, 'r-', t, x_observer_2(:,2), 'b-')
title('Measured velocity vs Predicted velocity (observer)') 
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('True velocity', 'Predicted velocity (observer)')

%% 3) Simulation using a Kalman filter to estimate the state 
% Answer: Using the Kalman filter the estimate for the velocity is less
% noisy. The position estimation is unchanged. We can also see that the
% tuned Q (Q_3_B = [100 0; 0 1];) is performing better, especially at later  
% time steps, then Q = identity (Q_3_A = [1 0; 0 1]), 
% however the tuned Q is a bit noisy. 

% Create kalman filter parameters using 5.7 in the book (p. 162 - 163) 
Q_3_A = [1 0; 0 1];         % Q = I to compare with the tuned Q
Q_3_B = [100 0; 0 1];       % This was tuned by trial and error 
R_3 = 1*10^3; 

[~, ~, ~, M_3_A, ~] = kalman(sys_discrete_1, Q_3_A, R_3);
[~, ~, ~, M_3_B, ~] = kalman(sys_discrete_1, Q_3_B, R_3); 

K_3_A = M_3_A;  % Kalman gain K based on Q is identity matrix
K_3_B = M_3_B;  % Kalman gain K based on tuned Q 

% Discrete state-space kalman filter model matrices
A_3_A = A_1 - K_3_A * C_1;  % System matrix A with Q is identity matrix
A_3_B = A_1 - K_3_B * C_1;  % System matrix A with tuned Q
B_3_A = [B_1 K_3_A];
B_3_B = [B_1 K_3_B];
C_3 = C_1;
D_3 = 0; 

% Create state-space kalman filter model
sys_kalman_3_A = ss(A_3_A, B_3_A, C_3, D_3, dt); % System with Q is identity 
sys_kalman_3_B = ss(A_3_B, B_3_B, C_3, D_3, dt); % System with tuned Q

% Simulate state-space kalman filter model
u_3 = [u y]; 
[sim_kalman_3_A, t, x_kalman_3_A] = lsim(sys_kalman_3_A, u_3, t);
[sim_kalman_3_B, t, x_kalman_3_B] = lsim(sys_kalman_3_B, u_3, t); 

% Plot the results of the simulation
figure
subplot(2,1,1)
plot(t, ytrue, 'r-', t, sim_kalman_3_A, 'b-',...
     t, sim_kalman_3_B, 'c-')
title('Measured position vs Predicted position with Kalman filter') 
xlabel('Time [s]')
ylabel('Position [m]')
legend('True position',...
    'Predicted position with Kalman filter and Q = I',...
    'Predicted position with Kalman filter and Q = tuned')

subplot(2,1,2)
plot(t, ydottrue, 'r-', t, x_kalman_3_A(:,2), 'b-',...
     t, x_kalman_3_B(:,2), 'c-')
title('Measured velocity vs Predicted velocity with kalman filter') 
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('True velocity',...
    'Predicted velocity with Kalman fitler and Q = 1',...
    'Predicted velocity with Kalman fitler and Q = tuned')


%% 4) Another way to estimate rocket velocity
% Answer: Another way to estimate the rocket velocity would be by differentiating
% the position y. However, this method does not produce a good result due
% to the noisy data. 

%Plot of velocity estimation by differntiating
figure 
plot(t, ydottrue, 'r-', t, [0; diff(y)], 'b-')
title('Measured velocity vs Predicted velocity by differentiation of y') 
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('True velocity', ... 
       'Predicted velocity by differentiating of the position y')
   
%% 5) Drag estimation
% Answer: The position comparison is still near perfect. The
% velocity prediction is slightly improved in comparison with the question
% 3, the kalman fitler without drag. The drag prediction is decent 
% incomparison with the actual drag. 

% Append drag model to the state-space model using Eq. 5.7 and ch 5.8
% in the book (p. 166 - 167) 
A_drag = [1 dt -1*dt^2/(2*m); 0 1 -dt/m; 0 0 1]; 
B_drag = [dt^2/(2*m) -1*dt^2/2 0; dt/m -dt 0; 0 0 1]; 
C_drag = [1 0 0];
D_drag = 0;
sys_drag_5 = ss(A_drag, B_drag, C_drag, D_drag, dt);

% Create kalman filter parameters including drag
Q_5 = [1 0 0; 0 1 0; 0 0 20];   % This was tuned by trial and error 
R_5 = 1*10^3; 

[~,~, ~, M_5, ~] = kalman(sys_drag_5, Q_5, R_5); 
K_5 = M_5;

% Discrete state-space kalman filter model matrices including drag
A_5 = A_drag - K_5 * C_drag; 
B_5 = [B_drag K_5];              
C_5 = C_drag;
D_5 = D_drag; 

% Create state-space kalman filter model including drag
sys_kalman_5 = ss(A_5, B_5, C_5, D_5, dt);

% Simulate state-space kalman filter model including drag
x_0 = [0 0 u(1,3)];
u_5 = [u(:,1:2) zeros(length(y), 1) y]; 
[sim_kalman_5, t, x_kalman_5] = lsim(sys_kalman_5, u_5, t, x_0); 

% Plot the results of the simulation
figure
subplot(3,1,1)
plot(t, ytrue, 'r-', t, sim_kalman_5, 'b-')
title('Measured position vs Predicted postion with kalman filter and drag') 
xlabel('Time [s]')
ylabel('Position [m')
legend('True position', 'Predicted velocity with Kalman fitler for drag')

subplot(3,1,2)
plot(t, ydottrue, 'r-', t, x_kalman_5(:,2), 'b-')
title('Measured velocity vs Predicted velocity with kalman filter for') 
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('True velocity', 'Predicted velocity with Kalman fitler for drag')

subplot(3,1,3)
plot(t, dtrue, 'r-', t, x_kalman_5(:,3), 'b-')
title('Measured drag vs Predicted drag with kalman filter') 
xlabel('Time [s]')
ylabel('Drag [N]')
legend('True drag', 'Predicted drag with kalman filter')

%% 6) Root-mean-square errors
% Answer: We can conclude that the results become better by including more
% estimations. 

% Q1) RMSE of Measurements
y_1_rmse = rms(ytrue - y);
v_1_rmse = rms(ydottrue - [0; diff(y)]);
d_1_rmse = rms(dtrue - u_1(:,3));

% Q2) RMSE of Asymptotic state observer
y_2_rmse = rms(ytrue - sim_observer_2);
v_2_rmse = rms(ydottrue - x_observer_2(:,2));
d_2_rmse = rms(dtrue - u_2(:,3));

% Q3) RMSE of Kalman filter for the tuned Q
y_3_rmse = rms(ytrue - sim_kalman_3_B);
v_3_rmse = rms(ydottrue - x_kalman_3_B(:,2));
d_3_rmse = rms(dtrue - u_3(:,3));

% Q5) RMSE of kalman filter including drag 
y_5_rmse = rms(ytrue - sim_kalman_5);
v_5_rmse = rms(ydottrue - x_kalman_5(:,2));
d_5_rmse = rms(dtrue - x_kalman_5(:,3));

% Produce table with results 
positions = [y_1_rmse, y_2_rmse, y_3_rmse, y_5_rmse]';
velocities = [v_1_rmse, v_2_rmse, v_3_rmse, v_5_rmse]';
drags = [d_1_rmse, d_2_rmse, d_3_rmse, d_5_rmse]';

Name = {'Measurements','Asymptotic state observer','Kalman filter',...
        'Kalman filter with drag'};
Var = {'Position', 'Velocity', 'Drag',};
Table = table( positions, velocities, drags, 'VariableNames', Var,...
        'RowNames', Name)
