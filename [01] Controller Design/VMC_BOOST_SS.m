clc
clear all
close all

% Defining given parameters and matrices in MATLAB
syms d R RC Vg IL VC real

% Defining state vector X and input vector U
X = [IL; VC];
U = [Vg];

% Defining the matrices from the given state-space representation
A = [-(1-d)* (R * RC / (R + RC)), -(1-d) * R / (R + RC); (1-d)* R / (R + RC), -1 / (R + RC)];
B = [1; 0];


% Defining the output matrices (assuming output equation is similar)
C = [d * (R * RC / (R + RC)), R / (R + RC)];
D = 0;

% System parameters
% d = 0.5; % Duty cycle
% R = 1; % Load resistance
% RC = 10e-3; % Capacitor resistance
% Vg = 10; % Input voltage

D_num = double(subs(D, [d, R, RC], [0.5, 1, 0]))
A_num = double(subs(A, [d, R, RC], [0.5, 1, 0]))
B_num = double(subs(B, [d, R, RC], [0.5, 1, 0]))
C_num = double(subs(C, [d, R, RC], [0.5, 1, 0]))


% Defining numeric state-space representation
sys = ss(A_num, B_num, C_num, D_num);

% Displaying the state-space system
sys

% Now let's design a PID controller for the state-space system
% First, let's check the controllability of the system
Co = ctrb(sys);
rank(Co);

% The system is controllable, so we can design a state feedback controller
% Let's design a PID controller
Kp = 0;
Ki = 0;
Kd = 0;
C = pid(Kp, Ki, Kd);

% Now let's design a state feedback controller
sys_cl = feedback(sys, C);

% Displaying the closed-loop system
sys_cl

% Plotting the step response of the system with and without controller
figure;
step(sys_cl); % Displaying the step response information of the system
hold on;
step(sys); % Displaying the step response information of the system without controller


% Plotting the step response of each state
figure;
[Y, T, X] = step(sys_cl);
for i = 1:size(X, 2)
    subplot(size(X, 2), 1, i);
    plot(T, X(:, i));
    title(['Step Response of State ', num2str(i)]);
    xlabel('Time (s)');
    ylabel(['State ', num2str(i)]);
end