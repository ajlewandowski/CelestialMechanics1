% This function runs solution for HW5 Problem 1
%
% Alec Lewandowski
% 11/12/25
clear variables
close all
clc
format long

% Defining the ICs
r1 = [1; 0; 0];
v1 = [0; 1; 0];
STM1 = eye(6);
TOF = (2*pi)*(5/4);
GM = 1;
athrust = 5e-2;
fig = 1;

% Setting up and integrating the EOMs in cartesian coords
tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',tol);
[t1, y1] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM,athrust),[0 TOF],[r1; v1; STM1(:)],options);

%% Part A
figure(fig)
hold on
plot(y1(:,1), y1(:,2),'r')
plot(y1(1,1), y1(1,2),'b*')
plot(y1(end,1), y1(end,2),'ko')
hold off
xlabel('X')
ylabel('Y');
legend('Trajectory','Start Pt','End Pt')
title('Trajectory of the Thrust Propagation');
grid on;
fig = fig + 1;

sma1 = zeros(length(t1),1);
ecc1 = zeros(length(t1),1);
for i = 1:length(y1(:,1))
    [oe] = eci2oe(y1(i,1:3), y1(i,4:6), GM);
    sma1(i) = oe(1);
    ecc1(i) = oe(2);
end

%% Part B
figure(fig)
subplot(2,1,1)
plot(t1,sma1,'r')
xlabel('Time [TU]')
ylabel('Semi-Major Axis [LU]')
subplot(2,1,2)
plot(t1,ecc1,'b')
xlabel('Time [TU]')
ylabel('Eccentricity')
sgtitle('Problem 1.b: Semi-Major Axis and Eccentricity')
fig = fig + 1;

%% Part C
QSTMFinal = 0;
PhiFinal = reshape(y1(end,7:end), [6,6]);
for row = 1:6
    for col = 1:6
        QSTMFinal = QSTMFinal + (row*col + row + col)*PhiFinal(row,col);
    end
end
QSTMFinal

%% Part D
STM_oneone = y1(:,7);
figure(fig)
plot(t1, STM_oneone,'r')
xlabel('Time [TU]');
ylabel('STM(1,1)');
title('State Transition Matrix Element STM(1,1) Over Time');
grid on;
fig = fig + 1;

%% Part E
r1 = y1(1,1);
r2 = y1(1,2);
r3 = y1(1,3);
GInit = [3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r1 ^ 2 - GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r1 ^ 2 + athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.1e1 / 0.2e1) 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r2 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r2 * r1 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r1; 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r2 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r2 * r1 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r2 ^ 2 - GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r2 ^ 2 + athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.1e1 / 0.2e1) 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r2 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r2; 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r1 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r2 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r2 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 ^ 2 - GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 ^ 2 + athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.1e1 / 0.2e1);]

r1 = y1(end,1);
r2 = y1(end,2);
r3 = y1(end,3);
GFinal = [3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r1 ^ 2 - GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r1 ^ 2 + athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.1e1 / 0.2e1) 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r2 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r2 * r1 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r1; 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r2 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r2 * r1 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r2 ^ 2 - GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r2 ^ 2 + athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.1e1 / 0.2e1) 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r2 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r2; 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r1 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r2 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r2 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 ^ 2 - GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 ^ 2 + athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.1e1 / 0.2e1);]

QGi = 0;
QGf = 0;
for row = 1:3
    for col = 1:3
        QGi = QGi + (row*col + row + col)*GInit(row,col);
        QGf = QGf + (row*col + row + col)*GFinal(row,col);
    end
end
QGi
QGf