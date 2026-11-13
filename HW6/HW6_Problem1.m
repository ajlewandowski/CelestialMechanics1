% This file runs solution for HW6 Problem 1
%
% Alec Lewandowski
% 12/1/25
clear variables
close all
clc
format long

tau = 0.63071926935591928D+01;
fig = 1;

IC = [0.29166723152821866D+00;
    0.72969718788751359D-01;
    0.00000000000000000D+00;
    0.12873925508436218D+00;
    -0.47477120552157931D+00;
    0.18669082434116795D+01];
ICSTM = eye(6);

% Setting up and integrating the EOMs in cartesian coords
tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',tol);
[tvec, traj] = ode45(@(t,traj)HillThreeBodyPropNorm(t,traj),[0 tau],[IC; ICSTM(:)],options);
phiFinal = reshape(traj(end,7:end),[6,6]);

figure(fig)
hold on
plot3(traj(:,1),traj(:,2),traj(:,3),'r')
plot3(IC(1),IC(2),IC(3),'b*')
plot3(traj(end,1),traj(end,2),traj(end,3),'ks')
hold off
xlabel('X Position')
ylabel('Y Position')
zlabel('Z Position')
legend('Trajectory','Start Pt','End Pt')
title('Trajectory of the Hill Three-Body Problem')
view(3)
grid on
fig = fig + 1;

%% Problem 1.ii
% Calculating error in final state to see if orbit is periodic
finalState = traj(end, 1:3);
error = norm(finalState - IC(1:3));
disp(['Error in final state: ', num2str(error)]);

%% Problem 1.iii
fprintf('\nEigenvalues of the Monodromy Matrix:\n')
eigMono = eigs(phiFinal)

k1 = eigMono(1) + 1/eigMono(2);
magk1 = abs(k1);
k2 = eigMono(3) + 1/eigMono(4);
magk2 = abs(k2);

disp(['K1 stability value: ', num2str(k1)]);
disp(['K1 stability magnitude: ', num2str(magk1)]);
disp(['K2 stability value: ', num2str(k2)]);
disp(['K2 stability magnitude: ', num2str(magk2)]);

fprintf('\nBecause k1 and k2 are not both real and less than 2, the orbit is NOT stable!\n\n')

%% Problem 1.iv
figure(fig)
subplot(1,3,1)
hold on
plot(traj(:,1),traj(:,2),'r')
plot(IC(1),IC(2),'b*')
plot(traj(end,1),traj(end,2),'ks')
hold off
axis equal
xlabel('X-Pos')
ylabel('Y-Pos')
legend('Traj','Start','End')
subplot(1,3,2)
hold on
plot(traj(:,1),traj(:,3),'r')
plot(IC(1),IC(3),'b*')
plot(traj(end,1),traj(end,3),'ks')
hold off
axis equal
xlabel('X-Pos')
ylabel('Z-Pos')
legend('Traj','Start','End')
subplot(1,3,3)
hold on
plot(traj(:,2),traj(:,3),'r')
plot(IC(2),IC(3),'b*')
plot(traj(end,2),traj(end,3),'ks')
hold off
axis equal
xlabel('Y-Pos')
ylabel('Z-Pos')
legend('Traj','Start','End')
fig = fig + 1;

%% Problem 1.v
Cj = zeros(length(tvec),1);
for i = 1:length(tvec)
    x = traj(i,1);
    y = traj(i,2);
    z = traj(i,3);
    vx = traj(i,4);
    vy = traj(i,5);
    vz = traj(i,6);
    v = sqrt(vx^2 + vy^2 + vz^2);
    r = sqrt(x^2 + y^2 + z^2);
    Cj(i) = 3*x^2 - z^2 + 2/r - v^2;
end

cMean = mean(Cj);

figure(fig)
plot(tvec, Cj-cMean,'r')
xlabel('Time [TU]')
ylabel('Jacobi Int')
title(['Jacobi Integral vs Time, Cmean = ',num2str(cMean),'e-13'])
fig = fig + 1;

%% Problem 1.vi
semiMaj = zeros(length(tvec),1);
ecc = zeros(length(tvec),1);
inc = zeros(length(tvec),1);
LAN = zeros(length(tvec),1);
aop = zeros(length(tvec),1);
anom = zeros(length(tvec),1);

for i = 1:length(tvec)
    r = traj(i,1:3);
    v = traj(i,4:6);
    [oe] = eci2oe(r, v, 1);
    semiMaj(i) = oe(1);
    ecc(i) = oe(2);
    inc(i) = oe(3);
    aop(i) = oe(4);
    LAN(i) = oe(5);
    anom(i) = oe(6);
end

% Creating a subplot to show all the orbital elements
figure(fig)
subplot(6, 1, 1)
plot(tvec, semiMaj, 'k')
title('Semi-Major Axis')
ylabel('Length Units')
subplot(6, 1, 2)
plot(tvec, ecc, 'k')
title('Eccentricity')
ylabel('Unitless')
subplot(6, 1, 3)
plot(tvec, (inc), 'k')
title('Inclination')
ylabel('Degrees')
subplot(6, 1, 4)
plot(tvec, (aop), 'k')
title('Argument of Periapsis')
ylabel('Degrees')
subplot(6, 1, 5)
plot(tvec, (LAN), 'k')
title('Right Ascending Node')
ylabel('Degrees)')
subplot(6, 1, 6)
plot(tvec, (anom), 'k')
title('Mean Anomaly')
xlabel('Time [s]')
ylabel('Degrees')
sgtitle('Orbital Elements of Hill-Sim rel to Body 2')
fig = fig + 1;

%% Problem 1.vii
dPrimaries = 3.844000D+05;
gmE = 3.986004415D+05;
gmM = 4902.801076D+0;

nMoon = sqrt(gmE/dPrimaries^3);
moonHillRad = dPrimaries*(gmM/(gmE))^(1/3);

% Using the earth-moon distancing
LU = moonHillRad;
TU = 1/nMoon;

fprintf('Initial conditions in Earth-moon\n')
% Now recalculate the ICs
ICNew = [IC(1)*LU; IC(2)*LU; IC(3)*LU; IC(4)*(LU/TU); IC(5)*(LU/TU); IC(6)*(LU/TU)]
tauNew = tau*TU

%% Problem 1.viii
% Converting from rotating to intertial frame, pos same, vel diff
ICNew(4:6) = [ICNew(4) + nMoon; ICNew(5) + nMoon; ICNew(6)];

rIC = [ICNew(1); ICNew(2); ICNew(3)];
vIC = [ICNew(4); ICNew(5); ICNew(6)];

[oes] = eci2oe(rIC, vIC, gmM);

fprintf('Orbital Element ICs in Moon Hill EOMs \n')
fprintf('[a, ecc, inc, RAAN, aop, anom]\n')
oes'