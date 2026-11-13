% This file simulates the relative three body EOMs for the Jupiter,
% Spacecraft, and Europa system.
% 
% HW 2 Problem 3
%
% Alec Lewandowski
% 16 Sep 2025

function HW2_problem3()

% ICs and simulation parameters for the Europa/Spacecraft/Jupiter
% simulation
eps = 1e-12;
M3 = 2e27; M1 = 5e22; M2 = 1000; G = 6.67e-11;
Re = 6.7e8; Rsc = 6e6;
tSim = [0 8]*24*3600;
IC03 = [0; 0; 0; 0; 0; 0];
IC01 = [Re; 0; 0; 0; sqrt(G*M3/Re); 0];
RelIC = [Rsc; 1; 0; 0; 0; sqrt(G*M1/Rsc)];
IC02 = IC01 + RelIC;
init = [IC01; IC02; IC03];
initrel = [IC01; IC02; IC03; RelIC];


% Same trajectory validation for both ode models
options = odeset('RelTol',eps,'AbsTol',eps);
[t,y] = ode45(@(t,y)RelThreeBodyProp(t,y,G,M1,M2,M3),tSim,initrel,options);
r12rel = [y(:,19) y(:,20) y(:,21)];
r1 = [y(:,1) y(:,2) y(:,3)];
r2 = [y(:,7) y(:,8) y(:,9)];
r12 = r2 - r1;
v12rel = [y(:,22) y(:,23) y(:,24)];
v1 = [y(:,4) y(:,5) y(:,6)];
v2 = [y(:,10) y(:,11) y(:,12)];
v12 = v2 - v1;

figure(1)
hold on
plot3(r12(:,1), r12(:,2), r12(:,3),'r')
plot3(r12rel(:,1), r12rel(:,2), r12rel(:,3),'--k')
hold off
view(-40,30)
xlabel('X-Direction [km]')
ylabel('Y-Direction [km]')
zlabel('Z-Direction [km]')
legend('Relative Trajectory','Relative State','Location','northwest')
title('Europa and Spacecraft Relative Motion Three-Body Sim: Problem 3.i')


%% 3.i
rel = r12rel-r12;
normrel = sqrt(rel(:,1).^2 + rel(:,2).^2 + rel(:,3).^2);

figure(2)
plot(t, normrel, 'r')
yscale log
xlabel('Time [s]')
ylabel('Norm of the Relative Trajectory Errors')
title('Relative Difference Norm Between Relative Trajectory and State')


%% 3.ii
vrel = v12rel-v12;
normvrel = sqrt(vrel(:,1).^2 + vrel(:,2).^2 + vrel(:,3).^2);

figure(3)
plot(t, normvrel, 'r')
yscale log
xlabel('Time [s]')
ylabel('Norm of the Relative Trajectory Errors')
title('Relative Difference Norm Between Velocities')