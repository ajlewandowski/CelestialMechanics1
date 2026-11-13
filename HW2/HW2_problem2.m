% This file simulates a general 3-body problem assuming generic units. The
% problem continues into the Europa problem and analyzing periodic orbits.
% 
% HW 2 Problem 2
%
% Alec Lewandowski
% 16 Sep 2025

function HW2_problem2()

%%%%%%%%%%%%%%% PART A %%%%%%%%%%%%%%%

% ICs and the properties of the orbit progogator
eps = 1e-8;
tSim = [0 6];
M1 = 1; M2 = 2; M3 = 1; G = 1;
IC01 = [0.5; 0; 0; 0; 0; 0];
IC02 = [2; 0; 0; 0; 0; 0];
IC03 = [3.5; 0; 0; 0; 1; 0];
init = [IC01; IC02; IC03];

% Using ODE45 to propogate the orbit
options = odeset('RelTol',eps,'AbsTol',eps);
[t,y] = ode45(@(t,y)ThreeBodyProp(t,y,G,M1,M2,M3),tSim,init,options);
r1 = [y(:,1) y(:,2) y(:,3)];
r2 = [y(:,7) y(:,8) y(:,9)];
r3 = [y(:,13) y(:,14) y(:,15)];
v1 = [y(:,4) y(:,5) y(:,6)];
v2 = [y(:,10) y(:,11) y(:,12)];
v3 = [y(:,16) y(:,17) y(:,18)];

% Calculating the vector of the center of mass of the system
comx = (M1*y(:,1) + M2*y(:,7) + M3*y(:,13)) ./ (M1+M2+M3);
comy = (M1*y(:,2) + M2*y(:,8) + M3*y(:,14)) ./ (M1+M2+M3);

%% 2.a.i
% Plotting the trajectories of the three body simulation given in (a)
figure(1)
hold on
plot(y(:,1), y(:,2), 'r')
plot(y(:,7), y(:,8), 'b')
plot(y(:,13), y(:,14), 'k')
plot(comx, comy, 'g')
plot([IC01(1) IC02(1) IC03(1)],[IC01(2) IC02(2) IC03(2)],'*c')
hold off
xlabel('X-Direction [LU]')
ylabel('Y-Direction [LU]')
legend('Body 1','Body 2','Body 3','CoM','Starting Pts','Location','northwest')
title('Three-Body Simulation: Problem 2.a')

%% 2.a.ii
% Calculating the angular momentum of the system and plotting the
% components of the angular momentum
H1 = M1*cross(r1, v1);
H2 = M2*cross(r2, v2);
H3 = M3*cross(r3, v3);

Hsysx = H1(:,1) + H2(:,1) + H3(:,1);
Hsysy = H1(:,2) + H2(:,2) + H3(:,2);
Hsysz = H1(:,3) + H2(:,3) + H3(:,3);

Hsysxavg = mean(Hsysx);
Hsysyavg = mean(Hsysy);
Hsyszavg = mean(Hsysz);

figure(2)
subplot(3,1,1)
plot(t, Hsysx-Hsysxavg,'r')
title(['X-Angular Momentum of System: ' num2str(Hsysxavg)])
xlabel('Time [TU]')
ylabel('Hx-Sys Deviation')
subplot(3,1,2)
plot(t, Hsysy-Hsysyavg,'r')
title(['Y-Angular Momentum of System: ' num2str(Hsysyavg)])
xlabel('Time [TU]')
ylabel('Hy-Sys Deviation')
subplot(3,1,3)
plot(t, Hsysz-Hsyszavg,'r')
title(['Z-Angular Momentum of System: ' num2str(Hsyszavg)])
xlabel('Time [TU]')
ylabel('Hz-Sys Deviation')

%% 2.a.iii
% Calculating the energy of the system. In order to do so the magnitude
% values of the position and velocity vectors are needed
r1mag = sqrt(r1(:,1).^2 + r1(:,2).^2 + r1(:,3).^2);
r2mag = sqrt(r2(:,1).^2 + r2(:,2).^2 + r2(:,3).^2);
r3mag = sqrt(r3(:,1).^2 + r3(:,2).^2 + r3(:,3).^2);
v1mag = sqrt(v1(:,1).^2 + v1(:,2).^2 + v1(:,3).^2);
v2mag = sqrt(v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2);
v3mag = sqrt(v3(:,1).^2 + v3(:,2).^2 + v3(:,3).^2);

Tsys = 0.5*M1*v1mag.^2 + 0.5*M2*v2mag.^2 + 0.5*M3*v3mag.^2;
Usys = G*M1*M2./(abs(r2mag-r1mag)) + G*M1*M3./(abs(r3mag-r1mag)) + G*M2*M3./(abs(r3mag-r2mag));
Esys = Tsys - Usys;
Edev = Esys - mean(Esys);

figure(3)
subplot(3,1,1)
plot(t, Tsys,'r')
title('Total Kinetic Energy of the System')
xlabel('Time [TU]')
ylabel('Kinetic Energy')
subplot(3,1,2)
plot(t, Usys,'r')
title('Total Potential Energy of the System')
xlabel('Time [TU]')
ylabel('Potential Energy')
subplot(3,1,3)
plot(t, Edev,'r')
title(['Total Energy of the System E - avg(E): ' num2str(mean(Esys))])
xlabel('Time [TU]')
ylabel('System Energy')



%%%%%%%%%%%%%%% PART B %%%%%%%%%%%%%%%

%% 2.b.i & 2.b.ii
% ICs - by pure guess. Making initial conditions symmetric with eachother
% to view behavior and find periodic orbits. In order to check going to
% plot each component on a subplot. This orbit is a symmetric triangle
% where each vertex is 2 units from origin
eps = 1e-12;
M1 = 100; M2 = 100; M3 = 100; G = 1;
tSim = [0 50];
IC01 = [-2*cosd(120); -2*sind(120); 0; 5*cosd(30); -5*sind(30); 0];
IC02 = [0; 2; 0; -5; 0; 0];
IC03 = [2*cosd(120); -2*sind(120); 0; 5*sind(30); 5*cosd(30); 0];
init = [IC01; IC02; IC03];

comxi = (M1*IC01(1) + M2*IC02(1) + M3*IC03(1)) ./ (M1+M2+M3);
comyi = (M1*IC01(2) + M2*IC02(2) + M3*IC03(2)) ./ (M1+M2+M3);
comzi = (M1*IC01(3) + M2*IC02(3) + M3*IC03(3)) ./ (M1+M2+M3);
comi = [comxi comyi comzi];

% Running the simulation for the symmetric triangular case
options = odeset('RelTol',eps,'AbsTol',eps);
[t,y] = ode45(@(t,y)ThreeBodyProp(t,y,G,M1,M2,M3),tSim,init,options);
r1 = [y(:,1) y(:,2) y(:,3)];
r2 = [y(:,7) y(:,8) y(:,9)];
r3 = [y(:,13) y(:,14) y(:,15)];
v1 = [y(:,4) y(:,5) y(:,6)];
v2 = [y(:,10) y(:,11) y(:,12)];
v3 = [y(:,16) y(:,17) y(:,18)];

% Calculating the center-of-mass of the symmetric-triangular system
comx = (M1*y(:,1) + M2*y(:,7) + M3*y(:,13)) ./ (M1+M2+M3);
comy = (M1*y(:,2) + M2*y(:,8) + M3*y(:,14)) ./ (M1+M2+M3);
comz = (M1*y(:,3) + M2*y(:,9) + M3*y(:,15)) ./ (M1+M2+M3);
com = [comx comy comz];

% Plotting the orbits relative to the center-of-mass to check if it is
% symmetric and repeating
figure(4)
subplot(3,1,1)
hold on
plot(t, r1(:,1)-comx, 'r')
plot(t, r2(:,1)-comx, 'b')
plot(t, r3(:,1)-comx, 'k')
hold off
xlabel('Time [TU]')
ylabel('Distance [LU]')
legend('Body 1','Body 2','Body 3','Location','northeast')
title('X-Component of Three-Body System')
subplot(3,1,2)
hold on
plot(t, r1(:,2)-comy, 'r')
plot(t, r2(:,2)-comy, 'b')
plot(t, r3(:,2)-comy, 'k')
hold off
xlabel('Time [TU]')
ylabel('Distance [LU]')
legend('Body 1','Body 2','Body 3','Location','northeast')
title('Y-Component of Three-Body System')
subplot(3,1,3)
hold on
plot(t, r1(:,3)-comz, 'r')
plot(t, r2(:,3)-comz, 'b')
plot(t, r3(:,3)-comz, 'k')
hold off
xlabel('Time [TU]')
ylabel('Distance [LU]')
legend('Body 1','Body 2','Body 3','Location','northeast')
title('Z-Component of Three-Body System')
sgtitle('Relative to CoM Three-Body Periodic Orbit: Problem 2.b.i Symmetric')

figure(5)
hold on
plot(r1(:,1), r1(:,2), 'r')
plot(r2(:,1), r2(:,2), 'b')
plot(r3(:,1), r3(:,2), 'k')
hold off
xlabel('X [LU]')
ylabel('Y [LU]')
legend('Body 1','Body 2','Body 3','Location','northwest')
title('XY-Plane of Triangular-Symmetric Three-Body Orbit')

% Finding the first two instances of a local minimum in the data where the
% state position vectors most closely match the initial positions.
IC12 = IC02 - IC01;
IC32 = IC02 - IC03;
IC31 = IC01 - IC03;
r12 = r2 - r1;
r32 = r2 - r3;
r31 = r1 - r3;
for var = 2:length(t)
    tem1 = norm((r12(var,:)-IC12).^2);
    tem2 = norm((r32(var,:)-IC32).^2);
    tem3 = norm((r31(var,:)-IC31).^2);
    temmin(var-1) = norm([tem1 tem2 tem3]);
    indtemmin(var-1) = var-1;
end
[~, indmins] = findpeaks(1./temmin);
taustc1 = t(indmins(1) + 1);
taustc2 = t(indmins(2) + 1);

figure(20)
plot(indtemmin, temmin, 'r')
xlabel('Index')
ylabel('Value')
title('Norm of Difference Squared between State and ICs')

% Calculating the error in position after one periodic orbit. 
% Symmetric Triangular
errorstc1 = temmin(indmins(1) + 1);
errorstc2 = temmin(indmins(2) + 1);

% Source for a periodic orbit:
% Observable HQ Planar Three-Body Orbits
% https://observablehq.com/@rreusser/periodic-planar-three-body-orbits
% Using the intial conditions from the website above for a periodic
% three-body orbit
eps = 1e-12;
M1 = 1; M2 = 1; M3 = 1; G = 1;
tSim = [0 50];
IC01 = [-1; 0; 0; 0.306893; 0.125507; 0];
IC02 = [1; 0; 0; 0.306893; 0.125507; 0];
IC03 = [0; 0; 0; -0.613786; -0.251014; 0];
init = [IC01; IC02; IC03];

% Running the simulation
options = odeset('RelTol',eps,'AbsTol',eps);
[t,y] = ode45(@(t,y)ThreeBodyProp(t,y,G,M1,M2,M3),tSim,init,options);
r1 = [y(:,1) y(:,2) y(:,3)];
r2 = [y(:,7) y(:,8) y(:,9)];
r3 = [y(:,13) y(:,14) y(:,15)];
v1 = [y(:,4) y(:,5) y(:,6)];
v2 = [y(:,10) y(:,11) y(:,12)];
v3 = [y(:,16) y(:,17) y(:,18)];

% Calculating the center of mass of the system
comx = (M1*y(:,1) + M2*y(:,7) + M3*y(:,13)) ./ (M1+M2+M3);
comy = (M1*y(:,2) + M2*y(:,8) + M3*y(:,14)) ./ (M1+M2+M3);
comz = (M1*y(:,3) + M2*y(:,9) + M3*y(:,15)) ./ (M1+M2+M3);

% Plotting the individual components to check for the orbits periodic
% properties
figure(6)
subplot(3,1,1)
hold on
plot(t, r1(:,1)-comx, 'r')
plot(t, r2(:,1)-comx, 'b')
plot(t, r3(:,1)-comx, 'k')
hold off
xlabel('Time [TU]')
ylabel('Distance [LU]')
legend('Body 1','Body 2','Body 3','Location','northeast')
title('X-Component of Three-Body System')
subplot(3,1,2)
hold on
plot(t, r1(:,2)-comy, 'r')
plot(t, r2(:,2)-comy, 'b')
plot(t, r3(:,2)-comy, 'k')
hold off
xlabel('Time [TU]')
ylabel('Distance [LU]')
legend('Body 1','Body 2','Body 3','Location','northeast')
title('Y-Component of Three-Body System')
subplot(3,1,3)
hold on
plot(t, r1(:,3)-comz, 'r')
plot(t, r2(:,3)-comz, 'b')
plot(t, r3(:,3)-comz, 'k')
hold off
xlabel('Time [s]')
ylabel('Distance [km]')
legend('Body 1','Body 2','Body 3','Location','northeast')
title('Z-Component of Three-Body System')
sgtitle('Three-Body Periodic Orbit: Problem 2.b.i From Literature')

figure(7)
hold on
plot(r1(:,1), r1(:,2), 'r')
plot(r2(:,1), r2(:,2), 'b')
plot(r3(:,1), r3(:,2), 'k')
hold off
xlabel('X [LU]')
ylabel('Y [LU]')
legend('Body 1','Body 2','Body 3','Location','northwest')
title('XY-Plane of Symmetric Three-Body Orbit (from Literature)')

% Finding the first two instances of a local minimum in the data where the
% state position vectors most closely match the initial positions.
IC12 = IC02 - IC01;
IC32 = IC02 - IC03;
IC31 = IC01 - IC03;
r12 = r2 - r1;
r32 = r2 - r3;
r31 = r1 - r3;
temmin = zeros(1, length(t)-1);
indtemmin = zeros(1, length(1)-1);
for var = 2:length(t)
    tem1 = norm((r12(var,:)-IC12).^2);
    tem2 = norm((r32(var,:)-IC32).^2);
    tem3 = norm((r31(var,:)-IC31).^2);
    temmin(var-1) = norm([tem1 tem2 tem3]);
    indtemmin(var-1) = var-1;
end
[~, indmins] = findpeaks(1./temmin);
tauptb1 = t(indmins(10) + 1);
tauptb2 = t(indmins(20) + 1);

figure(21)
plot(indtemmin, temmin, 'r')
xlabel('Index')
ylabel('Value')
title('Norm of Difference Squared between State and ICs')

% Calculating the error in position after one periodic orbit. 
errorptb1 = temmin(indmins(10) + 1);
errorptb2 = temmin(indmins(20) + 1);

% Printing the results to the command window
fprintf('\n\n')
fprintf('The error and time from the initial state after one orbit in the triangular symmetric case are: \n')
fprintf('          %1.12f km     and     %5.12f sec \n',errorstc1, taustc1)
fprintf('The error from the initial state after two orbits in the triangular symmetric case is: \n')
fprintf('          %1.12f km     and     %5.12f sec \n',errorstc2, taustc2)


fprintf('\n\n')
fprintf('The error and time from the initial state after one orbit in symmetric case are: \n')
fprintf('          %1.12f km     and     %5.12f sec \n',errorptb1, tauptb1)
fprintf('The error from the initial state after two orbits in the symmetric case is: \n')
fprintf('          %1.12f km     and     %5.12f sec \n',errorptb2, tauptb2)



%%%%%%%%%%%%%%% PART C %%%%%%%%%%%%%%%

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

% Running the simulation
options = odeset('RelTol',eps,'AbsTol',eps);
[t,y] = ode45(@(t,y)ThreeBodyProp(t,y,G,M1,M2,M3),tSim,init,options);
r1 = [y(:,1) y(:,2) y(:,3)];
r2 = [y(:,7) y(:,8) y(:,9)];
r3 = [y(:,13) y(:,14) y(:,15)];
v1 = [y(:,4) y(:,5) y(:,6)];
v2 = [y(:,10) y(:,11) y(:,12)];
v3 = [y(:,16) y(:,17) y(:,18)];

% Calcuating the center-of-mass of the system
comx = (M1*y(:,1) + M2*y(:,7) + M3*y(:,13)) ./ (M1+M2+M3);
comy = (M1*y(:,2) + M2*y(:,8) + M3*y(:,14)) ./ (M1+M2+M3);
comz = (M1*y(:,3) + M2*y(:,9) + M3*y(:,15)) ./ (M1+M2+M3);

%% 2.c.i
% Plotting the 3D trajectories of the simulation and the XY projection
figure(8)
hold on
plot3(r1(:,1), r1(:,2), r1(:,3), 'r')
plot3(r2(:,1), r2(:,2), r2(:,3), 'k')
plot3(r3(:,1), r3(:,2), r3(:,3), 'b')
plot3(comx, comy, comz, 'g')
plot3([IC01(1) IC02(1) IC03(1)],[IC01(2) IC02(2) IC03(2)],[IC01(3) IC02(3) IC03(3)],'*c')
hold off
view(-40,30)
xlabel('X-Direction [km]')
ylabel('Y-Direction [km]')
zlabel('Z-Direction [km]')
legend('Europa','Spacecraft','Jupiter','CoM','Starting Pts','Location','northwest')
title('Jupiter-Centered Europa/Spacecraft Three-Body Sim: Problem 2.c.i')

figure(9)
hold on
plot(r1(:,1), r1(:,2), 'r')
plot(r2(:,1), r2(:,2), 'k')
plot(r3(:,1), r3(:,2), 'b')
plot(comx, comy, 'g')
plot([IC01(1) IC02(1) IC03(1)],[IC01(2) IC02(2) IC03(2)],'*c')
hold off
xlabel('X-Direction [LU]')
ylabel('Y-Direction [LU]')
legend('Europa','Spacecraft','Jupiter','CoM','Starting Pts','Location','northwest')
title('Jupiter-Centered Europe/Spacecraft XY Projection: Problem 2.c.i')

%% 2.c.ii
% Initializing arrays that store the orbital elements and calculating the
% relative motion vectors
r12 = r2 - r1;
v12 = v2 - v1;
lengthsim = length(r12);
semimaj = zeros(lengthsim, 1);
ecc = zeros(lengthsim, 1);
inc = zeros(lengthsim, 1);
argp = zeros(lengthsim, 1);
ascn = zeros(lengthsim, 1);
anom = zeros(lengthsim, 1);

% Iterating through the r and v vectors and using the eci2oe algorithm the
% determine the orbital elements
for var = 1:lengthsim
    matoe = eci2oe(r12(var,:), v12(var,:), G*(M1 + M3));
    % output is semi-major axis, eccentricity, inclination, arg of periapsis, 
    % ascending node, true anomaly
    semimaj(var) = matoe(1);
    ecc(var) = matoe(2);
    inc(var) = matoe(3);
    argp(var) = matoe(4);
    ascn(var) = matoe(5);
    anom(var) = matoe(6);
end

% Plotting the orbital elements of the spacecraft relative to Europa
figure(10)
subplot(2, 3, 1)
plot(t, semimaj, 'k')
title('Semi-Major Axis')
xlabel('Time [s]')
ylabel('Distance [km]')
subplot(2, 3, 2)
plot(t, ecc, 'k')
title('Eccentricity')
xlabel('Time [s]')
ylabel('Unitless')
subplot(2, 3, 3)
plot(t, inc, 'k')
title('Inclination')
xlabel('Time [s]')
ylabel('Degrees')
subplot(2, 3, 4)
plot(t, argp, 'k')
title('Argument of Periapsis')
xlabel('Time [s]')
ylabel('Degrees')
subplot(2, 3, 5)
plot(t, ascn, 'k')
title('Right Ascending Node')
xlabel('Time [s]')
ylabel('Degrees')
subplot(2, 3, 6)
plot(t, anom, 'k')
title('True Anomaly')
xlabel('Time [s]')
ylabel('Degrees')
sgtitle('Spacecraft Orbital Elements Relative to Europa')

%% 2.c.iii
% Finding the magnitude of the radii of the spacecraft and Europa relative
% to Jupiter. Plotting the magnitudes as a function of time
r1mag = sqrt(r1(:,1).^2 + r1(:,2).^2 + r1(:,3).^2);
r2mag = sqrt(r2(:,1).^2 + r2(:,2).^2 + r2(:,3).^2);

figure(11)
hold on
plot(t, r1mag,'b')
plot(t, r2mag,'k')
hold off
xlabel('Time [s]')
ylabel('Radial Magnitude [km]')
title('Radial Distance Magnitudes from Jupiter')
legend('Europa','Spacecraft')

%% 2.c.iv
% Plotting the relative radius magnitude to determine if the spacecraft
% hits Europa by its position vector being less than the radius of Europa
r12mag = sqrt(r12(:,1).^2 + r12(:,2).^2 + r12(:,3).^2);

figure(12)
hold on
plot(t, r12mag, 'r')
yline(1565e3)
hold off
xlabel('Time [s]')
ylabel('Radial Magnitude [km]')
legend('Spacecraft Radius to Europa','Europa Radius','Location','northwest')
title('Relative Radial Magnitude of the Spacecraft to Europa')

for k = 1:length(r12mag)
    if r12mag(k) <= 1565e3
        intind = k;
        break
    end
end

fprintf('The intersection point happens at time: %8.4f seconds \n\n', t(intind))