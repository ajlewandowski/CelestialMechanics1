% This function runs solution for HW3 Problem 1
%
% Alec Lewandowski
% 10/6/25
clear variables
close all
clc

% Defining initial conditions and parameters
dhat = [0.2; 0.7; sqrt(0.47)];
gammaSRP = 0.03;
mu = 1;
a = 1.1; e = 0.1; i = 40; aop = 10; LAN = 100; nu = 75;
oeVec = [a, e, deg2rad(i), deg2rad(aop), deg2rad(LAN), deg2rad(nu)];

ti = 0; tf = 300;

% Part i: reporting initial r and v
[ri,vi] = oe2eciNU(oeVec, mu)

tol = 1e-8;
options = odeset('RelTol',tol,'AbsTol',tol);
[ta, Y] = ode45(@(t,y)TwoBodySRP(t,y,gammaSRP,mu,dhat),[ti tf],[ri; vi],options);
Y = real(Y);

% Plotting the projections
figure(1)
hold on 
plot3(Y(:,1),Y(:,2),Y(:,3),'r')
plot3(Y(1,1),Y(1,2),Y(1,3),'ko')
plot3(Y(end,1),Y(end,2),Y(end,3),'bo')
hold off
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
legend('Trajectory','Start Pt','End Pt')
title('Trajectory')

figure(2)
subplot(3,1,1)
hold on
plot(Y(:,1),Y(:,2),'c--')
plot(0,0,'ko')
plot(Y(1,1),Y(1,2),'ro')
hold off
title('XY')
ylabel('Y')
xlabel('X')
legend('Projection','Center','Starting Pt')
subplot(3,1,2)
hold on
plot(Y(:,1),Y(:,3),'c--')
plot(0,0,'ko')
plot(Y(1,1),Y(1,3),'ro')
hold off
title('XZ')
ylabel('Z')
xlabel('X')
legend('Projection','Center','Starting Pt')
subplot(3,1,3)
hold on
plot(Y(:,2),Y(:,3),'c--')
plot(0,0,'ko')
plot(Y(1,2),Y(1,3),'ro')
hold off
title('YZ')
ylabel('Z')
xlabel('Y')
legend('Projection','Center','Starting Pt')

% Plotting the orbital elements a, e, i, aop, LAN, M as function of time
oes = zeros(length(ta),6);
meanAnomA = zeros(length(ta),1);
for i = 1:length(ta)
    [oe] = eci2oeM([Y(i,1) Y(i,2) Y(i,3)], [Y(i,4) Y(i,5) Y(i,6)], mu);
    oes(i,:) = oe;
    meanAnomA(i) = oe(6);
end

oes(:,3:6) = wrapTo180(oes(:,3:6));

% Creating a subplot to show all the orbital elements
figure(3)
subplot(6, 1, 1)
plot(ta, oes(:,1), 'k')
title('Semi-Major Axis')
ylabel('Length Units')
subplot(6, 1, 2)
plot(ta, oes(:,2), 'k')
title('Eccentricity')
ylabel('Unitless')
subplot(6, 1, 3)
plot(ta, oes(:,3), 'k')
title('Inclination')
ylabel('Degrees')
subplot(6, 1, 4)
plot(ta, oes(:,4), 'k')
title('Argument of Periapsis')
ylabel('Degrees')
subplot(6, 1, 5)
plot(ta, oes(:,5), 'k')
title('Right Ascending Node')
ylabel('Degrees)')
subplot(6, 1, 6)
plot(ta, oes(:,6), 'k')
title('Mean Anomaly')
xlabel('Time [s]')
ylabel('Degrees')
sgtitle('Classical Orbital Elements of Simulation')

% Plotting Rsp and Rbarsp as function of time
Rsp = zeros(length(ta),1);
Rbarsp = zeros(length(ta),1);
for i = 1:length(ta)
    LAN = oes(i,5);
    aop = oes(i,4);
    inc = oes(i,3);
    ehat = [cosd(LAN)*cosd(aop) - sind(LAN)*sind(aop)*cosd(inc);
    sind(LAN)*cosd(aop) + cosd(LAN)*sind(aop)*cosd(inc);
    sind(aop)*sind(inc)];


    Rsp(i) = gammaSRP*dot(dhat,[Y(i,1) Y(i,2) Y(i,3)]);
    Rbarsp(i) = -1.5*oes(i,1)*gammaSRP*dot(dhat,ehat);
end

figure(4)
hold on
plot(ta, Rsp, 'r')
plot(ta, Rbarsp, 'k')
hold off
legend('Rsp','RBarsp')
xlabel('Time [TU]')
ylabel('R')
title('R and averaged R as functions of time')

% Plot the eccentricity vector projections
etraj = zeros(length(ta),3);
for i = 1:length(ta)
    h = cross(Y(i,1:3),Y(i,4:6));
    ecc = cross(Y(i,4:6), h) - Y(i,1:3)/norm(Y(i,1:3));
    etraj(i,:) = ecc;
end

figure(5)
plot3(etraj(:,1),etraj(:,2),etraj(:,3),'r', etraj(1,1),etraj(1,2),etraj(1,3),'ko', etraj(end,1),etraj(end,2),etraj(end,3),'bo')
xlabel('X')
ylabel('Y')
zlabel('Z')
legend('Ecc Trajectory','Start Pt','End Pt')
title('Eccentricity Trajectory')

figure(6)
subplot(3,1,1)
hold on
plot(etraj(:,1),etraj(:,2),'r')
plot(0,0,'ro')
plot(etraj(1,1),etraj(1,2),'go')
plot(etraj(end,1),etraj(end,2),'ko')
hold off
title('XY')
ylabel('Y')
xlabel('X')
legend('Projection','Center','Starting Pt','Ending Pt')
subplot(3,1,2)
hold on
plot(etraj(:,1),etraj(:,3),'r')
plot(0,0,'ro')
plot(etraj(1,1),etraj(1,3),'go')
plot(etraj(end,1),etraj(end,3),'ko')
hold off
title('XZ')
ylabel('Z')
xlabel('X')
legend('Projection','Center','Starting Pt','Ending Pt')
subplot(3,1,3)
hold on
plot(etraj(:,2),etraj(:,3),'r')
plot(0,0,'ro')
plot(etraj(1,2),etraj(1,3),'go')
plot(etraj(end,2),etraj(end,3),'ko')
hold off
title('YZ')
ylabel('Z')
xlabel('Y')
legend('Projection','Center','Starting Pt','Ending Pt')

%% PART B

% Defining initial conditions and parameters
E = 2*atan(tan(deg2rad(nu)/2)*(sqrt((1-e))/(1+e)));
M = E - e*sin(E);
oeVec = [a, e, deg2rad(i), deg2rad(aop), deg2rad(LAN), M];

tol = 1e-7;
options = odeset('RelTol',tol,'AbsTol',tol);
[tb, Y] = ode15s(@(t,y)TwoBodySRPinLPE(t,y,gammaSRP,mu,dhat),[ti tf],oeVec,options);
Y = real(Y);
Y(:,3:6) = wrapTo180(rad2deg(Y(:,3:6)));

for i = 1:length(tb)
    Y(i,2) = max(Y(i,2),tol);
end

% Part i: numerical answers
fprintf('dOE/dt at t0: \n')
[dydt] = TwoBodySRPinLPE(ti,Y(1,:),gammaSRP,mu,dhat)
fprintf('\n dOE/dt at tf: \n')
[dydt] = TwoBodySRPinLPE(tf,Y(end,:),gammaSRP,mu,dhat)

% Part ii overlay plot
% Creating a subplot to show all the orbital elements
figure(7)
subplot(6, 1, 1)
hold on
plot(ta, oes(:,1), 'k')
plot(tb, Y(:,1),'r')
hold off
title('Semi-Major Axis')
ylabel('Length Units')
legend('Cartesian Int','LPEs')
subplot(6, 1, 2)
hold on
plot(ta, oes(:,2), 'k')
plot(tb, Y(:,2),'r')
hold off
title('Eccentricity')
ylabel('Length Units')
legend('Cartesian Int','LPEs')
subplot(6, 1, 3)
hold on
plot(ta, oes(:,3), 'k')
plot(tb, Y(:,3),'r')
hold off
title('Inclination')
ylabel('Degrees')
legend('Cartesian Int','LPEs')
subplot(6, 1, 4)
hold on
plot(ta, oes(:,4), 'k')
plot(tb, Y(:,4),'r')
hold off
title('Arg of Periapsis')
ylabel('Degrees')
legend('Cartesian Int','LPEs')
subplot(6, 1, 5)
hold on
plot(ta, oes(:,5), 'k')
plot(tb, Y(:,5),'r')
hold off
title('LAN')
ylabel('Degrees')
legend('Cartesian Int','LPEs')
subplot(6, 1, 6)
hold on
plot(ta, oes(:,6), 'k')
plot(tb, Y(:,6),'r')
hold off
title('Mean Anomaly')
xlabel('Time [s]')
ylabel('Degrees')
legend('Cartesian Int','LPEs')

% Part iv: Showing rbarsp computed with LPEs is constant
RbarLPE = zeros(length(tb),1);
for i = 1:length(tb)
    LAN = Y(i,5);
    aop = Y(i,4);
    inc = Y(i,3);
    ehat = [cosd(LAN)*cosd(aop) - sind(LAN)*sind(aop)*cosd(inc);
    sind(LAN)*cosd(aop) + cosd(LAN)*sind(aop)*cosd(inc);
    sind(aop)*sind(inc)];

    RbarLPE(i) = -1.5*Y(i,1)*gammaSRP*dot(dhat,ehat);
end

figure(8)
hold on
plot(ta,abs(meanAnomA-mean(meanAnomA)),'r')
plot(tb,abs(Y(:,6)-mean(Y(:,6))),'b')
hold off
xlabel('Time [TU]')
ylabel('Mean Anomaly [deg]')
legend('Part a','Part b')
title('Difference from mean anomaly and its mean for parts a and b')


figure(9)
hold on
plot(ta, Rsp, 'r')
plot(ta, Rbarsp, 'k')
plot(tb, RbarLPE, 'b')
hold off
legend('Rsp','RBarsp','Rbarsp from LPEs')
xlabel('Time [TU]')
ylabel('SRP')
title('SRP and its mean as functions of time')

Rbarmean = mean(RbarLPE)


% Part v: Overlaying the eccentricity computed with LPEs and 2av
etrajnew = zeros(length(tb),3);
for i = 1:length(tb)
    [r,v] = oe2eciM([Y(i,1) Y(i,2) deg2rad(Y(i,3)) deg2rad(Y(i,4)) deg2rad(Y(i,5)) deg2rad(Y(i,6))],mu);
    h = cross(r,v);
    ecc = cross(v, h) - r/norm(r);
    etrajnew(i,:) = ecc;
end

figure(10)
subplot(3,1,1)
hold on
plot(etraj(:,1),etraj(:,2),'r')
plot(etrajnew(:,1),etrajnew(:,2),'b--')
plot(0,0,'ro')
plot(etraj(1,1),etraj(1,2),'go')
plot(etraj(end,1),etraj(end,2),'ko')
hold off
title('XY')
ylabel('Y')
xlabel('X')
legend('Cartesian','LPEs','Center','Starting Pt','Ending Pt')
subplot(3,1,2)
hold on
plot(etraj(:,1),etraj(:,3),'r')
plot(etrajnew(:,1),etrajnew(:,3),'b--')
plot(0,0,'ro')
plot(etraj(1,1),etraj(1,3),'go')
plot(etraj(end,1),etraj(end,3),'ko')
hold off
title('XZ')
ylabel('Z')
xlabel('X')
legend('Cartesian','LPEs','Center','Starting Pt','Ending Pt')
subplot(3,1,3)
hold on
plot(etraj(:,2),etraj(:,3),'r')
plot(etrajnew(:,2),etrajnew(:,3),'b--')
plot(0,0,'ro')
plot(etraj(1,2),etraj(1,3),'go')
plot(etraj(end,2),etraj(end,3),'ko')
hold off
title('YZ')
ylabel('Z')
xlabel('Y')
legend('Cartesian','LPEs','Center','Starting Pt','Ending Pt')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part vi: Rerun the simulation for a different set of ICs
