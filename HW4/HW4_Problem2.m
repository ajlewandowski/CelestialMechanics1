% This function runs solution for HW4 Problem 2
%
% Alec Lewandowski
% 10/22/25
clear variables
close all
clc

format long

% Defining the parameters for the simulation, LU, TU
r0 = [1; 1; 0];
v0 = [0; 0.5; 0.3];
tSpan = [0 60];
mu = 1;
aPert = [0.01; 0.01; -0.01];
fig = 1;

%% PART A
% Setting up and integrating the EOMs in cartesian coords
tol = 1e-14;
options = odeset('RelTol',tol,'AbsTol',tol);
[t1, M1] = ode45(@(t,y)TwoBodyaPert(t,y,mu,aPert),tSpan,[r0; v0],options);

X1 = M1(:,1); Y1 = M1(:,2); Z1 = M1(:,3);
VX1 = M1(:,4); VY1 = M1(:,5); VZ1 = M1(:,6);

% Plotting the projections of the simulations
figure(fig)
subplot(3,1,1)
hold on
plot(X1, Y1, 'r')
plot(0, 0, 'k*')
plot(X1(1), Y1(1), 'ob')
hold off
xlabel('X')
ylabel('Y')
legend('XY Trajectory','Center','Starting Pt')
subplot(3,1,2)
hold on
plot(X1, Z1, 'r')
plot(0, 0, 'k*')
plot(X1(1), Z1(1), 'ob')
hold off
xlabel('X')
ylabel('Z')
legend('XZ Trajectory','Center','Starting Pt')
subplot(3,1,3)
hold on
plot(Y1, Z1, 'r')
plot(0, 0, 'k*')
plot(Y1(1), Z1(1), 'ob')
hold off
xlabel('Y')
ylabel('Z')
legend('YZ Trajectory','Center','Starting Pt')
sgtitle('Problem 2.a.i: Projection Plots')
fig = fig + 1;

% Converting and plotting in OEs now
oes = zeros(length(t1),6);
aPertXYZ = zeros(length(t1), 3);
for i = 1:length(t1)
    oes(i,:) = eci2oeNU(M1(i,1:3),M1(i,4:6),mu);

    LAN = oes(i,5);
    aop = oes(i,4);
    inc = oes(i,3);
    nu = oes(i,6);

    RLAN = [cos(LAN) -sin(LAN) 0;
        sin(LAN)  cos(LAN) 0;
        0 0 1];
    Ri = [1 0 0;
        0 cos(inc) -sin(inc);
        0 sin(inc)  cos(inc)];
    Raop = [cos(aop) -sin(aop) 0;
        sin(aop)  cos(aop) 0;
        0 0 1];
    Fnu = [cos(nu) -sin(nu) 0;
        sin(nu)  cos(nu) 0;
        0 0 1];

    % Perturbing acceleration in inertial frame from RTN components
    Q = RLAN*Ri*Raop*Fnu;
    aPertXYZ(i,:) = Q * aPert;
end

% Creating a subplot to show all the orbital elements
figure(fig)
subplot(6, 1, 1)
plot(t1, oes(:,1), 'k')
title('Semi-Major Axis')
ylabel('Length Units')
subplot(6, 1, 2)
plot(t1, oes(:,2), 'k')
title('Eccentricity')
ylabel('Unitless')
subplot(6, 1, 3)
plot(t1, rad2deg(oes(:,3)), 'k')
title('Inclination')
ylabel('Degrees')
subplot(6, 1, 4)
plot(t1, rad2deg(oes(:,4)), 'k')
title('Argument of Periapsis')
ylabel('Degrees')
subplot(6, 1, 5)
plot(t1, rad2deg(oes(:,5)), 'k')
title('Right Ascending Node')
ylabel('Degrees)')
subplot(6, 1, 6)
plot(t1, rad2deg(oes(:,6)), 'k')
title('True Anomaly')
xlabel('Time [s]')
ylabel('Degrees')
sgtitle('Problem 2.a.ii: Classical OEs')
fig = fig + 1;

% Plotting the XYZ of aPert as function of time
figure(fig)
subplot(3,1,1)
plot(t1, aPertXYZ(:,1), 'r')
xlabel('Time')
ylabel('aPert - X')
subplot(3,1,2)
plot(t1, aPertXYZ(:,2), 'r')
xlabel('Time')
ylabel('aPert - Y')
subplot(3,1,3)
plot(t1, aPertXYZ(:,3), 'r')
xlabel('Time')
ylabel('aPert - Z')
sgtitle('Problem 2.a.iii: aPert Plots')
fig = fig + 1;

%% PART B
% Setting up and integrating the EOMs in Gauss coords
oe0 = eci2oeNU(r0, v0, mu);

tol = 1e-14;
options = odeset('RelTol',tol,'AbsTol',tol);
[t2, M2] = ode45(@(t,y)TwoBodyaPertGauss(t,y,mu,aPert),t1,oe0,options);
a = M2(:,1); ecc = M2(:,2); inc = (M2(:,3)); aop = (M2(:,4)); LAN = (M2(:,5)); nu = (M2(:,6));

% Converting results to r and v and comparing with part a
rv = zeros(6,length(t2));
for i = 1:length(t2)
    [r,v] = oe2eciNU(M2(i,:),mu);
    rv(1:3,i) = r;
    rv(4:6,i) = v;
end

X2 = rv(1,:); Y2 = rv(2,:); Z2 = rv(3,:);
R2 = [X2; Y2; Z2];
VX2 = rv(4,:); VY2 = rv(5,:); VZ2 = rv(6,:);
V2 = [VX2; VY2; VZ2];

R1 = interp1(t1, M1(:,1:3), t2)';
V1 = interp1(t1, M1(:,4:6), t2)';

figure(fig)
subplot(3,1,1)
hold on
plot(X1, Y1, 'r')
plot(X2, Y2, 'b--')
plot(0, 0, 'k*')
plot(X1(1), Y1(1), 'ob')
hold off
xlabel('X')
ylabel('Y')
legend('Cartesian','Gauss','Center','Starting Pt')
subplot(3,1,2)
hold on
plot(X1, Z1, 'r')
plot(X2, Z2, 'b--')
plot(0, 0, 'k*')
plot(X1(1), Z1(1), 'ob')
hold off
xlabel('X')
ylabel('Z')
legend('Cartesian','Gauss','Center','Starting Pt')
subplot(3,1,3)
hold on
plot(Y1, Z1, 'r')
plot(Y2, Z2, 'b--')
plot(0, 0, 'k*')
plot(Y1(1), Z1(1), 'ob')
hold off
xlabel('Y')
ylabel('Z')
legend('Cartesian','Gauss','Center','Starting Pt')
sgtitle('Problem 2.b.ii: Projection Plots Comparison')
fig = fig + 1;

% Plotting the semilogy norms in the difference between position and
% velocity of the two simulation to compare
Rmat = R2 - R1;
Vmat = V2 - V1;
diffR = vecnorm(Rmat);
diffV = vecnorm(Vmat);

figure(fig)
hold on
semilogy(t2,diffR,'r')
semilogy(t2,diffV,'b')
hold off
xlabel('Time')
ylabel('Norm Difference Error')
legend('R-Diff','V-Diff')
title('Problem 2.b.iii: SemilogY Normed Differences')

% Printing doedt vector at ti and tf
fprintf('\ndOE/dt at ti:\n');  
disp(TwoBodyaPertGauss(t2(1), M2(1,:), mu, aPert'))
fprintf('\ndOE/dt at tf:\n');  
disp(TwoBodyaPertGauss(t2(end), M2(end,:), mu, aPert)')