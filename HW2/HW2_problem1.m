% This file simulates the relative two body EOMs with a flight time of 1
% period. Initial conditions are stated below but are arbitrary in
% length/time units. 
% 
% HW 2 Problem 1
%
% Alec Lewandowski
% 16 Sep 2025

function HW2_problem1()

% Initial Two Body Parameters
mu = 1;
r0 = [1; 0; 0];
v0 = [1; 0; sqrt(0.5)];
init = [r0; v0];

% ODE Parameters and Tolerances
mags = -2:-1:-13;
tols = 1*10.^(mags);
tSpan = [0 1e3];
options = odeset('RelTol',tols(7),'AbsTol',tols(7),'Events', @CompletedPeriodEventFuncHW2P1);

% Calculating the specific energy, initial semi-major axis, and analytical
% solution for the time period
specE = 0.5*norm(v0)^2 - (mu / norm(r0));
smai = -mu/(2*specE);
tau = sqrt((4*pi^2)*smai^3/mu);

% Using ODE45 to porpogate orbit
[t,r,tp,~,~] = ode45(@(t,r)RelTwoBodyProp(t,r,mu),tSpan,init,options);

%% 1.i and 1.ii
% Plotting the 3D orbit and its XY projection
figure(1)
hold on
plot3(r(:,1), r(:,2), r(:,3), 'r')
plot3(0, 0, 0, 'ob')
plot3(r0(1), r0(2), r0(3), 'or')
plot3(r(end,1), r(end,2), r(end,3),'*b')
hold off
view(-40,30)
xlabel('X-Direction [LU]')
ylabel('Y-Direction [LU]')
zlabel('Z-Direction [LU]')
legend('Orbit','Origin','Start Pt','End Pt')
title('3D Orbit Plot: Problem 1.ii')

figure(2)
plot(r(:,1), r(:,2), 'r')
xlabel('X-Direction [LU]')
ylabel('Y-Direction [LU]')
title('3D Orbit XY Projection Plot: Problem 1.ii')

figure(3)
plot(r(:,1), r(:,3), 'r')
xlabel('X-Direction [LU]')
ylabel('Z-Direction [LU]')
title('3D Orbit XZ Projection Plot: Problem 1.ii')

% Reporting the analytical and event function solutions for time period
fprintf('\nProblem 1.i \n')
fprintf('The time of one period for this two-body orbit is: %5.16f (event function)\n', tp)
fprintf('The time of one period for this two-body orbit is: %5.16f (analytical equation)\n\n', tau)


%% 1.iii
% Initializing arrays that store the orbital elements
lengthsim = length(r(:,1));
semimaj = zeros(lengthsim, 1);
ecc = zeros(lengthsim, 1);
inc = zeros(lengthsim, 1);
argp = zeros(lengthsim, 1);
ascn = zeros(lengthsim, 1);
anom = zeros(lengthsim, 1);

% Iterating through all of the r and v vectors and converting them to
% orbital elements with the eci2oe algorithm
for var = 1:lengthsim
    matoe = eci2oe([r(var,1) r(var,2) r(var,3)], [r(var,4) r(var,5) r(var,6)], mu);
    % output is semi-major axis, eccentricity, inclination, arg of periapsis, 
    % ascending node, true anomaly
    semimaj(var) = matoe(1);
    ecc(var) = matoe(2);
    inc(var) = matoe(3);
    argp(var) = matoe(4);
    ascn(var) = matoe(5);
    anom(var) = matoe(6);
end

% Creating a subplot to show all the orbital elements
figure(4)
subplot(2, 3, 1)
plot(t, semimaj, 'k')
title('Semi-Major Axis')
xlabel('Time [s]')
ylabel('Length Units')
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
ylabel('Degrees)')
subplot(2, 3, 6)
plot(t, anom, 'k')
title('True Anomaly')
xlabel('Time [s]')
ylabel('Degrees')
sgtitle('Classical Orbital Elements of Simulation')

%% 1.iv and 1.v
% Initializing vectors to store the position and semi-major axis error.
% Following that the simulation is run again.
poserrortot = zeros(length(tols), 1);
smaerrortot = zeros(length(tols), 1);
for var = 1:length(tols)
    options = odeset('RelTol',tols(var),'AbsTol',tols(var),'Events', @CompletedPeriodEventFuncHW2P1);
    [~,r,~,~,~] = ode45(@(t,r)RelTwoBodyProp(t,r,mu),tSpan,init,options);

    rvec = [r(end, 1), r(end, 2), r(end, 3)];
    vvec = [r(end, 4), r(end, 5), r(end, 6)];

    specE = 0.5*norm(vvec)^2 - (mu / norm(rvec));
    sma = -mu/(2*specE);

    poserrortot(var) = sqrt((r0(1)-rvec(1))^2 + (r0(2)-rvec(2))^2 + (r0(3)-rvec(3))^2);
    smaerrortot(var) = abs(sma-smai);
end

% Plotting the position and semi-major axis error as functions of the ODE
% tolerance
figure(5)
loglog(tols, poserrortot,'r')
title('Log-Log Plot of Error in Position of 1 Period')
xlabel('ODE Tolerance')
ylabel('Error in Final Positon [LU]')

figure(6)
loglog(tols, smaerrortot,'r')
title('Log-Log Plot of Error in Semi-Major Axis after 1 Period')
xlabel('ODE Tolerance')
ylabel('Error in Semi-Major Axis [LU]')

%% 1.vi
% Running the simulation again to produce semi-major axis change as a
% function of time of flight.
mags = -4:1:4;
tSpan = 1*10.^(mags);
options = odeset('RelTol',tols(end),'AbsTol',tols(end));

smaerrortot = zeros(length(tSpan),1);
for var = 1:length(tSpan)
    [~,r] = ode45(@(t,r)RelTwoBodyProp(t,r,mu),[0 tSpan(var)],init,options);

    rvec = [r(end, 1), r(end, 2), r(end, 3)];
    vvec = [r(end, 4), r(end, 5), r(end, 6)];

    specE = 0.5*norm(vvec)^2 - (mu / norm(rvec));
    sma = -mu/(2*specE);

    smaerrortot(var) = abs(sma-smai);
end

% Plotting the semi-major axis error as a function of time of flight
figure(7)
loglog(tSpan, smaerrortot,'r')
title('Log-Log Plot of Change in Semi-Major Axis as Function of Time of Flight')
xlabel('Time of Flight [TU]')
ylabel('Change in Semi-Major Axis [LU]')

%% 1.vii
% Running the simulation to find the number of steps necessary to
% completely integrate each ODE tolerance
tSpan = [0 1e3];
inttots = zeros(length(tols),1);
for var = 1:length(tols)
    options = odeset('RelTol',tols(var),'AbsTol',tols(var),'Events', @CompletedPeriodEventFuncHW2P1);
    [t,~,~,~,~] = ode45(@(t,r)RelTwoBodyProp(t,r,mu),tSpan,init,options);

    inttots(var) = length(t);
end

figure(8)
loglog(tols, inttots,'r')
title('Log-Log Plot of Integration Steps for Each Tolerance')
xlabel('ODE Tolerance')
ylabel('Integration Steps')