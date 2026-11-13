% Celestial Mechanics HW1
% Alec Lewandowski
% September 4th, 2025
clear variables
close all
clc

%% PROBLEM 3ii

% Defining vector values for mesh and constants.
vec = [10:10:500];
G = 6.67e-11;
M = 5.972e24;
m = 4.2e5;

% Iterating through vector "mesh" values for case (a), the orbiting mass on
% the inside of the shell.
r = 600e3;
R = 6378e3;
ftot = zeros(length(vec),1);
for i = 1:length(vec)
    th = linspace(0, pi, vec(i)+1);

    fnew = (G*M*m*0.5*sin(th).*(r-R*cos(th)))./((r^2 + R^2 - 2*r*R*cos(th)).^1.5);

    fmid = cumtrapz(fnew, th);

    ftot(i) = fmid(end);
end


figure(1)
hold on
plot(vec, ftot, 'r-')
yline(0)
hold off
xlabel('Mesh Size')
ylabel('Total Force Value [N]')
legend('Mascon Model','Analytical (F = 0)')
title('Mesh Model of an Obriting Body Inside a Spherical Shell')




% Iterating through vector "mesh" values for case (a), the orbiting mass on
% the inside of the shell.
r = 6778e3;
ftot = zeros(length(vec),1);
for k = 1:length(vec)
    th = linspace(0, pi, vec(k)+1);

    fnew = (G*M*m*0.5*sin(th).*(r-R*cos(th)))./((r^2 + R^2 - 2*r*R*cos(th)).^1.5);

    fmid = cumtrapz(fnew, th);

    ftot(k) = fmid(end);
end

analytical = -G*m*M/r^2;


figure(2)
hold on
plot(vec, ftot, 'r-')
yline(analytical)
hold off
xlabel('Mesh Size')
ylabel('Total Force Value [N]')
legend('Mascon Model','Analytical')
title('Mesh Model of an Obriting Body Outside a Spherical Shell')