% This function runs solution for HW4 Problem 4
%
% Alec Lewandowski
% 10/22/25
clear variables
close all
clc
format long

% Adding Mice path and grabbing those data files
addpath('C:\Users\alecl\OneDrive\Desktop\CM1\HW4\mice\src\mice')
addpath('C:\Users\alecl\OneDrive\Desktop\CM1\HW4\mice\lib')

cspice_furnsh('de421.bsp'); 
cspice_furnsh('pck00010.tpc'); 
cspice_furnsh('naif0011.tls');

% Defining the parameters given for the problem
t0 = 'Oct 14, 2020 11:45 AM CST';
et0 = cspice_str2et(t0);
etf = et0 + 3*30*24*3600;
fig = 1;

%% 4.i
fprintf('t0 = %g seconds past J2000 \n \n', et0);

%% 4.ii
tSpan = linspace(et0, etf, 2000);
relPos = zeros(length(tSpan),3);
relX = zeros(length(tSpan),1);
relY = zeros(length(tSpan),1);
relZ = zeros(length(tSpan),1);

% Iterate through to get position for each time step
for i = 1:length(tSpan)
    [pos] = mice_spkezr('EARTH', tSpan(i), 'ECLIPJ2000', 'NONE', 'MOON');
    relPos(i, :) = pos.state(1:3)';
    relX(i) = pos.state(1);
    relY(i) = pos.state(2);
    relZ(i) = pos.state(3);
end

figure(fig)
subplot(3,1,1)
hold on
plot(relX, relY)
plot(relX(1), relY(1), '*')
plot(relX(end), relY(end), 'o')
legend('Traj','Start Pt','End Pt')
hold off
xlabel('X Position (km)')
ylabel('Y Position (km)')
subplot(3,1,2)
hold on
plot(relX, relZ)
plot(relX(1), relZ(1), '*')
plot(relX(end), relZ(end), 'o')
legend('Traj','Start Pt','End Pt')
hold off
xlabel('X Position (km)')
ylabel('Z Position (km)')
subplot(3,1,3)
hold on
plot(relY, relZ)
plot(relY(1), relZ(1), '*')
plot(relY(end), relZ(end), 'o')
legend('Traj','Start Pt','End Pt')
hold off
xlabel('Y Position (km)')
ylabel('Z Position (km)')
sgtitle('Plot of the Earth Position rel to Moon in Inertial')
fig = fig + 1;

%% 4.iii
relPos = zeros(length(tSpan),3);
relX = zeros(length(tSpan),1);
relY = zeros(length(tSpan),1);
relZ = zeros(length(tSpan),1);

% Iterate through to get position for each time step
for i = 1:length(tSpan)
    [pos] = mice_spkezr('EARTH', tSpan(i), 'IAU_MOON', 'NONE', 'MOON');
    relPos(i, :) = pos.state(1:3)';
    relX(i) = pos.state(1);
    relY(i) = pos.state(2);
    relZ(i) = pos.state(3);
end

figure(fig)
subplot(3,1,1)
hold on
plot(relX, relY)
plot(relX(1), relY(1), '*')
plot(relX(end), relY(end), 'o')
legend('Traj','Start Pt','End Pt')
hold off
xlabel('X Position (km)')
ylabel('Y Position (km)')
subplot(3,1,2)
hold on
plot(relX, relZ)
plot(relX(1), relZ(1), '*')
plot(relX(end), relZ(end), 'o')
legend('Traj','Start Pt','End Pt')
hold off
xlabel('X Position (km)')
ylabel('Z Position (km)')
subplot(3,1,3)
hold on
plot(relY, relZ)
plot(relY(1), relZ(1), '*')
plot(relY(end), relZ(end), 'o')
legend('Traj','Start Pt','End Pt')
hold off
xlabel('Y Position (km)')
ylabel('Z Position (km)')
sgtitle('Plot of the Earth Position rel to Moon in IAU Moon')
fig = fig + 1;

%% 4.iv
rMoon = 1738;
lon = zeros(length(tSpan),1);
lat = zeros(length(tSpan),1);
dist = zeros(length(tSpan),1);

% Mice gives lon\lat in earth frame, not moon, need to transform
for i = 1:length(tSpan)
    [radius, a, b] = cspice_reclat([relX(i); relY(i); relZ(i)]);
    lon(i) = rad2deg(a);
    lat(i) = rad2deg(b);
    dist(i) = radius;
end

figure(fig)
hold on
plot(lon, lat)
plot(lon(1), lat(1), '*')
plot(lon(end), lat(end), 'o')
legend('Traj','Start Pt','End Pt')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
title('Ground Track of Earth in Moon Frame')
axis equal
fig = fig + 1;

figure(fig)
plot(tSpan/86400, dist)
xlabel('Days since J2000')
ylabel('Distance to Earth (km)')
title('Distance from Moon to Earth')
fig = fig + 1;

minDist = min(dist);
maxDist = max(dist);
avgDist = mean(dist);
fprintf('Min dist = %.2f km, Max dist = %.2f km, Avg dist = %.2f km\n \n', minDist, maxDist, avgDist);

%% Part v
etf = et0 + 18*365*24*3600;
tSpan = linspace(et0, etf, 4000);

poleVec = zeros(length(tSpan),3);

for i = 1:length(tSpan)
    temp = cspice_pxform('IAU_MOON','ECLIPJ2000',tSpan(i));
    poleVec(i,:) = (temp * [0;0;1])';
end

posPole = poleVec * rMoon;
poleX = posPole(:,1);
poleY = posPole(:,2);
poleZ = posPole(:,3);

figure(fig)
hold on
plot3(poleX, poleY, poleZ)
plot3(poleX(1), poleY(1), poleZ(1),'*')
plot3(poleX(end), poleY(end), poleZ(end),'o')
hold off
legend('Traj','Start Pt','End Pt')
xlabel('X Position (km)')
ylabel('Y Position (km)')
zlabel('Z Position (km)')
title('Cartesian Evolution of the Moon North Pole')
fig = fig + 1;


figure(fig)
subplot(3,1,1)
hold on
plot(poleX, poleY)
plot(poleX(1), poleY(1), '*')
plot(poleX(end), poleY(end), 'o')
legend('Traj','Start Pt','End Pt')
hold off
xlabel('X Position (km)')
ylabel('Y Position (km)')
subplot(3,1,2)
hold on
plot(poleX, poleZ)
plot(poleX(1), poleZ(1), '*')
plot(poleX(end), poleZ(end), 'o')
legend('Traj','Start Pt','End Pt')
hold off
xlabel('X Position (km)')
ylabel('Z Position (km)')
subplot(3,1,3)
hold on
plot(poleY, poleZ)
plot(poleY(1), poleZ(1), '*')
plot(poleY(end), poleZ(end), 'o')
legend('Traj','Start Pt','End Pt')
hold off
xlabel('Y Position (km)')
ylabel('Z Position (km)')
sgtitle('Cartesian Evolution of the Moon North Pole')
fig = fig + 1;

temp = poleVec * [0;0;1];
tempAng = acosd(temp);
tempOblate = mean(tempAng);
fprintf('Mean Oblateness Angle = %.2f degrees\n \n', tempOblate);