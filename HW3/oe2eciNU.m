% This function converts the orbital elements to the classic cartesian
% coordinates given nu as 6th element
%
% Alec Lewandowski
% 10/6/25

function [r,v] = oe2eciNU(oeVec,mu)

a = oeVec(1);
e = oeVec(2);
i = oeVec(3);
aop = oeVec(4);
LAN = oeVec(5);
nu = oeVec(6);

% Calculating the radius vector in the orbital plane
r = a * (1 - e^2) / (1 + e * cos(nu));
rp = zeros(3,1);
rp(1) = r*cos(nu);
rp(2) = r*sin(nu);

% Calculating the velocity vector in the orbital plane
p = a*(1-e^2);
vp = zeros(3,1);
vp(1) = -sqrt(mu/p)*sin(nu);
vp(2) = sqrt(mu/p)*(e+cos(nu));

% Construct the rotation matrix from PQW to ECI
% using the 3-1-3 Euler rotation sequence
RLAN = [cos(LAN), -sin(LAN), 0; sin(LAN), cos(LAN), 0; 0, 0, 1];
Ri = [1, 0, 0; 0, cos(i), -sin(i); 0, sin(i), cos(i)];
Raop = [cos(aop), -sin(aop), 0; sin(aop), cos(aop), 0; 0, 0, 1];

% The full rotation matrix is the product of the individual rotations
R = RLAN * Ri * Raop;
r = R*rp;
v = R*vp;