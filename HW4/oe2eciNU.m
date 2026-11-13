% This function converts the orbital elements to the classic cartesian
% coordinates given nu as 6th element
%
% Alec Lewandowski
% 10/6/25

function [r,v] = oe2eciNU(oeVec,mu)

a = oeVec(1);   % Semi-major axis
e = oeVec(2);   % Eccentricity
i = oeVec(3);   % Inclination
aop = oeVec(4); % Argument of periapsis
LAN = oeVec(5); % RAAN (Right Ascending Node)
nu = oeVec(6);  % True anomaly

% Compute the semi-latus rectum
p = a * (1 - e^2);

% Position vector in perifocal coordinates (r_pf)
r_pf = [p * cos(nu) / (1 + e * cos(nu));
        p * sin(nu) / (1 + e * cos(nu));
        0];

% Velocity vector in perifocal coordinates (v_pf)
v_pf = sqrt(mu / p) * [-sin(nu);
                        e + cos(nu);
                        0];

% Define the rotation matrices
RLAN = [cos(LAN), -sin(LAN), 0;
        sin(LAN),  cos(LAN), 0;
        0,         0,        1];

Ri = [1, 0, 0;
      0, cos(i), -sin(i);
      0, sin(i),  cos(i)];

Raop = [cos(aop), -sin(aop), 0;
        sin(aop),  cos(aop), 0;
        0,         0,        1];

% The combined rotation matrix: 
% First rotate by argument of periapsis, then inclination, then RAAN
Q = RLAN * Ri * Raop;

% Convert the position and velocity vectors to ECI frame
r = Q * r_pf;
v = Q * v_pf;