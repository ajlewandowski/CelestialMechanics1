% This function is the ode function for two-body field with perturbation
% acceleration a = Sr + TTheta + WH
%
% Alec Lewandowski
% 10/22/25

function [dOE] = TwoBodyaPertGauss(t, y, mu, aPert)

% Orbital elements
a = y(1);    % semi-major axis
e = y(2);    % eccentricity
i = y(3);    % inclination
aop = y(4);  % argument of periapsis
LAN = y(5);  % RAAN (Right Ascending Node)
nu = y(6);   % true anomaly

% Perturbation parameters
S = aPert(1);
T = aPert(2);
W = aPert(3);

% Semi-latus rectum
p = a * (1 - e^2);

% Radius
r = p / (1 + e * cos(nu));

% Specific angular momentum
h = sqrt(mu * p);

n = sqrt(mu/a^3);

% Acceleration perturbations in RTN frame (S, T, W)
u = nu + aop;  % argument of latitude

% Compute derivatives of orbital elements:

% Semi-major axis (da/dt)
dadt = (2 * a^2 / h) * (e * sin(nu) * S + (p / r) * T);

% Eccentricity (de/dt)
dedt = (1 / h) * (p * sin(nu) * S + ((p + r) * cos(nu) + r * e) * T);

% Inclination (di/dt)
didt = (r * cos(u) / h) * W;

% RAAN (dΩ/dt)
dOmegadt = (r * sin(u) / (h * sin(i))) * W;

% Argument of periapsis (dω/dt)
domegadt = (1 / (h * e)) * (-p * cos(nu) * S + (p + r) * sin(nu) * T) - (r * sin(u) * cos(i) / (h * sin(i))) * W;

% True anomaly (dν/dt)
dnudt = (h / r^2) + (1 / (h * e)) * (p * cos(nu) * S - (p + r) * sin(nu) * T);

% Return the derivatives of the orbital elements
dOE = [dadt; dedt; didt; domegadt; dOmegadt; dnudt];

 