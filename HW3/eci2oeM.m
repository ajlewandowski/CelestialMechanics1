% This function convert the r and v state vectors to the classical orbital
% elements. 
%
% output is semi-major axis, eccentricity, inclination, arg of periapsis, 
% ascending node, mean anomaly in units of degrees
%
% Alec Lewandowski
% 10/6/25

function [oe] = eci2oeM(r, v, mu)

rmag = norm(r);
vmag = norm(v);

ivec = [1, 0, 0];
kvec = [0, 0, 1];

specE = 0.5*vmag^2 - (mu / rmag);
sma = -mu/(2*specE);

hvec = cross(r, v);
h = norm(hvec);

evec = cross(v, hvec) / mu - r / rmag;
ecc = sqrt(1 + (2*specE*h^2)/mu^2);

inc = acos(dot(kvec, hvec)/h)*(180/pi);

nvec = cross(kvec, hvec) / norm(cross(kvec, hvec));

ascn = (2*pi - acos(dot(ivec, nvec)));

argp = (2*pi - acos(dot(nvec, evec)/ecc));

anom = acos(dot(evec, r)/(ecc*rmag));

E = 2*atan(tan(anom/2)*(sqrt((1+ecc))/(1-ecc))^-1);

M = E - ecc*sin(E);

oe = [sma, ecc, rad2deg(inc), rad2deg(argp), rad2deg(ascn), rad2deg(M)];


end