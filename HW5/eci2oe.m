% This function convert the r and v state vectors to the classical orbital
% elements. 
%
% output is semi-major axis, eccentricity, inclination, arg of periapsis, 
% ascending node, true anomaly in units of degrees
%
% Alec Lewandowski
% 10 Nov 2025

function [oe] = eci2oe(r, v, mu)

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

ascn = (2*pi - acos(dot(ivec, nvec)))*(180/pi);

argp = (2*pi - acos(dot(nvec, evec)/ecc))*(180/pi);

anom = acos(dot(evec, r)/(ecc*rmag))*(180/pi);

oe = [sma, ecc, inc, argp, ascn, anom];

end