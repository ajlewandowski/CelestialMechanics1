% This function is the ode function for a two-body solar radiation pressure problem.
%
% Alec Lewandowski
% 10/6/25

function [dydt] = TwoBodySRP(t, y, gammaSRP, mu, dhat)

rx = y(1);
ry = y(2);
rz = y(3);
rxdot = y(4);
rydot = y(5);
rzdot = y(6);

rmag = sqrt(rx^2 + ry^2 + rz^2);

rddot = (-mu/rmag^3).*[rx; ry; rz] + gammaSRP*dhat;

dydt = [rxdot; rydot; rzdot; rddot];