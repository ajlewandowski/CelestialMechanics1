% This function is the ode function for a two-body problem.
%
% Alec Lewandowski
% 16 Sep 2025

function dydt = RelTwoBodyProp(t, y, mu)

rx = y(1);
ry = y(2);
rz = y(3);
rxdot = y(4);
rydot = y(5);
rzdot = y(6);

rmag = sqrt(rx^2 + ry^2 + rz^2);

rddot = (-mu/rmag^3).*[rx; ry; rz];

dydt = [rxdot; rydot; rzdot; rddot];

end