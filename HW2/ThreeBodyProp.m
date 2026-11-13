% This function is the ode function for a three-body problem.
%
% Alec Lewandowski
% 16 Sep 2025

function dydt = ThreeBodyProp(t, y, G, M1, M2, M3)

r1 = y(1:3);
r1dot = y(4:6);

r2 = y(7:9);
r2dot = y(10:12);

r3 = y(13:15);
r3dot = y(16:18);

r12 = r2 - r1;
r13 = r3 - r1;
r21 = r1 - r2;
r23 = r3 - r2;
r31 = r1 - r3;
r32 = r2 - r3;

r1ddot = G * ((M2*r12/norm(r12)^3) + (M3*r13/norm(r13)^3));
r2ddot = G * ((M1*r21/norm(r21)^3) + (M3*r23/norm(r23)^3));
r3ddot = G * ((M1*r31/norm(r31)^3) + (M2*r32/norm(r32)^3));

dydt = [r1dot; r1ddot; r2dot; r2ddot; r3dot; r3ddot];

end