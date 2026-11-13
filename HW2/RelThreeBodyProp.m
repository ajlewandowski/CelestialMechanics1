% This function is the ode function for a three-body problem.
%
% Alec Lewandowski
% 16 Sep 2025

function dydt = RelThreeBodyProp(t, y, G, M1, M2, M3)

mu1 = G*M1;
mu2 = G*M2;
mu3 = G*M3;

r1 = y(1:3);
r1dot = y(4:6);

r2 = y(7:9);
r2dot = y(10:12);

r3 = y(13:15);
r3dot = y(16:18);

r12 = y(19:21);
r12dot = y(22:24);

r13 = r3 - r1;
r21 = -r12;
r23 = r13 - r12;
r31 = -r13;
r32 = -r23;

r1ddot = G * ((M2/norm(r12)^3)*r12 + (M3/norm(r13)^3)*r13);
r2ddot = G * ((M1/norm(r21)^3)*r21 + (M3/norm(r23)^3)*r23);
r3ddot = G * ((M1/norm(r31)^3)*r31 + (M2/norm(r32)^3)*r32);
r12ddot = -(mu1 + mu2)*(r12/norm(r12)^3) + mu3*(r23/norm(r23)^3 - r13/norm(r13)^3);

dydt = [r1dot; r1ddot; r2dot; r2ddot; r3dot; r3ddot; r12dot; r12ddot];

end