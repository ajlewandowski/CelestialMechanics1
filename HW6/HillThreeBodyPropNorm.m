% This function propogates the STM for two body with radial thrust
%
% Alec Lewandowski
% 11/12/25

function dtrajdt = HillThreeBodyPropNorm(t, traj)

x = traj(1);
y = traj(2);
z = traj(3);
u = traj(4);
v = traj(5);
w = traj(6);
STM = reshape(traj(7:end), [6,6]);

r = sqrt(x^2 + y^2 + z^2);

xddot = 2*v - x/r^3 + 3*x;
yddot = -2*u - y/r^3;
zddot = -z/r^3 - z;

deriv = [u; v; w; xddot; yddot; zddot];


F = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; -(x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.3e1 / 0.2e1) + 3 * x ^ 2 * (x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1) + 3 3 * y * (x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1) * x 3 * z * (x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1) * x 0 0 0; 3 * y * (x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1) * x -(x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.3e1 / 0.2e1) + 3 * y ^ 2 * (x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1) 3 * z * (x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1) * y 0 0 0; 3 * z * (x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1) * x 3 * z * (x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1) * y -(x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.3e1 / 0.2e1) + 3 * z ^ 2 * (x ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1) - 1 0 0 0;];

STMdot = F * STM;

% derivative vector
dtrajdt = zeros(42,1);
dtrajdt(1:6) = deriv;
dtrajdt(7:end) = STMdot(:);