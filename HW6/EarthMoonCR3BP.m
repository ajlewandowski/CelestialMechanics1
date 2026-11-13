% This function propgates the Earth-Moon CR3BP

function deriv = EarthMoonCR3BP(t, traj, mu)

x = traj(1);
y = traj(2);
z = traj(3);
xdot = traj(4);
ydot = traj(5);
zdot = traj(6);

r1 = sqrt((x+mu)^2 + y^2 + z^2);
r2 = sqrt((x-1+mu)^2 + y^2 + z^2);

var1 = -(1-mu)/r1^3;
var2 = -mu/r2^3;

xddot = 2*ydot + x + var1*(x+mu) + var2*(x-1+mu);
yddot = -2*xdot + y + var1*y + var2*y;
zddot = var1*z + var2*z;

deriv = [xdot; ydot; zdot; xddot; yddot; zddot];