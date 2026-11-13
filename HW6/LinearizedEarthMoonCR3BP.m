% This function propgates the Linearized Earth-Moon CR3BP

function deriv = LinearizedEarthMoonCR3BP(t, traj, mu)

x = traj(1);
y = traj(2);
xdot = traj(3);
ydot = traj(4);

Amat = getAMatCR3BP(x, y, 0, mu);

deriv = Amat*[x; y; xdot; ydot];