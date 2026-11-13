% This function propogates the STM for two body with radial thrust
%
% Alec Lewandowski
% 11/12/25

function dydt = RelTwoBodyThrustProp(t, y, GM, athrust)

r   = y(1:3);
v   = y(4:6);
STM = reshape(y(7:end), [6,6]);
rnorm = norm(r);

aVec = -GM * r / rnorm^3;
aThrust = athrust * r / rnorm;
aTot = aVec + aThrust;

r1 = r(1);
r2 = r(2);
r3 = r(3);

% getting F, see maple file for derivation
F = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r1 ^ 2 - GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r1 ^ 2 + athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.1e1 / 0.2e1) 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r2 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r2 * r1 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r1 0 0 0; 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r2 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r2 * r1 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r2 ^ 2 - GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r2 ^ 2 + athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.1e1 / 0.2e1) 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r2 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r2 0 0 0; 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r1 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r1 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 * r2 - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 * r2 3 * GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.5e1 / 0.2e1) * r3 ^ 2 - GM * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) - athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.3e1 / 0.2e1) * r3 ^ 2 + athrust * (r1 ^ 2 + r2 ^ 2 + r3 ^ 2) ^ (-0.1e1 / 0.2e1) 0 0 0;];

STMdot = F * STM;

% derivative vector
dydt = zeros(42,1);
dydt(1:3) = v;
dydt(4:6) = aTot;
dydt(7:end) = STMdot(:);