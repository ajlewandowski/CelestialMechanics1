% This function takes in the x, y, z, and mu values for a planar CR3BP and 
% produces A matrix for stability to make main code cleaner.

function [A] = getAMatCR3BP(x, y, z, mu)

Amat = [0 0 1 0; 
    0 0 0 1; 
    0.1e1 + 0.3e1 / 0.4e1 * (1 - mu) * (((x + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1)) * ((2 * x + 2 * mu) ^ 2) - ((1 - mu) * ((x + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.3e1 / 0.2e1)) + 0.3e1 / 0.4e1 * mu * (((x - 1 + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1)) * ((2 * x - 2 + 2 * mu) ^ 2) - (mu * ((x - 1 + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.3e1 / 0.2e1)) 0.3e1 / 0.2e1 * (1 - mu) * (((x + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1)) * (2 * x + 2 * mu) * y + 0.3e1 / 0.2e1 * mu * (((x - 1 + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1)) * (2 * x - 2 + 2 * mu) * y 0 2; 0.3e1 / 0.2e1 * (1 - mu) * (((x + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1)) * (2 * x + 2 * mu) * y + 0.3e1 / 0.2e1 * mu * (((x - 1 + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1)) * (2 * x - 2 + 2 * mu) * y 1 + 3 * (1 - mu) * ((x + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1) * y ^ 2 - (1 - mu) * ((x + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.3e1 / 0.2e1) + 3 * mu * ((x - 1 + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.5e1 / 0.2e1) * y ^ 2 - mu * ((x - 1 + mu) ^ 2 + y ^ 2 + z ^ 2) ^ (-0.3e1 / 0.2e1) -2 0;];
A = Amat;