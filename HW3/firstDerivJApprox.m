% This function finds a first order approx of the jacobian of the
% transformation.
%
% Alec Lewandowski
% 10/6/25

function [J] = firstDerivJApprox(oeVec, tpe, mu, eps)

% Initializing the J matrix
J = zeros(6, 6);
[mJ, ~] = size(J);

oeVec(2) = min(max(oeVec(2), 0), 1-eps);

for k = 1:mJ
    deltaOE = zeros(1,mJ);
    deltaOE(k) = eps;

    % Plugging in the plus part of first deriv
    [rPlusE, vPlusE] = oe2eci(oeVec + deltaOE, tpe, mu);
    xPlus = [rPlusE; vPlusE];

    % Plugging in the minus part of first deriv
    [rMinusE, vMinusE] = oe2eci(oeVec - deltaOE, tpe, mu);
    xMinus = [rMinusE; vMinusE];

    % First deriv approx into J matrix
    J(:,k) = (xPlus - xMinus) / (2*eps);
end
