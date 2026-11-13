% % This function is the ode function for a two-body solar radiation
% pressure problem using the LPEs to numerically integrate.
%
% Alec Lewandowski
% 10/6/25

function [dydt] = TwoBodySRPinLPE(t, y, gammaSRP, mu, dhat)

a = y(1);
e = y(2);
i = y(3);
aop = y(4);
LAN = y(5);
M = y(6);

eps = 1e-10;
e = max(e,eps);
sin_i = sin(i);
if abs(sin_i) < eps
    sin_i = eps;
end
sqrt_1_e2 = max(sqrt(1-e^2),eps);
eccvar = max(1-e^2,eps);


ehat = [cos(LAN)*cos(aop) - sin(LAN)*sin(aop)*cos(i);
    sin(LAN)*cos(aop) + cos(LAN)*sin(aop)*cos(i);
    sin(aop)*sin_i];

n = sqrt(mu/a^3);

dehatdi = [sin(LAN)*sin(aop)*sin_i; -cos(LAN)*sin(aop)*sin_i; sin(aop)*cos(i)];
dehatdaop = [-cos(LAN)*sin(aop) - sin(LAN)*cos(aop)*cos(i); -sin(LAN)*sin(aop) + cos(LAN)*cos(aop)*cos(i); cos(aop)*sin_i];
dehatdLAN = [-sin(LAN)*cos(aop) - cos(LAN)*sin(aop)*cos(i); cos(LAN)*cos(aop) - sin(LAN)*sin(aop)*cos(i); 0];

dRda = -1.5*gammaSRP*dot(dhat,ehat);
dRde = -1.5*gammaSRP*sum(dhat);
dRdi = -1.5*gammaSRP*a*dot(dhat,dehatdi);
dRdaop = -1.5*gammaSRP*a*dot(dhat,dehatdaop);
dRdLAN = -1.5*gammaSRP*a*dot(dhat,dehatdLAN);
dRdM = 0;

dadt = 2*dRdM/(n*a);
dedt = -(sqrt_1_e2/e*n*a^2)*dRdaop + (eccvar)*dRdM/(n*e*a^2);
didt = (1/(n*sin(i)*sqrt_1_e2*a^2))*(cos(i)*dRdaop - dRdLAN);
daopdt = dRde*sqrt(eccvar)/n*e*a^2 - dRdi*cos(i)/sin_i*(sqrt_1_e2*n*a^2);
dLANdt = dRdi/(sin_i*sqrt_1_e2*n*a^2);
dMdt = n + (-dRde*(eccvar)/n*e*a^2 - dRda*2/n*a);

dydt = [dadt; dedt; didt; daopdt; dLANdt; dMdt];