% This function is the ode function for two-body field with perturbation
% acceleration a = Sr + TTheta + WH
%
% Alec Lewandowski
% 10/22/25

function [dydt] = TwoBodyaPert(t, y, mu, aPert)

r = y(1:3);
v = y(4:6);

r_norm = norm(r);

h = cross(r, v);
h_norm = norm(h);

r_hat = r / r_norm;    
h_hat = h / h_norm;    
t_hat = cross(h_hat, r_hat); 

S = aPert(1);
T = aPert(2);
W = aPert(3);

aPert_RTN = S * r_hat + T * t_hat + W * h_hat;  

aPert_ECI = aPert_RTN;

a_grav = -mu * r / r_norm^3; 

a_total = a_grav + aPert_ECI;  

dydt = zeros(6,1);

dydt(1:3) = v;
dydt(4:6) = a_total;