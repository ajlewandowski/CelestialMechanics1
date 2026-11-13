% This file is the event function for ode45 when the propogated orbit
% completes one period
%
% Alec Lewandowski
% 16 Sep 2025

function [val, timeStop, direction] = CompletedPeriodEventFuncHW2P1(t, y)

r0 = [1; 0; 0];

val = max(r0 - y(1:3));

timeStop = 1;
direction = -1;

end