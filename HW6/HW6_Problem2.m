% This file runs solution for HW6 Problem 2
%
% Alec Lewandowski
% 12/1/25
clear variables
close all
clc
format long

fig = 1;

%% Problem 2.B
% Defining the Lagrange Points
L1 = [0.8369180073; 0; 0];
L2 = [1.155679913; 0; 0];
L3 = [-1.005062402; 0; 0];
L4 = [0.48785; 0.8660254038; 0];
L5 = [0.48785; -0.8660254038; 0];
mu = 0.01215;

L1eigs = eigs(getAMatCR3BP(L1(1), L1(2), 0, mu))
L2eigs = eigs(getAMatCR3BP(L2(1), L2(2), 0, mu))
L3eigs = eigs(getAMatCR3BP(L3(1), L3(2), 0, mu))
L4eigs = eigs(getAMatCR3BP(L4(1), L4(2), 0, mu))
L5eigs = eigs(getAMatCR3BP(L5(1), L5(2), 0, mu))

fprintf('Therefore because L1-L3 have real components they are unstable saddles. \n')
fprintf('L4-L5 are stable because all eigs are imaginary only.\n\n')

%% Problem 2.C
tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',tol);
eps = 1e-3*[1;1;0];
velos = [0; 0; 0];

% Testing L1
IC = [L1+eps; velos];
[tvec1, traj1] = ode45(@(t,traj)EarthMoonCR3BP(t,traj,mu),[0 100],IC,options);

% Testing L2
IC = [L2+eps; velos];
[tvec2, traj2] = ode45(@(t,traj)EarthMoonCR3BP(t,traj,mu),[0 100],IC,options);

% Testing L3
IC = [L3+eps; velos];
[tvec3, traj3] = ode45(@(t,traj)EarthMoonCR3BP(t,traj,mu),[0 100],IC,options);

% Testing L4
IC = [L4+eps; velos];
[tvec4, traj4] = ode45(@(t,traj)EarthMoonCR3BP(t,traj,mu),[0 100],IC,options);

% Testing L5
IC = [L5+eps; velos];
[tvec5, traj5] = ode45(@(t,traj)EarthMoonCR3BP(t,traj,mu),[0 100],IC,options);

figure(fig)
hold on
plot(traj1(:,1),traj1(:,2),'r')
plot(traj1(1,1),traj1(1,2),'bo')
plot(traj1(end,1),traj1(end,2),'ks')
hold off
legend('Traj','Start','End')
xlabel('X [LU]')
ylabel('Y [LU]')
title('L1 Stability Test')
fig = fig + 1;

figure(fig)
hold on
plot(traj2(:,1),traj2(:,2),'r')
plot(traj2(1,1),traj2(1,2),'bo')
plot(traj2(end,1),traj2(end,2),'ks')
hold off
legend('L2 Traj','Start','End')
xlabel('X [LU]')
ylabel('Y [LU]')
title('L2 Stability Test')
fig = fig + 1;

figure(fig)
hold on
plot(traj3(:,1),traj3(:,2),'r')
plot(traj3(1,1),traj3(1,2),'bo')
plot(traj3(end,1),traj3(end,2),'ks')
hold off
legend('L3 Traj','Start','End')
xlabel('X [LU]')
title('L3 Stability Test')
fig = fig + 1;ylabel('Y [LU]')

figure(fig)
hold on
plot(traj4(:,1),traj4(:,2),'r')
plot(traj4(1,1),traj4(1,2),'bo')
plot(traj4(end,1),traj4(end,2),'ks')
hold off
legend('L4 Traj','Start','End')
xlabel('X [LU]')
ylabel('Y [LU]')
title('L4 Stability Test')
fig = fig + 1;

figure(fig)
hold on
plot(traj5(:,1),traj5(:,2),'r')
plot(traj5(1,1),traj5(1,2),'bo')
plot(traj5(end,1),traj5(end,2),'ks')
hold off
legend('L5 Traj','Start','End')
xlabel('X [LU]')
ylabel('Y [LU]')
title('L5 Stability Test')
fig = fig + 1;

%% Problem 2.D
% Testing L5
IC = [eps(1); eps(2); 0; 0];
A = getAMatCR3BP(L5(1), L5(2), 0, mu);
[~, lintraj51] = ode45(@(t,ds)A*ds,tvec5,IC,options);

figure(fig)
hold on
plot(lintraj51(:,1)+L5(1),lintraj51(:,2)+L5(2),'b')
plot(traj5(:,1),traj5(:,2),'r--')
plot(traj5(1,1),traj5(1,2),'bo')
plot(traj5(end,1),traj5(end,2),'ks')
hold off
legend('L5 Lin Traj','L5 Traj','Start','End','location','northwest')
xlabel('X [LU]')
ylabel('Y [LU]')
title('L5 Linearized vs Nonlinear Test')
fig = fig + 1;

normDiff = zeros(length(tvec5),1);
for i = 1:length(normDiff)
    normDiff(i) = norm(lintraj51(i,1:2)'+L5(1:2) - traj5(i,1:2)');
end

figure(fig)
semilogy(tvec5, normDiff, 'r')
xlabel('Time [TU]')
ylabel('Position Difference Norm')
title('Linear vs Nonlinear Norm Difference')
fig = fig + 1;

%% Problem 2.E
% Repeating C and D with multiple epsilons, going to create graph that
% shows the error in

epsvec = [-1e-3; -1e-4; -1e-5; -1e-7; 0; 1e-7; 1e-5; 1e-4; 1e-3];
errvec = zeros(length(epsvec),1);

for i = 1:length(epsvec)
    IC = [L5+epsvec(i); velos];
    [tvec5, traj5] = ode45(@(t,traj)EarthMoonCR3BP(t,traj,mu),[0 100],IC,options);

    IC = [epsvec(i); epsvec(i); 0; 0];
    A = getAMatCR3BP(L5(1), L5(2), 0, mu);
    [~, lintraj] = ode45(@(t,ds)A*ds,[0 100],IC,options);

    errvec(i) = norm([lintraj(end,1)+L5(1); lintraj(end,2)+L5(2)] - [traj5(end,1); traj5(end,2)]);
end

figure(fig)
semilogy(epsvec, errvec, 'r')
xlabel('Eps Value')
ylabel('Error in Final Pos')
title('Epsilon Comparison')
fig = fig + 1;

%% Problem 2.F
% Proving mu star
muvec = [0:0.00001:0.05];
dataMat = zeros(11,length(muvec));
maxreal4 = zeros(length(muvec),1);
maxreal5 = zeros(length(muvec),1);

for i = 1:length(muvec)
    mu = muvec(i);
    A4 = getAMatCR3BP(-mu+0.5, sqrt(3)/2, 0, muvec(i));
    A5 = getAMatCR3BP(-mu+0.5, -sqrt(3)/2, 0, muvec(i));
    
    eigs4 = eigs(A4);
    eigs5 = eigs(A5);

    dataMat(1,i) = mu;
    
    maxreal4(i) = round(max(real(eigs4)),8);
    maxreal5(i) = round(max(real(eigs5)),8);

    if maxreal4(i) <= 0
        dataMat(2,i) = -1;
    else
        dataMat(2,i) = 1;
    end

    dataMat(3:6) = eigs4;

    if maxreal5(i) <= 0
        dataMat(7,i) = -1;
    else
        dataMat(7,i) = 1;
    end

    dataMat(8:11) = eigs5;
end

signVec4 = sign(dataMat(2,:));
signVec5 = sign(dataMat(7,:));

sD4 = diff(signVec4);
sD5 = diff(signVec5);

tranInd4 = find(sD4 == 2);
tranInd5 = find(sD5 == 2);

mutransition4 = dataMat(1,tranInd4 + 1)
mutransition5 = dataMat(1,tranInd5 + 1)

mustar = (1-sqrt(23/27))/2

% Before case
IC = [eps(1); eps(2); 0; 0];
A = getAMatCR3BP(L5(1), L5(2), 0, mustar-0.01);
[~, lintraj5Before] = ode45(@(t,ds)A*ds,[0 100],IC,options);

% After case
IC = [eps(1); eps(2); 0; 0];
A = getAMatCR3BP(L5(1), L5(2), 0, mustar+0.01);
[~, lintraj5After] = ode45(@(t,ds)A*ds,[0 20],IC,options);

figure(fig)
hold on
plot(lintraj5Before(:,1)+L5(1),lintraj5Before(:,2)+L5(2),'b')
plot(lintraj5After(:,1)+L5(1),lintraj5After(:,2)+L5(2),'r')
plot(L5(1),L5(2),'ko')
hold off
legend('L5 Before','L5 After','Start','location','northwest')
xlabel('X [LU]')
ylabel('Y [LU]')
title('L5 MuStar Test')
fig = fig + 1;
