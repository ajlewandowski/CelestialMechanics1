% This function runs solution for HW3 Problem 1
%
% Alec Lewandowski
% 10/6/25
clear variables
close all
clc

% Defining the variable for part C
eps = 1e-5;
mu = 1;
tpe = 2;
a = 1.5;
e = 0.7;
i = 70;
aop = -30;
LAN = 160;
M0 = -90;
oeVec = [a, e, deg2rad(i), deg2rad(aop), deg2rad(LAN), deg2rad(M0)];

[JfirstDeriv] = firstDerivJApprox(oeVec, tpe, mu, eps)
[Janalytic] = partialJ(oeVec, tpe, mu)

W = JfirstDeriv - Janalytic

% Computing the performance metric in part iii
Q = 0;
for i = 1:6
    for j = 1:6
        temp = (i*j + i + j)*Janalytic(i,j);
        Q = Q + temp;
    end
end
fprintf('Q = \n')
fprintf('%1.16f \n \n',Q)

% Repeating approximate computation with different epsilon
eps = [1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2];
errorJ = zeros(length(eps),1);
for k = 1:length(eps)
    [JfirstDeriv] = firstDerivJApprox(oeVec, tpe, mu, eps(k));

    Wtemp = Janalytic - JfirstDeriv;

    errorJ(k) = norm(Wtemp);
end

figure(1)
loglog(eps,errorJ)
xlabel('Epsilon Values')
ylabel('W Norm Error')
title('Error in Approximate versus Analytic J')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART D
tpe = [1:0.1:12];
LB = zeros(6,6,length(tpe));
temp = [zeros(3), eye(3); -eye(3), zeros(3)];
for k = 1:length(tpe)
    J = partialJ(oeVec, tpe(k), mu);
   
    LB(:,:,k) = J'*temp*J;
end

figure(2)
for i = 1:6
    for j = 1:6
        if i>=j
            data = squeeze(LB(i,j,:));
            exactmean(i,j) = mean(data);
            exactspread(i,j) = max(data) - min(data);
            subplot(6,6,j+(i-1)*6);
            plot(tpe, data - mean(data), 'r')
            title([num2str(exactmean(i,j), '%.2e') ' , ' num2str(exactspread(i,j), '%.2e')])
            xlabel('Time [TU]')
            ylabel(['\Delta LB[' num2str(i) ',' num2str(j) ']'])
            sgtitle('(Mean, Spread)')
        end
    end
end

exactmean

exactspread


tpe = [1:0.1:12];
LB = zeros(6,6,length(tpe));
temp = [zeros(3), eye(3); -eye(3), zeros(3)];
for k = 1:length(tpe)
    J = firstDerivJApprox(oeVec, tpe(k), mu, 1e-8);
   
    LB(:,:,k) = J'*temp*J;
end

figure(3)
for i = 1:6
    for j = 1:6
        if i>=j
            data = squeeze(LB(i,j,:));
            exactmean(i,j) = mean(data);
            exactspread(i,j) = max(data) - min(data);
            subplot(6,6,j+(i-1)*6);
            plot(tpe, data - mean(data), 'r')
            title([num2str(exactmean(i,j), '%.2e') ' , ' num2str(exactspread(i,j), '%.2e')])
            xlabel('Time [TU]')
            ylabel(['\Delta LB[' num2str(i) ',' num2str(j) ']'])
            sgtitle('(Mean, Spread)')
        end
    end
end

exactmean

exactspread