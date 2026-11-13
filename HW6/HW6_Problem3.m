% This file runs solution for HW6 Problem 3
%
% Alec Lewandowski
% 12/1/25
clear variables
close all
clc
format long

fig = 1;

mu = 0.01215;
Cvec = [3, 3.17, 3.18, 3.19, 3.2];

%% Problem 3.i
x = linspace(-1.5,1.5,500);
y = linspace(-1.5,1.5,500);
[X,Y] = meshgrid(x,y);
V = zeros(length(x),length(y));

for i = 1:length(Cvec)
    C = Cvec(i);
    for j = 1:length(x)
        for k = 1:length(y)
            r1 = sqrt((x(j)+mu)^2 + y(k)^2);
            r2 = sqrt((x(j)-1+mu)^2 + y(k)^2);
            V(j,k) = x(j)^2 + y(k)^2 + 2*(1-mu)/r1 + 2*mu/r2 - C;
            if V(j,k) < 0
                V(j,k) = NaN;
            end
        end
    end
    figure(fig)
    surf(Y, X, log10(V),'EdgeColor','none');
    view(0,90)
    title(['Log10 V2 Surface for C = ', num2str(C)]);
    xlabel('X [LU]');
    ylabel('Y [LU]');
    colorbar;
    fig = fig + 1;
end

%% Problem 3.ii

x = linspace(-1.5,1.5,500);
y = linspace(-1.5,1.5,500);
[X,Y] = meshgrid(x,y);
V = zeros(length(x),length(y));

for i = 1:length(Cvec)
    C = Cvec(i);
    for j = 1:length(x)
        for k = 1:length(y)
            r1 = sqrt((x(j)+mu)^2 + y(k)^2);
            r2 = sqrt((x(j)-1+mu)^2 + y(k)^2);
            V(j,k) = x(j)^2 + y(k)^2 + 2*(1-mu)/r1 + 2*mu/r2 + 2*0.01*x(j) - C;
            if V(j,k) < 0
                V(j,k) = NaN;
            end
        end
    end
    figure(fig)
    surf(Y, X, log10(V),'EdgeColor','none');
    view(0,90)
    title(['Log10 V2 Surface with Const Thrust for C = ', num2str(C)]);
    xlabel('X [LU]');
    ylabel('Y [LU]');
    colorbar;
    fig = fig + 1;
end