% This function runs solution for HW5 Problem 2
%
% Alec Lewandowski
% 11/12/25
clear variables
close all
clc
format long

% Defining the ICs
ri = [1; 0; 0];
rf = [0; 1; 0];
v0 = [0; 1; 0];
STM1 = eye(6);
TOF = (2*pi)*(5/4);
GM = 1;
athrust = 5e-2;
fig = 1;
TR = 0.01;
maxIter = 100;
costNow = 1;
i = 1;

tol = 1e-12;
tolConv = 1e-8;
options = odeset('RelTol',tol,'AbsTol',tol);

while costNow > tolConv
    [t1, y1] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM,athrust),[0 TOF],[ri; v0; STM1(:)],options);
    PhiFinal = reshape(y1(end,7:end),[6,6]);
    B = PhiFinal(1:3,4:6);
    rfnow = y1(end,1:3)';
    rmiss = rfnow - rf;
    rfPosX(i) = rfnow(1);
    rfPosY(i) = rfnow(2);

    cost(i) = norm(rmiss)/norm(rf);
    costNow = cost(i);
    costInd(i) = i;
    
    delV0 = B\rmiss;
    dV0 = min([norm(delV0),TR])*delV0/norm(delV0);
    v0 = v0 - dV0;
    i = i+1;
end

%% Results
vFinal = y1(end,4:6)
totalDelV = norm(v0 - [0; 1; 0])
numIter = i
convPos = y1(end,1:3)
convVel = y1(end,4:6)

figure(fig)
hold on
plot(y1(:,1),y1(:,2),'r')
plot(y1(1,1),y1(1,2),'b*')
plot(y1(end,1),y1(end,2),'bs')
plot(rfPosX,rfPosY,'-ko')
plot(rfPosX(1),rfPosY(1),'g*')
plot(rfPosX(end),rfPosY(end),'gs')
hold off
axis equal
grid on
xlabel('X-Pos [LU]')
ylabel('Y-Pos [LU]')
legend('Trajectory','Start Pt','End Pt','Final Positions','Final Pos 1','Final Pos End')
title('Problem 2: Perturbed Lambert Shooting')
fig = fig + 1;

figure(fig)
plot(costInd,cost,'r')
xlabel('Iterations')
ylabel('Norm(rmiss)/Norm(r0)')
title('Cost as Function of Time')
fig = fig + 1;






%% Experimenting with different TR
TRs = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005];
iters = zeros(length(TRs),1);
for j = 1:length(TRs)
    i = 1;ri = [1; 0; 0];
    rf = [0; 1; 0];
    v0 = [0; 1; 0];
    STM1 = eye(6);
    TOF = (2*pi)*(5/4);
    GM = 1;
    athrust = 5e-2;
    costNow = 1;
    while costNow > tolConv
        [t1, y1] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM,athrust),[0 TOF],[ri; v0; STM1(:)],options);
        PhiFinal = reshape(y1(end,7:end),[6,6]);
        B = PhiFinal(1:3,4:6);
        rfnow = y1(end,1:3)';
        rmiss = rfnow - rf;
        rfPosX(i) = rfnow(1);
        rfPosY(i) = rfnow(2);
    
        cost(i) = norm(rmiss)/norm(rf);
        costNow = cost(i);
        costInd(i) = i;
        
        delV0 = B\rmiss;
        dV0 = min([norm(delV0),TRs(j)])*delV0/norm(delV0);
        v0 = v0 - dV0;
        i = i+1;
    end
    if j == 1
        costFinal1 = cost;
        costFinalInd1 = costInd;
    end
    if j == 2
        costFinal2 = cost;
        costFinalInd2 = costInd;
    end
    if j == 3
        costFinal3 = cost;
        costFinalInd3 = costInd;
    end
    if j == 4
        costFinal4 = cost;
        costFinalInd4 = costInd;
    end
    if j == 5
        costFinal5 = cost;
        costFinalInd5 = costInd;
    end
    if j == 6
        costFinal6 = cost;
        costFinalInd6 = costInd;
    end
end

figure(fig)
hold on
plot(costFinalInd1,costFinal1,'r')
plot(costFinalInd2,costFinal2,'g')
plot(costFinalInd3,costFinal3,'b')
plot(costFinalInd4,costFinal4,'c')
plot(costFinalInd5,costFinal5,'m')
plot(costFinalInd6,costFinal6,'k')
hold off
xlabel('Iterations')
ylabel('Norm(rmiss)/Norm(r0)')
legend('TR = 0.1','TR = 0.05','TR = 0.01','TR = 0.005','TR = 0.001','TR = 0.0005')
title('Cost as Function of Time')
fig = fig + 1;