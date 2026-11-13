% This function runs solution for HW5 Problem 3
%
% Alec Lewandowski
% 11/12/25
clear variables
close all
clc
format long

% Defining the ICs for reference orbit
r1 = [1; 0; 0];
v1 = [0; 1; 0];
STM1 = eye(6);
TOF = 14;
GM = 1;
athrust = 0.1;
fig = 1;

% Setting up and integrating the EOMs in cartesian coords
tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',tol);
[tref, yref] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM,athrust),0:1:TOF,[r1; v1; STM1(:)],options);
[trefs, yrefs] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM,athrust),[0 TOF],[r1; v1; STM1(:)],options);
refPhiFinal = reshape(yref(end,7:end),[6,6]);
finalPosRef = yrefs(end,1:3)'
finalVelRef = yrefs(end,4:6)';

% Now propogating the guidance matrix from tfinal to t0
[tGuidance, yGuidance] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM,athrust),TOF:-1:0,[finalPosRef; finalVelRef; STM1(:)],options); 
PhiatT8 = reshape(yGuidance(7,7:end),[6,6]);
BatT8 = PhiatT8(1:3,4:6)

phiGuide = zeros(6,6,length(tGuidance));
B = zeros(3,3,length(tGuidance));
D = zeros(3,3,length(tGuidance));
for i = 1:length(tGuidance)
    phiTemp = reshape(yGuidance(i,7:end),[6,6]);
    B(:,:,i) = phiTemp(1:3, 4:6);
    D(:,:,i) = phiTemp(4:6, 4:6);
    phiGuide(:,:,i) = phiTemp;
end

B = flip(B,3);
D = flip(D,3);

% Flying truth trajectory with Gm = 1.001
GM = 1.001;
TOFvec = 0:1:14;
refPos = zeros(length(TOFvec),3);
refVel = zeros(length(TOFvec),3);
delVTot = zeros(length(TOFvec),3);
maneuverCounter = 0;

rInit = r1;
vInit = v1;
tInit = 0;
STMInit = phiGuide(:,:,1);
posTruthX = [];
posTruthY = [];
truthXPts = [];
truthYPts = [];
relPos = zeros(length(TOFvec),3);
relPos(1,:) = [0; 0; 0];
relPos(end,:) = [0; 0; 0];
magRelPos = zeros(length(tref),3);
magManeuver = zeros(length(tref),1);

for i = 1:length(TOFvec)-1
    tInit = TOFvec(i);
    tFinal = TOFvec(i+1);
    applyManeuver = false;
    [tTruth, yTruth] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM,athrust),[tInit tFinal],[rInit; vInit; STMInit(:)],options);
    phiTempTruth = reshape(yTruth(end,7:end),[6,6]);
    posTemp = yTruth(end,1:3)';
    velTemp = yTruth(end,4:6)';
    delR = posTemp - yref(i+1,1:3)';
    delV = velTemp - yref(i+1,4:6)';
    relPos(i+1,:) = delR;
    magRelPos(i+1,:) = norm(delR);
    normErr = norm(delR);
    if normErr > 0.01
        applyManeuver = true;
    end
    if applyManeuver && i ~= 14
        delVTot(i+1,:) = D(:,:,i+1)/B(:,:,i+1)*delR - delV;
        maneuverCounter = maneuverCounter + 1;
        STMInit = phiTempTruth;
        rInit = posTemp;
        vInit = velTemp + delVTot(i+1,:)';
        magManeuver(i+1) = norm(delVTot(i+1,:));
    else
        STMInit = phiTempTruth;
        rInit = posTemp;
        vInit = velTemp;
    end
    posTruthX = [posTruthX, yTruth(:,1)'];
    posTruthY = [posTruthY, yTruth(:,2)'];
    truthXPts = [truthXPts, yTruth(end,1)];
    truthYPts = [truthYPts, yTruth(end,2)];
end

DAtT6 = D(:,:,7)

delVtotal = sum(sum(abs(delVTot)))
maneuverCounter

figure(fig)
hold on
plot(yrefs(:,1), yrefs(:,2),'r')
plot(posTruthX, posTruthY,'b--')
plot(truthXPts,truthYPts,'k*')
plot(1,0,'bs','MarkerSize',25)
quiver(truthXPts, truthYPts, delVTot(2:end,1)', delVTot(2:end,2)','Color', 'g')
hold off
axis equal
grid on
xlabel('X [LU]')
ylabel('Y [LU]')
title('Feedback Guidance with STM')
legend('Ref Traj','Truth Traj','Maneuver Pts','Start Pt','Delta Vs')
fig = fig + 1;

figure(fig)
plot(relPos(:,1),relPos(:,2),'r')
xlabel('X [LU]')
ylabel('Y [LU]')
title('Relative Trajectory')
fig = fig + 1;

figure(fig)
plot(tref,magRelPos,'r')
xlabel('Time [TU]')
ylabel('deltaR Magnitude [LU]')
title('Magnitude of the Delta Position vs Time')
fig = fig + 1;

figure(fig)
plot(tref,magManeuver,'r')
xlabel('Time [TU]')
ylabel('Maneuver deltaV Magnitude [LU]')
title('Magnitude of the Delta Vs Maneuvers vs Time')
fig = fig + 1;

finalPosTruth = yTruth(end,1:3)'
finalVelTruth = yTruth(end,4:6)'










%% Different strategies and GMs

% Testing different GMs
GM = [0.75 0.85 0.95 1 1.05 1.15 1.25];
TOFvec = 0:1:14;

figure(fig)
hold on
for k = 1:length(GM)
    refPos = zeros(length(TOFvec),3);
    refVel = zeros(length(TOFvec),3);
    delVTot = zeros(length(TOFvec),3);
    maneuverCounter = 0;
    
    rInit = r1;
    vInit = v1;
    tInit = 0;
    STMInit = phiGuide(:,:,1);
    posTruthX = [];
    posTruthY = [];
    truthXPts = [];
    truthYPts = [];
    relPos = zeros(length(TOFvec),3);
    relPos(1,:) = [0; 0; 0];
    relPos(end,:) = [0; 0; 0];
    magRelPos = zeros(length(tref),3);
    magManeuver = zeros(length(tref),1);
    
    
    for i = 1:length(TOFvec)-1
        tInit = TOFvec(i);
        tFinal = TOFvec(i+1);
        applyManeuver = false;
        [tTruth, yTruth] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM(k),athrust),[tInit tFinal],[rInit; vInit; STMInit(:)],options);
        phiTempTruth = reshape(yTruth(end,7:end),[6,6]);
        posTemp = yTruth(end,1:3)';
        velTemp = yTruth(end,4:6)';
        delR = posTemp - yref(i+1,1:3)';
        delV = velTemp - yref(i+1,4:6)';
        relPos(i+1,:) = delR;
        magRelPos(i+1,:) = norm(delR);
        normErr = norm(delR);
        if normErr > 0.01
            applyManeuver = true;
        end
        if applyManeuver && i ~= 14
            delVTot(i+1,:) = D(:,:,i+1)/B(:,:,i+1)*delR - delV;
            maneuverCounter = maneuverCounter + 1;
            STMInit = phiTempTruth;
            rInit = posTemp;
            vInit = velTemp + delVTot(i+1,:)';
            magManeuver(i+1) = norm(delVTot(i+1,:));
        else
            STMInit = phiTempTruth;
            rInit = posTemp;
            vInit = velTemp;
        end
        posTruthX = [posTruthX, yTruth(:,1)'];
        posTruthY = [posTruthY, yTruth(:,2)'];
        truthXPts = [truthXPts, yTruth(end,1)];
        truthYPts = [truthYPts, yTruth(end,2)];
    end
    
    delVtotal(k) = sum(sum(abs(delVTot)));

    plot(tref,magRelPos)
end
hold off
xlabel('Time [TU]')
ylabel('deltaR Magnitude [LU]')
title('Magnitude of the Delta Position vs Time')
legend('GM=0.75','GM=0.85','GM=0.95','GM=1','GM=1.05','GM=1.15','GM=1.25','location','northwest')
fig = fig + 1;



%% Additionally can change interval times
intTimes = [0.001 0.01 0.1 1];
GM = 1.001;

figure(fig)
hold on
for k = 1:length(intTimes)
    % Setting up and integrating the EOMs in cartesian coords
    tol = 1e-12;
    options = odeset('RelTol',tol,'AbsTol',tol);
    [tref, yref] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM,athrust),0:intTimes(k):TOF,[r1; v1; STM1(:)],options);
    [trefs, yrefs] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM,athrust),[0 TOF],[r1; v1; STM1(:)],options);
    refPhiFinal = reshape(yref(end,7:end),[6,6]);
    finalPosRef = yrefs(end,1:3)';
    finalVelRef = yrefs(end,4:6)';

    % Now propogating the guidance matrix from tfinal to t0
    [tGuidance, yGuidance] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM,athrust),TOF:-intTimes(k):0,[finalPosRef; finalVelRef; STM1(:)],options); 

    phiGuide = zeros(6,6,length(tGuidance));
    B = zeros(3,3,length(tGuidance));
    D = zeros(3,3,length(tGuidance));
    for i = 1:length(tGuidance)
        phiTemp = reshape(yGuidance(i,7:end),[6,6]);
        B(:,:,i) = phiTemp(1:3, 4:6);
        D(:,:,i) = phiTemp(4:6, 4:6);
        phiGuide(:,:,i) = phiTemp;
    end

    B = flip(B,3);
    D = flip(D,3);

    % Flying truth trajectory with Gm = 1.001
    GM = 1.001;
    TOFvec = 0:intTimes(k):14;
    refPos = zeros(length(TOFvec),3);
    refVel = zeros(length(TOFvec),3);
    delVTot = zeros(length(TOFvec),3);
    maneuverCounter = 0;

    rInit = r1;
    vInit = v1;
    tInit = 0;
    STMInit = phiGuide(:,:,1);
    posTruthX = [];
    posTruthY = [];
    truthXPts = [];
    truthYPts = [];
    relPos = zeros(length(TOFvec),3);
    relPos(1,:) = [0; 0; 0];
    relPos(end,:) = [0; 0; 0];
    magRelPos = zeros(length(tref),3);
    magManeuver = zeros(length(tref),1);

    for i = 1:length(TOFvec)-1
        tInit = TOFvec(i);
        tFinal = TOFvec(i+1);
        applyManeuver = false;
        [tTruth, yTruth] = ode45(@(t,y)RelTwoBodyThrustProp(t,y,GM,athrust),[tInit tFinal],[rInit; vInit; STMInit(:)],options);
        phiTempTruth = reshape(yTruth(end,7:end),[6,6]);
        posTemp = yTruth(end,1:3)';
        velTemp = yTruth(end,4:6)';
        delR = posTemp - yref(i+1,1:3)';
        delV = velTemp - yref(i+1,4:6)';
        relPos(i+1,:) = delR;
        magRelPos(i+1,:) = norm(delR);
        normErr = norm(delR);
        if normErr > 0.01
            applyManeuver = true;
        end
        if applyManeuver && i ~= 14
            delVTot(i+1,:) = D(:,:,i+1)/B(:,:,i+1)*delR - delV;
            maneuverCounter = maneuverCounter + 1;
            STMInit = phiTempTruth;
            rInit = posTemp;
            vInit = velTemp + delVTot(i+1,:)';
            magManeuver(i+1) = norm(delVTot(i+1,:));
        else
            STMInit = phiTempTruth;
            rInit = posTemp;
            vInit = velTemp;
        end
        posTruthX = [posTruthX, yTruth(:,1)'];
        posTruthY = [posTruthY, yTruth(:,2)'];
        truthXPts = [truthXPts, yTruth(end,1)];
        truthYPts = [truthYPts, yTruth(end,2)];
    end
    plot(tref,magRelPos)
end
hold off
xlabel('Time [TU]')
ylabel('deltaR Magnitude [LU]')
title('Magnitude of the Delta Position vs Time')
legend('intTime = 0.001','intTime = 0.01','intTime = 0.1','intTime = 1','location','northwest')
fig = fig + 1;