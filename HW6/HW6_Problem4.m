% This file runs solution for HW6 Problem 4
%
% Alec Lewandowski
% 12/1/25
clear variables
close all
clc
format long
fig = 1;

mu = 0.01215;
C = 3.18;

%% Problem 4.i
% Defining ICs
% xcross = 5000;
% x0 = 1.0016;
% y0 = 0; %x-axis crossings
% u0 = 0.53;
% 
% r1IC = x0 + mu; % bc y0 is 0
% r2IC = x0 - 1 + mu;
% UIC = 0.5*x0^2 + (1-mu)/r1IC + mu/r2IC;
% v0 = sqrt(2*UIC - C - u0^2);
% 
% options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12,'Events', @(t,state)xAxisCross(t, state));
% % From part (c) search
% x0 = 0.849288;
% u0 = 0.009483;
% v0 = 0.100315;
% state0 = [x0; y0; 0; u0; v0; 0];
% pcx = [];
% pcu = [];
% pcv = [];
% pct = [];
% xtraj = x0;
% ytraj = y0;
% utraj = u0;
% vtraj = v0;
% ttraj = 0;
% 
% state = state0;
% crosscount = 0;
% tNow = 0;
% 
% while crosscount < xcross
%     [tvec, traj, te, ye, ie] = ode45(@(t,state)EarthMoonCR3BP(t,state,mu),[tNow tNow+100],state0,options);
% 
%     xtraj = [xtraj; traj(2:end,1)];
%     ytraj = [ytraj; traj(2:end,2)];
%     utraj = [utraj; traj(2:end,4)];
%     vtraj = [vtraj; traj(2:end,5)];
%     ttraj = [ttraj; tvec(2:end)];
% 
%     if ~isempty(te)
%         for i = 1:length(te)
%             if (ye(i,4)) > 0
%                 crosscount = crosscount + 1;
%                 pcx = [pcx; ye(i,1)];
%                 pcu = [pcu; ye(i,4)];
%                 pcv = [pcv; ye(i,5)];
%                 pct = [pct; te(i)];
% 
%                 if crosscount >= xcross
%                     break
%                 end
%             end
%         end
%     end
% 
%     state = traj(end,:)';
%     tNow = tvec(end);
% end
% 
% [t_out, states_out, poincare] = pcProp(x0, u0, v0, mu, xcross);
% 
% xtraj = states_out(:,1);
% utraj = states_out(:,4);
% vtraj = states_out(:,5);
% 
% pcx = poincare(:,1);
% pcu = poincare(:,3);
% pcv = poincare(:,4);
% 
% figure(fig)
% plot3(xtraj,utraj,vtraj)
% xlabel('X Pos');
% ylabel('U Velo');
% zlabel('V Velo');
% title('Trajectory in CR3BP');
% grid on;
% fig = fig + 1;
% 
% figure(fig)
% plot(xtraj,utraj)
% xlabel('X Pos [LU]')
% ylabel('U Velo [LU/TU]')
% title('{x,u} View of Trajectory')
% grid on
% fig = fig + 1;
% 
% figure(fig)
% plot(xtraj,vtraj)
% xlabel('X Pos [LU]')
% ylabel('V Velo [LU/TU]')
% title('{x,u} View of Trajectory')
% grid on
% fig = fig + 1;
% 
% figure(fig);
% scatter3(pcx, pcu, pcv, 10, 1:length(pcx), 'filled');
% colormap(jet);
% colorbar;
% grid on;
% xlabel('X Pos [LU]');
% ylabel('U Velo [LU/TU]');
% zlabel('V Velo [LU/TU]');
% title('3D Poincaré Map: Ascending X-Crossing');
% fig = fig + 1;
% 
% figure(fig);
% scatter(pcx,pcu)
% xlabel('X Pos [LU]')
% ylabel('U Velo [LU/TU]')
% title('{x,u} View of Poincare')
% grid on
% fig = fig + 1;
% 
% figure(fig);
% scatter(pcx,pcv)
% xlabel('X Pos [LU]')
% ylabel('V Velo [LU/TU]')
% title('{x,v} View of Poincare')
% grid on
% fig = fig + 1;
% 
% %% Part (c): Find periodic orbits from part (a)
% fprintf('=== Part (c): Periodic Orbit Search from Part (a) ===\n');
% % Look for crossings with u ≈ 0 (perpendicular crossing)
% idx_perp = find(abs(pcu) < 0.01);
% if ~isempty(idx_perp)
%     fprintf('Approximate perpendicular crossings found at:\n');
%     for i = 1:min(5, length(idx_perp))
%         fprintf('  x = %.6f, u = %.6f, v = %.6f\n', ...
%             pcx(idx_perp(i),1), pcu(idx_perp(i),2), pcv(idx_perp(i),3));
%     end
% 
%     % Refine search around first candidate
%     x0_refine = pcx(idx_perp(1),1);
%     fprintf('\nRefining search around x0 = %.6f\n', x0_refine);
% 
%     % Tighter grid
%     x0_refined = linspace(x0_refine-0.01, x0_refine+0.01, 20);
%     u0_refined = 0.0;
% 
%     poincare_c = [];
%     for x0 = x0_refined
%         r1IC = x0 + mu; % bc y0 is 0
%         r2IC = x0 - 1 + mu;
%         UIC = 0.5*x0^2 + (1-mu)/r1IC + mu/r2IC;
%         v0 = sqrt(2*UIC - C - u0_refined^2);
%         if isnan(v0), continue; end
%         [~, ~, pm] = pcProp(x0, u0_refined, v0, mu, 500);
%         poincare_c = [poincare_c; pm];
%     end
% 
%     % Plot refined search
%     figure('Position', [100 100 1200 400]);
%     subplot(1,3,1);
%     plot(poincare_c(:,1), poincare_c(:,2), 'b.', 'MarkerSize', 3);
%     grid on; xlabel('x [LU]'); ylabel('u [LU/TU]');
%     title('Part (c): Refined Poincaré Map (x,u)');
% 
%     % Find best periodic orbit candidate
%     idx_best = find(abs(poincare_c(:,2)) < 0.001);
%     if ~isempty(idx_best)
%         x0_best = poincare_c(idx_best(1),1);
%         fprintf('Improved estimate: x0 = %.8f, u0 = 0.0\n', x0_best);
% 
%         % Propagate periodic orbit
%         r1IC = x0_best + mu; % bc y0 is 0
%         r2IC = x0_best - 1 + mu;
%         UIC = 0.5*x0_best^2 + (1-mu)/r1IC + mu/r2IC;
%         v0_best = sqrt(2*UIC - C);
% 
% 
%         [t_po, states_po, pm_po] = pcProp(x0_best, 0, v0_best, mu, 10);
% 
%         if size(pm_po,1) >= 2
%             period = pm_po(2,4);
%             fprintf('Estimated period: %.6f TU\n', period);
% 
%             subplot(1,3,2);
%             plot(states_po(:,1), states_po(:,2), 'r-', 'LineWidth', 2);
%             hold on;
%             plot3(-mu, 0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
%             plot3(1-mu, 0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
%             grid on; xlabel('x [LU]'); ylabel('y [LU]');
%             title(sprintf('Periodic Orbit (T=%.3f TU)', period));
%             axis equal;
% 
%             subplot(1,3,3);
%             plot(poincare_c(:,1), poincare_c(:,2), 'b.', 'MarkerSize', 3);
%             hold on;
%             plot(x0_best, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
%             grid on; xlabel('x [LU]'); ylabel('u [LU/TU]');
%             title('Periodic Orbit Location on P.M.');
%         end
%     end
% end
% 
% 
% %% Problem 4.b
% % Grid of initial conditions
% 
% x0 = 0.849288;
% u0 = 0.009483;
% v0 = 0.100315;


x0_grid = 1.05:0.005:1.1;     % x0 range
u0_grid = 0.0:0.07:0.7;       % u0 range

pcxtot = [];
pcutot = [];
pcvtot = [];
ictot = [];      % Track which IC each point came from

% ODE options with event detection
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'Events', @(t,state)xAxisCross(t, state));

xcross = 1000;   % Target per trajectory
count = 0;
totICs = length(x0_grid) * length(u0_grid);

% Progress tracking
tic;
valid_trajectories = 0;
skipped_trajectories = 0;

for i = 1:length(x0_grid)
    for j = 1:length(u0_grid)
        count = count + 1;
        x0 = x0_grid(i);
        u0 = u0_grid(j);
        y0 = 0;  

        r1 = sqrt((x0 + mu)^2 + y0^2);
        r2 = sqrt((x0 - 1 + mu)^2 + y0^2);
        U = (x0^2 + y0^2)/2 + (1-mu)/r1 + mu/r2 + mu*(1-mu)/2;
        v0_squared = 2*U - C - u0^2;

        % % Skip if initial conditions are incompatible with Jacobi constant
        % if v0_squared < 0
        %     skipped_trajectories = skipped_trajectories + 1;
        %     if mod(count, 10) == 0
        %         fprintf('[%d/%d] Skipped IC (x0=%.3f, u0=%.3f) - incompatible with C\n', ...
        %                 count, totICs, x0, u0);
        %     end
        %     continue;
        % end

        v0 = sqrt(v0_squared);  % Positive v0
        state0 = [x0; y0; 0; u0; v0; 0];

        [t_out, states_out, poincare] = pcProp(x0, u0, v0, mu, xcross);

        xtraj = states_out(:,1);
        utraj = states_out(:,4);
        vtraj = states_out(:,5);

        pcx = poincare(:,1);
        pcu = poincare(:,3);
        pcv = poincare(:,4);

        % Integration parameters
        state = state0;
        tNow = 0;
        crosscount = 0;
        tmax = 500;  % Maximum integration tim

        pcxtot = [pcxtot; pcx];
        pcutot = [pcutot; pcu];
        pcvtot = [pcvtot; pcv];
        ictot = [ictot; count * ones(length(pcx), 1)];

    end
end

% Plot 1: 3D Poincare Map {x, u, v}
figure(fig);
scatter3(pcxtot, pcutot, pcvtot, 3, ictot, 'filled');
colormap(jet);
cb = colorbar;
ylabel(cb, 'Initial Condition Index');
grid on;
xlabel('X Pos [LU]');
ylabel('U Velo [LU/TU]');
zlabel('V Velo [LU/TU]');
title('3D Poincaré Map');
view(45, 25);
rotate3d on;
fig = fig + 1;


figure(fig);
plot(pcxtot, pcutot, 'b.','MarkerSize',10);
grid on;
xlabel('X Pos [LU]');
ylabel('U Velo [LU/TU]');
title('Poincaré Map {x, u} - Full View');
fig = fig + 1;


figure(fig);
plot(pcxtot, pcvtot, 'r.','MarkerSize',10);
grid on;
xlabel('X Pos [LU]');
ylabel('V Velo [LU/TU]');
title('Poincaré Map {x, v} - Full View');
fig = fig + 1;








