% 16-714 Advanced Control for Robotics
% Script for Lecture 2: Vehicle Models
%
% The script contains the following examples
% 1. Comparison between CT and DT
% 2. Comparison between ZOH and Euler and RK4
% 3. Simulate and visualize vehicle trajectories in 2D

%% 1. Comparison between CT and DT (add RK4)
name = "single_integrator";
x0 = 0; Tmax = 5;

% controller (state feedback)
u = @(x,t) 5 - x;

% Simulate Trajectory in Continuous Time
simCT = struct('type','CT','T',Tmax,'integrator','ode45');
[tlist_ct, xlist_ct, ulist_ct] = roll_out(name, u, x0, simCT);

% Simulate Trajectory in Discrete Time
dt = 0.1; 
K  = round(Tmax/dt);
simDT_ZOH    = struct('type','DT','K',K,'dt',dt,'integrator','ZOH');
simDT_Euler  = struct('type','DT','K',K,'dt',dt,'integrator','Euler');
simDT_RK4    = struct('type','DT','K',K,'dt',dt,'integrator','RK4');

[tlist_dt,   xlist_dt,   ulist_dt ]  = roll_out(name, u, x0, simDT_ZOH);
[tlist_dte,  xlist_dte,  ulist_dte]  = roll_out(name, u, x0, simDT_Euler);
[tlist_rk4,  xlist_rk4,  ulist_rk4]  = roll_out(name, u, x0, simDT_RK4);

% Plot
figure(1); clf;
subplot(2,1,1); hold on; box on;
plot(tlist_ct,  xlist_ct,  'k', 'LineWidth',1.25);
plot(tlist_dt,  xlist_dt,  'r');
plot(tlist_dte, xlist_dte, 'b');
plot(tlist_rk4, xlist_rk4, 'g');
legend('Continuous Time', ...
       ['ZOH dt = ',num2str(dt)], ...
       ['Euler dt = ',num2str(dt)], ...
       ['RK4 dt = ',num2str(dt)], 'Location','best');
title("State Trajectory");

subplot(2,1,2); hold on; box on;
plot(tlist_ct(1:end-1),  ulist_ct,  'k', 'LineWidth',1.25);
plot(tlist_dt(1:end-1),  ulist_dt,  'r');
plot(tlist_dte(1:end-1), ulist_dte, 'b');
plot(tlist_rk4(1:end-1), ulist_rk4, 'g');
legend('Continuous Time', ...
       ['ZOH dt = ',num2str(dt)], ...
       ['Euler dt = ',num2str(dt)], ...
       ['RK4 dt = ',num2str(dt)], 'Location','best');
title("Control Signal");

%% 2. Comparison between ZOH, Euler, and RK4
name = "double_integrator";
x0 = [0;0]; Tmax = 5;

% controller
u = @(x,t) 5 - x(1) - 2*x(2); % sin(t); % 5 - x(1) - 2*x(2);

% Simulate Trajectory in Continuous Time
simCT = struct('type','CT','T',Tmax,'integrator','ode45');
[tlist_ct, xlist_ct, ulist_ct] = roll_out(name, u, x0, simCT);

% Simulate Trajectory in Discrete Time
dt = 1; 
K  = round(Tmax/dt);
simDT_ZOH    = struct('type','DT','K',K,'dt',dt,'integrator','ZOH');
simDT_Euler  = struct('type','DT','K',K,'dt',dt,'integrator','Euler');
simDT_RK4    = struct('type','DT','K',K,'dt',dt,'integrator','RK4');

[tlist_dt,   xlist_dt,   ulist_dt ]  = roll_out(name, u, x0, simDT_ZOH);
[tlist_dte,  xlist_dte,  ulist_dte]  = roll_out(name, u, x0, simDT_Euler);
[tlist_rk4,  xlist_rk4,  ulist_rk4]  = roll_out(name, u, x0, simDT_RK4);

% Plot
figure(2); clf;

% Position trajectory
subplot(3,1,1); hold on; box on;
plot(tlist_ct,  xlist_ct(1,:),  'k','LineWidth',1.25);
plot(tlist_dt,  xlist_dt(1,:),  'r');
plot(tlist_dte, xlist_dte(1,:), 'b');
plot(tlist_rk4, xlist_rk4(1,:), 'g');
legend('CT','ZOH','Euler','RK4','Location','best');
title("Position Trajectory");
ylabel("x_1");

% Velocity trajectory
subplot(3,1,2); hold on; box on;
plot(tlist_ct,  xlist_ct(2,:),  'k','LineWidth',1.25);
plot(tlist_dt,  xlist_dt(2,:),  'r');
plot(tlist_dte, xlist_dte(2,:), 'b');
plot(tlist_rk4, xlist_rk4(2,:), 'g');
legend('Continuous Time', ...
       ['ZOH dt = ',num2str(dt)], ...
       ['Euler dt = ',num2str(dt)], ...
       ['RK4 dt = ',num2str(dt)], 'Location','best');
title("Velocity Trajectory");
ylabel("x_2");

% Control input
subplot(3,1,3); hold on; box on;
plot(tlist_ct(1:end-1),  ulist_ct,  'k','LineWidth',1.25);
plot(tlist_dt(1:end-1),  ulist_dt,  'r');
plot(tlist_dte(1:end-1), ulist_dte, 'b');
plot(tlist_rk4(1:end-1), ulist_rk4, 'g');
legend('Continuous Time', ...
       ['ZOH dt = ',num2str(dt)], ...
       ['Euler dt = ',num2str(dt)], ...
       ['RK4 dt = ',num2str(dt)], 'Location','best');
title("Control Signal");
xlabel("time (s)");
ylabel("u");

%% 3a. Simulate and visualize unicycle3 trajectory
name = "unicycle3";
x0 = [0;0;0]; Tmax = 5;
x = x0;

% For plotting
figure(3); clf; hold on; axis equal; axis([-10 10 -10 10]); box on; title("Unicycle3");
vehicle.handle = plot([x(1) x(1)+cos(x(3))], [x(2) x(2)+sin(x(3))], ...
    'linewidth',3,'color','r','markersize',50);
set(vehicle.handle,'XDataSource','[x(1) x(1)+cos(x(3))]');
set(vehicle.handle,'YDataSource','[x(2) x(2)+sin(x(3))]');

% Controller
goal = [5;5]; kp = 1; ktheta = 1;
u = @(x, t) [ kp*dot(goal - x(1:2), [cos(x(3)), sin(x(3))]); ...
              ktheta*(atan2(goal(2)-x(2), goal(1)-x(1)) - x(3)) ];

% Simulate Trajectory (DT Euler)
dt = 0.1; 
K  = round(Tmax/dt);
simDT = struct('type','DT','K',K,'dt',dt,'integrator','Euler');
[tlist, xlist, ulist] = roll_out(name, u, x0, simDT);

% Visualize
for k = 2:length(tlist)
    x = xlist(:,k);
    pause(tlist(k)-tlist(k-1));
    refreshdata([vehicle.handle],'caller');
end
plot(xlist(1,:), xlist(2,:), 'k');

%% 3b. Simulate and visualize dynamic bicycle trajectory
name = "bicycle_dynamic";

% State: [X; Y; vx; vy; r; psi]
x0   = [0; 0; 0; 0; 0; 0];
Tmax = 5;

% Goal and controller gains
goal = [5; 5];
p    = par('bicycle_dynamic');

% Control u(x,t) = [Fx; delta_f]
u = @(x,t) bicycle_ctrl(x, goal, p);

% Simulate
dt = 0.01;
K  = round(Tmax/dt);
simDT = struct('type','DT','K',K,'dt',dt,'integrator','RK4');
simDT.stop = @(x,t,log) (norm(x(1:2) - goal) <= 0.1);
[tlist, xlist, ulist] = roll_out(name, u, x0, simDT);

% --- visualize (animated single-track with steering) ---
figure(3); clf; hold on; axis equal; axis([-2 12 -2 10]); box on; title('Bicycle (dynamic)');
x = xlist(:,1);  delta = ulist(2,1);
[rearC, frontC, rearRect, frontRect] = bicycle_geom(x, delta, p);
centerLine = plot([rearC(1) frontC(1)], [rearC(2) frontC(2)], 'b-', 'LineWidth',2);
rearPoly  = plot(rearRect(1,:),  rearRect(2,:),  'r-', 'LineWidth',2);
frontPoly = plot(frontRect(1,:), frontRect(2,:), 'g-', 'LineWidth',2);
plot(goal(1),goal(2),'ko','MarkerFaceColor','y');  % goal

for k = 2:length(tlist)
    x     = xlist(:,k);
    delta = ulist(2,k-1);
    [rearC, frontC, rearRect, frontRect] = bicycle_geom(x, delta, p);
    set(centerLine,'XData',[rearC(1) frontC(1)], 'YData',[rearC(2) frontC(2)]);
    set(rearPoly,  'XData',rearRect(1,:),  'YData',rearRect(2,:));
    set(frontPoly, 'XData',frontRect(1,:), 'YData',frontRect(2,:));
    pause(tlist(k)-tlist(k-1));
    drawnow;
end
plot(xlist(1,:), xlist(2,:), 'k');  % path trace

% ================= helper functions =================
function u = bicycle_ctrl(x, goal, p)
% x = [X;Y;vx;vy;r;psi]
X=x(1); Y=x(2); vx=x(3); vy=x(4); r=x(5); psi=x(6);

% --- steering (Stanley-style) ---
% heading error to goal
psi_des = atan2(goal(2)-Y, goal(1)-X);
e_psi   = wrapToPi_local(psi_des - psi);
% cross-track error in body frame
Rwb = [cos(psi) sin(psi); -sin(psi) cos(psi)]; % world->body
e_vec_b = Rwb * ([goal(1);goal(2)] - [X;Y]);
e_y     = e_vec_b(2);
% delta = heading error + arctan(k*e_y / v)
v_long  = max(abs(vx), p.eps_vx) * sign(vx + (vx==0));
delta   = e_psi + atan2(p.k_cte * e_y, v_long);
delta   = max(-p.delta_max, min(p.delta_max, delta));

% --- longitudinal force to track a target speed ---
dist    = norm([goal(1)-X, goal(2)-Y]);
v_ref   = min(p.v_ref_max, p.k_speed * dist);
Fx_cmd  = p.m * ( p.k_v*(v_ref - vx) - r*vy );     % simple feedforward/feedback
Fx      = max(-p.Fmax, min(p.Fmax, Fx_cmd));

u = [Fx; delta];
end

function [rearC, frontC, rearRect, frontRect] = bicycle_geom(x, delta, p)
% geometry for visualization
X=x(1); Y=x(2); psi=x(6);
R = [cos(psi) -sin(psi); sin(psi) cos(psi)];
Cg = [X; Y];
frontC = Cg + R*[ p.lf; 0];
rearC  = Cg + R*[-p.lr; 0];

Lw = getfield_default(p,'viz_wheel_len',0.8);
Ww = getfield_default(p,'viz_wheel_w',0.25);
rect = [-Lw/2 -Ww/2; Lw/2 -Ww/2; Lw/2 Ww/2; -Lw/2 Ww/2; -Lw/2 -Ww/2]';
Rr = [cos(psi)        -sin(psi);        sin(psi)        cos(psi)];
Rf = [cos(psi+delta)  -sin(psi+delta);  sin(psi+delta)  cos(psi+delta)];
rearRect  = Rr*rect  + rearC;
frontRect = Rf*rect  + frontC;

rearC  = rearC(:).';
frontC = frontC(:).';
end

function val = getfield_default(S, name, default)
if isfield(S,name), val = S.(name); else, val = default; end
end

function ang = wrapToPi_local(ang)
% robust wrap in case Mapping Toolbox isn't available
ang = mod(ang + pi, 2*pi) - pi;
end
