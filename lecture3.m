% 16-714 Advanced Control for Robotics
% Script for Lecture 3: Matlab Tutorial

% The script contains the following examples
% 1. comparison between CT and DT
% 2. Comparison between ZOH and Euler
% 3. Simulate and visualize vehicle trajectories in 2D
% 4. Simulate and visualize manipulator trajectories in 3D
% 5. Compare between joint space control and Cartesian space control


%% 1. Comparison between CT and DT
name = 'single_integrator';
x0 = 0; Tmax = 5;

%controller
u = @(x,t) 5-x;

% Simulate Trajectory in Continuous Time
[tlist_ct, xlist_ct, ulist_ct] = roll_out(name, u, x0, 'CT', Tmax);

% Simulate Trajectory in Discrete Time
dt = 0.1; 
[tlist_dt, xlist_dt, ulist_dt] = roll_out(name, u, x0, 'DT', Tmax/dt, dt, 'ZOH');
[tlist_dte, xlist_dte, ulist_dte] = roll_out(name, u, x0, 'DT', Tmax/dt, dt, 'Euler');

% Plot
figure(1);clf;hold on;
subplot(2,1,1); hold on;
plot(tlist_ct,xlist_ct,'k');
plot(tlist_dt,xlist_dt,'r');
plot(tlist_dte,xlist_dte,'b');
legend('Continuous Time',['Discrete Time dt = ',num2str(dt)]);
title("State Trajectory")
subplot(2,1,2); hold on;
plot(tlist_ct(1:end-1),ulist_ct,'k');
plot(tlist_dt(1:end-1),ulist_dt,'r');
plot(tlist_dte(1:end-1),ulist_dte,'b');
legend('Continuous Time',['Discrete Time dt = ',num2str(dt)]);
title("Control Signal")

%% 2. Comparison between ZOH and Euler
name = 'double_integrator';
x0 = [0;0]; Tmax = 5;

%controller
u = @(x,t) 1;%sin(t);%5-x(1)-2*x(2);

% Simulate Trajectory in Continuous Time
[tlist_ct, xlist_ct, ulist_ct] = roll_out(name, u, x0, 'CT', Tmax);

% Simulate Trajectory in Discrete Time
dt = 0.5; 
[tlist_dt, xlist_dt, ulist_dt] = roll_out(name, u, x0, 'DT', Tmax/dt, dt, 'ZOH');
[tlist_dte, xlist_dte, ulist_dte] = roll_out(name, u, x0, 'DT', Tmax/dt, dt, 'Euler');

% Plot
figure(1);clf;hold on;
subplot(2,1,1); hold on;
plot(tlist_ct,xlist_ct,'k');
plot(tlist_dt,xlist_dt,'r');
plot(tlist_dte,xlist_dte,'b');
legend('Continuous Time (position)','Continuous Time (velocity)', ...
    ['ZOH position dt = ',num2str(dt)], ['ZOH velocity dt = ',num2str(dt)], ...
    ['Euler position dt = ',num2str(dt)], ['Euler velocity dt = ',num2str(dt)]);
title("State Trajectory")
subplot(2,1,2); hold on;
plot(tlist_ct(1:end-1),ulist_ct,'k');
plot(tlist_dt(1:end-1),ulist_dt,'r');
plot(tlist_dte(1:end-1),ulist_dte,'b');
legend('Continuous Time',['ZOH dt = ',num2str(dt)], ...
    ['Euler dt = ',num2str(dt)]);
title("Control Signal")

%% 3. Simulate and visualize vehicle trajectories in 2D
name = 'unicycle3';
x0 = [0;0;0]; Tmax = 5;
x = x0;
% For plotting
figure(1);clf;hold on;axis([-10 10 -10 10]);box on;title("Unicycle3")
vehicle.handle = plot([x(1) x(1)+cos(x(3))],[x(2) x(2)+sin(x(3))], ...
    'linewidth',3,'color','r','markersize',50);
set(vehicle.handle,'XDataSource','[x(1) x(1)+cos(x(3))]');
set(vehicle.handle,'YDataSource','[x(2) x(2)+sin(x(3))]');

% Set up controller
goal = [5;5]; kp = 1;ktheta = 1;
u = @(x, t) [kp*dot(goal - x(1:2),[cos(x(3)),sin(x(3))]);...
    ktheta*(atan2(goal(2)-x(2),goal(1)-x(1)) - x(3))]; 

% Simulate Trajectory
dt = 0.1;
[tlist, xlist, ulist] = roll_out(name, u, x0, 'DT', Tmax/dt, dt, 'Euler');

% visualize
for k = 2:length(tlist)
    x = xlist(:,k);
    pause(tlist(k)-tlist(k-1));
    refreshdata([vehicle.handle],'caller');
end
plot(xlist(1,:), xlist(2,:),'k');


%% 4. Simulate and visualize manipulator trajectories in 3D
global robot;
robot = loadrobot('kinovaGen3','DataFormat','row','Gravity',[0 0 -9.81]);
x0 = homeConfiguration(robot);
name = 'single_integrator';
Tmax = 5;

figure(1); clf; axis([-1 1 -1 1 -0.1 1.5]);
show(robot,x0,'PreservePlot',false,'Frames','off');

u = @(x, t) randn(7,1);

% Simulate Trajectory
dt = 0.1; 
[tlist, xlist, ulist] = roll_out(name, u, x0', 'DT', Tmax/dt, dt, 'Euler');

% visualize
for k = 2:length(tlist)
    x = xlist(:,k);
    pause(tlist(k)-tlist(k-1));
    show(robot,x','PreservePlot',false,'Frames','off');
end
axis([-1 1 -1 1 -0.1 1.5]);

%% 5. Compare between joint space control and Cartesian space control
global robot;
robot = loadrobot('kinovaGen3','DataFormat','row','Gravity',[0 0 -9.81]);

goal = [0.4,0,0.6];
x0 = homeConfiguration(robot);
endEffector = "EndEffector_Link";
taskInit = getTransform(robot,x0,endEffector);
taskFinal = trvec2tform(goal)*axang2tform([0 1 0 pi]);
ik = inverseKinematics("RigidBodyTree",robot);
xT = ik(endEffector,taskFinal,[1 1 1 1 1 1],x0);xT = mod(xT,2*pi);

% 5.1. Joint space control
name = 'single_integrator';
Tmax = 5;
figure(1); clf; axis([-1 1 -1 1 -0.1 1.5]);
show(robot,x0,'PreservePlot',false,'Frames','off');

u = @(x, t) 0.7*(xT'-x);

% Simulate Trajectory
dt = 0.1; 
[tlist, xlist, ulist] = roll_out(name, u, x0', 'DT', Tmax/dt, dt, 'Euler');
clist = zeros(length(tform2trvec(taskInit)),length(tlist)); clist(:,1) = tform2trvec(taskInit)';
% visualize
for k = 2:length(tlist)
    x = xlist(:,k);
    pause(tlist(k)-tlist(k-1));
    show(robot,x','PreservePlot',false,'Frames','off');
    clist(:,k) = tform2trvec(getTransform(robot,x',endEffector))';
end
hold on;
plot3(clist(1,:),clist(2,:),clist(3,:),'k','LineWidth',2);
axis([-1 1 -1 1 -0.1 1.5]);

% 5.2. Cartesian space control
name = 'arm_first_order';
Tmax = 5;
figure(2); clf; 
show(robot,x0,'PreservePlot',false,'Frames','off');
c = tform2trvec(taskInit);
u = @(x, t) 0.5.*pinv(geometricJacobian(robot, ik("EndEffector_Link",trvec2tform(x'), ...
    [1 1 1 1 1 1],homeConfiguration(robot)), "EndEffector_Link"))*[goal' - x;zeros(3,1)];

% Simulate Trajectory
dt = 0.1; 
[tlist, xlist, ulist] = roll_out(name, u, c', 'DT', Tmax/dt, dt, 'Euler');

% visualize
x = x0;
for k = 2:length(tlist)
    pause(tlist(k)-tlist(k-1));
    x = ik(endEffector,trvec2tform(xlist(:,k)'),[1 1 1 1 1 1],x);
    show(robot,x,'PreservePlot',false,'Frames','off');
end
hold on;
plot3(xlist(1,:),xlist(2,:),xlist(3,:),'k','LineWidth',2);
axis([-1 1 -1 1 -0.1 1.5]);