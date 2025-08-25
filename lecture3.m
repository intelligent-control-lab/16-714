% 16-714 Advanced Control for Robotics
% Script for Lecture 3: Manipulator Models
%
% The script contains the following examples
% 1. Simulate and visualize manipulator trajectories in 3D
% 2. Compare between joint space control and Cartesian space control

%% 1. Simulate and visualize manipulator trajectories in 3D
global robot;
robot = loadrobot('kinovaGen3','DataFormat','row','Gravity',[0 0 -9.81]);
x0 = homeConfiguration(robot);
name = "single_integrator";
Tmax = 5;

figure(4); clf; axis([-1 1 -1 1 -0.1 1.5]);
show(robot, x0, 'PreservePlot', false, 'Frames', 'off');

u = @(x, t) randn(7,1);

% Simulate Trajectory (joint-space single-integrator with Euler)
dt = 0.1; 
K  = round(Tmax/dt);
simDT = struct('type','DT','K',K,'dt',dt,'integrator','Euler');
[tlist, xlist, ulist] = roll_out(name, u, x0', simDT);

% Visualize
for k = 2:length(tlist)
    x = xlist(:,k);
    pause(tlist(k)-tlist(k-1));
    show(robot, x', 'PreservePlot', false, 'Frames', 'off');
end
axis([-1 1 -1 1 -0.1 1.5]);

%% 2. Compare between joint space control and Cartesian space control
global robot;
robot = loadrobot('kinovaGen3','DataFormat','row','Gravity',[0 0 -9.81]);

goal = [0.4,0,0.6];
x0 = homeConfiguration(robot);
endEffector = "EndEffector_Link";
taskInit = getTransform(robot, x0, endEffector);
taskFinal = trvec2tform(goal) * axang2tform([0 1 0 pi]);
ik = inverseKinematics("RigidBodyTree", robot);
xT = ik(endEffector, taskFinal, [1 1 1 1 1 1], x0); xT = mod(xT, 2*pi);

% 2.1. Joint space control
name = "single_integrator";
Tmax = 5;
figure(5); clf; axis([-1 1 -1 1 -0.1 1.5]);
show(robot, x0, 'PreservePlot', false, 'Frames', 'off');

u = @(x, t) 0.7*(xT' - x);

dt = 0.1; 
K  = round(Tmax/dt);
simDT = struct('type','DT','K',K,'dt',dt,'integrator','Euler');
[tlist, xlist, ulist] = roll_out(name, u, x0', simDT);

clist = zeros(length(tform2trvec(taskInit)), length(tlist));
clist(:,1) = tform2trvec(taskInit)';

% Visualize
for k = 2:length(tlist)
    x = xlist(:,k);
    pause(tlist(k)-tlist(k-1));
    show(robot, x', 'PreservePlot', false, 'Frames', 'off');
    clist(:,k) = tform2trvec(getTransform(robot, x', endEffector))';
end
hold on;
plot3(clist(1,:), clist(2,:), clist(3,:), 'k', 'LineWidth', 2);
axis([-1 1 -1 1 -0.1 1.5]);

% 2.2. Cartesian space control
name = "arm_first_order";
Tmax = 5;
figure(6); clf;
show(robot, x0, 'PreservePlot', false, 'Frames', 'off');

c = tform2trvec(taskInit);
u = @(x, t) 0.5 .* pinv(geometricJacobian(robot, ...
        ik("EndEffector_Link", trvec2tform(x'), [1 1 1 1 1 1], homeConfiguration(robot)), ...
        "EndEffector_Link")) * [goal' - x; zeros(3,1)];

dt = 0.1; 
K  = round(Tmax/dt);
simDT = struct('type','DT','K',K,'dt',dt,'integrator','Euler');
[tlist, xlist, ulist] = roll_out(name, u, c', simDT);

% Visualize
x = x0;
for k = 2:length(tlist)
    pause(tlist(k)-tlist(k-1));
    x = ik(endEffector, trvec2tform(xlist(:,k)'), [1 1 1 1 1 1], x);
    show(robot, x, 'PreservePlot', false, 'Frames', 'off');
end
hold on;
plot3(xlist(1,:), xlist(2,:), xlist(3,:), 'k', 'LineWidth', 2);
axis([-1 1 -1 1 -0.1 1.5]);
