% 16-714 Advanced Control for Robotics
% Script for Lecture 19: Model Reference Adaptive Control

%% Problem Setup
n = 2; m = 2;
% x' = A x + B u
sys.A    = rand(n);
sys.B    = rand(n,m);
sys.f_dt    = @(x,u) sys.A*x + sys.B*u;

% target system x' = A* x (reference dynamics)
sys.Astar = zeros(n);
sys.Ahat  = eye(n);
sys.Bhat  = eye(n,m);

sys.x0 = [1;-1];
sys.F  = 10*eye(n+m) + ones(n+m);
sys.N  = 29;
sys.dt = 1;

% New roll_out sim/options (DT)
sim = struct('type','DT', 'K',sys.N, 'dt',sys.dt, 'integrator','Direct');
opts = struct();  % now the observation and parameter adaptation all happen inside the controller

%% Controller
c = synthesis('MRACx', sys);

% Roll_out signature:
% [tlist, xlist, ulist, log] = roll_out(model, ctrl, x0, sim, opts)
[~, xlist] = roll_out(sys, @(x,t) c.u(x,t), sys.x0, sim, opts);

% Post-analysis
c.analysis(sys.A, sys.B, sys.N);

%% Plot
figure(1); clf; hold on;
plot(xlist(:,:).','k')
plot(c.data.V,'r')
plot(c.data.error_pos','-*b')
box on;
legend(["x1","x2","value function","a posteriori error x1","a posteriori error x2"])
title('MRAC: States, Lyapunov Function, and A-Posteriori Error');
xlabel('k'); grid on;