% 16-714 Advanced Control for Robotics
% Script for Lecture 19: Model Reference Adaptive Control

%% Problem Setup
n = 2; m = 2;
% x' = Ax + Bu
sys.A = rand(n);
sys.B = rand(m);
sys.f = @(x,u) sys.A*x + sys.B*u;
% target system x' = Astar * x
sys.Astar = 0*eye(2);
sys.Ahat = eye(n);
sys.Bhat = eye(m);
sys.x0 = [1;-1];
sys.F = 10*eye(n+m) + ones(n+m);
sys.N = 29;
sys.dt = 1;
%% Controller
c = synthesis('MRACx', sys);
[~,xlist] = roll_out(sys, @(x,t) c.u(x,t), sys.x0, 'DT', sys.N, sys.dt, 'Direct');

c.analysis(sys.A, sys.B, sys.N);
%% Ploat

figure(1);clf;hold on;
plot(xlist(:,:)','k')
plot(c.data.V,'r')
%plot(c.data.error','--b')
plot(c.data.error_pos','-*b')
box on;
legend(["x1", "x2", "value function", "a posteriori error x1", "a posteriori error x2"])