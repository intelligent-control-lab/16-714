% 16-714 Advanced Control for Robotics
% Script for Lecture 12 ILC Frequency Domain

%% Problem Specification
sys.dt = 1; % Sampling time
sys.name = 'double_integrator';
[sys.A,sys.B] = getAB('double_integrator',2,'DT',sys.dt,'Euler');
sys.a = [1 -2 1];
sys.b = [sys.dt];
sys.C = [1 0]; % only measuring the position
sys.x0 = [0;0];

% Reference signal
ref = @(t) sin(0.2*t);
simhorizon = 100;
niter = 10;
sigma = 0; % control disturbance; 0 if no disturbance

%% Specify feedback control
kp = 0.1; kd = 0.5;
ufb = @(x,t) kp*(ref(t) - sys.C*x) + kd*((ref(t+sys.dt)-ref(t))/sys.dt - x(2)) + sigma * randn; % proportional error feedback
sys.c = [-kd/sys.dt, kd/sys.dt - kp]; % Transfer function of the controller; Note this controller is not causal

% First trajectory (Iteration 0)
[tlist, x_fb, u_fb] = roll_out(sys.name, ufb, sys.x0, 'DT', simhorizon, sys.dt, 'Euler');

% Plot
figure(1);clf;subplot(211); hold on
plot(0:simhorizon,ref(tlist),'--')
plot(0:simhorizon,x_fb(1,:),'k')
subplot(212)
plot(0:simhorizon-1,u_fb)
%% Specify the learning gain
cILC = synthesis('ILCw', sys);

%% ILC
error = zeros(niter,simhorizon+1);
ff = zeros(niter, simhorizon+10);
e = zeros(1,niter);
string = ["Reference", "Iter 0"]
for i = 1:niter
    error(i,:) = ref(tlist) - x_fb(1,:);
    % filter the error signal
    if i > 1
        ff(i,1:simhorizon+1) = cILC.u(error(i,:), ff(i-1,1:simhorizon+1));
    else
        ff(i,1:simhorizon+1) = cILC.u(error(i,:), zeros(1,simhorizon+1));
    end
    % Set the feedforward controller
    uff = @(x,t) ff(i,t/sys.dt+1) + ufb(x,t);
    % Simulate
    [tlist, x_fb, u_fb] = roll_out(sys.name, uff, sys.x0, 'DT', simhorizon, sys.dt, 'Euler');
    % plot
    subplot(211); hold on
    plot(0:simhorizon, x_fb(1,:), 'color',[i/niter,0,0])
    xlim([0,simhorizon])
    subplot(212); hold on
    plot(0:simhorizon-1, u_fb,'color',[i/niter,0,0])
    % Compute error norm
    e(i) = norm(error(i,:));
    string = [string,"Iter "+num2str(i)];
end
subplot(211);
legend(string);title('Output')
subplot(212);
title('Input')
figure(2);clf;
plot(0:niter-1, e);
title('Norm Error in Different ILC Iterations')
