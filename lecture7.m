% 16-714 Advanced Control for Robotics
% Script for Lecture 7: LQR 

clear; clc; close all;

%% Problem Parameters
dyn_name = 'double_integrator';       % key for resolve_dynamics
x0   = [10; 10; 1; 5];                % [pos; vel] for 2D double integrator
dim  = numel(x0);
m    = dim/2;                         % input dimension for double integrator
Tmax = 10;

% Discretization step for DT design/sim
dt     = 0.5;
Nsteps = round(Tmax/dt);

% LQR weights (penalize positions only in Q here)
Q = [eye(m) zeros(m); zeros(m) zeros(m)];
R = eye(m);

% Integrator
sim_integrator = 'ZOH';
model_integrator = 'ZOH';

%% Resolve dynamics (uniform interface)
dyn = resolve_dynamics('double_integrator');
[A,B] = dyn.getAB('DT', struct('dim',dim,'dt',dt,'dmode',model_integrator)); 
[Ac,Bc] = dyn.getAB('CT', struct('dim',dim));

%% Discrete-time LQR (design on DT model)
% (a) Finite-horizon Riccati recursion (terminal cost = I)
N = 100;
P = zeros(dim, dim, N); 
P(:,:,N) = eye(dim);
for i = N-1:-1:1
    S  = R + B' * P(:,:,i+1) * B;
    P(:,:,i) = Q + A' * (P(:,:,i+1) - P(:,:,i+1) * B / S * B' * P(:,:,i+1)) * A;
end
K_finite = (R + B' * P(:,:,1) * B) \ (B' * P(:,:,1) * A);
Pd_finite = P(:,:,1); 

% (b) Infinite-horizon (dlqr)
[K, Pd_dlqr] = dlqr(A, B, Q, R); 

% Plot Riccati entries vs. backward step
figure(1); clf; hold on;
labels = strings(0);
for r = 1:dim
    for c = 1:dim
        plot(permute(P(r,c,:), [3 2 1]), 'LineWidth', 1.0);
        labels(end+1) = sprintf('p_{%d%d}', r, c); 
    end
end
grid on; xlabel('Backward step i'); ylabel('Entry value');
title('Finite-horizon discrete Riccati recursion (P_i entries)');
legend(labels, 'Location', 'eastoutside');

fprintf('||K_finite - K_dlqr||_F = %.3e\n', norm(K_finite - K, 'fro'));

%% Continuous-time LQR (design on CT model)
[Kc, Pc, ~] = lqr(Ac, Bc, Q, R); 

%% Simulations with roll_out.m
% Continuous-time plant with CT LQR (u = -Kc x)
u_ct  = @(x, t) -Kc * x;
simCT = struct('type','CT','T',Tmax,'integrator','ode45');
[tlist_ct, xlist_ct, ulist_ct] = roll_out(dyn, u_ct, x0, simCT); 

% Discrete-time plant with DT LQR (x_{k+1} = A x_k + B u_k, u_k = -K x_k)
u_dt  = @(x, k) -K * x;  % DT controller
simDT = struct('type','DT','K',Nsteps,'dt',dt,'integrator',sim_integrator);
[tlist_dt, xlist_dt, ulist_dt] = roll_out(dyn, u_dt, x0, simDT);

%% Plot trajectories (position components) to compare CT vs DT LQR
figure(2); clf;
subplot(1,2,1); hold on;
plot(xlist_ct(1,:), xlist_ct(2,:), 'LineWidth', 2);      % CT
plot(xlist_dt(1,:), xlist_dt(2,:), 'LineWidth', 2);      % DT
axis equal; grid on; xlabel('x_1'); ylabel('x_2');
title('Position Trajectories');
legend('CT LQR (ode45)', ['DT LQR (',sim_integrator,], 'Location', 'best');

subplot(1,2,2); hold on;
plot(tlist_ct, xlist_ct.', 'LineWidth', 1.5);             % CT states
plot(tlist_dt, xlist_dt.', '.', 'MarkerSize', 8);         % DT states (markers)
grid on; xlabel('Time (s)'); ylabel('State');
title('States vs Time (CT solid, DT dots)');

%% Alternative implementation via synthesis()
sys.dt = dt; sys.A = A; sys.B = B; sys.Q = Q; sys.R = R;
c = synthesis('LQR', sys);  