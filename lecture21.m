% 16-714 Advanced Control for Robotics
% Script for Lecture 21: Value Approximation

%% Problem Setup
sys.A = 1; sys.B = 0.5;
sys.f_dt = @(x,u) sys.A*x + sys.B*u;   % DT dynamics
sys.Q = 1; sys.R = 1;
sys.h = @(x,u) (x'*sys.Q*x + u'*sys.R*u)/2;   % stage cost
sys.tolerance = 0.01;
sys.kmax = 100; sys.dt = 1; sys.x0 = 1;
sys.tcondition = @(x,t) norm(x) < sys.tolerance || t > sys.kmax * sys.dt;

% roll_out simulation config (DT)
sim = struct( ...
    'type','DT', ...
    'K', sys.kmax, ...
    'dt', sys.dt, ...
    'integrator','Direct', ...
    'stop', @(xk,tk,log) sys.tcondition(xk,tk) ...
);
opts = struct();   % keep empty unless you need extras

%% LQR
c = synthesis('LQR', sys);

[~, x_lqr] = roll_out(sys, @(x,t)c.u(x,t), sys.x0, sim, opts);

W_gt = [sys.Q+sys.A'*c.P*sys.A,  sys.A'*c.P*sys.B;
        sys.B'*c.P*sys.A,        sys.B'*c.P*sys.B + sys.R];

%% Learning parameters
sys.W0 = [4 1; 1 4];    % initial quadratic weights
sys.delta = 1;          % no discount
sys.epsilon = 0.1;      % epsilon-greedy
n_ep = 100;
alpha = 1;

%% Monte Carlo
rms_mc = zeros(1,n_ep);
figure(1); clf; subplot(311); hold on;
mc = RLagent(sys);
for episode = 1:n_ep
    x_list = mc.update('MC', alpha, sys);    % uses new roll_out inside
    plot(0:length(x_list)-1, x_list, 'color', [0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep])
    rms_mc(episode) = norm(mc.W - W_gt);
end
title("MC"); box on;

%% SARSA
rms_sarsa = zeros(1,n_ep);
subplot(312); hold on;
sarsa = RLagent(sys);
for episode = 1:n_ep
    x_list = sarsa.update('SARSA', alpha, sys);  % uses sys.f internally
    plot(0:length(x_list)-1, x_list, 'color', [0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep])
    rms_sarsa(episode) = norm(sarsa.W - W_gt);
end
title("Sarsa"); box on;

%% Q-Learning
rms_q = zeros(1,n_ep);
subplot(313); hold on;
qlearning = RLagent(sys);
for episode = 1:n_ep
    x_list = qlearning.update('QLearning', alpha, sys);  % uses sys.f internally
    plot(0:length(x_list)-1, x_list, 'color', [0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep])
    rms_q(episode) = norm(qlearning.W - W_gt);
end
title("Q-learning"); box on;

%% Visualize SARSA vs Q-Learning vs MC
figure; hold on;
plot(rms_sarsa); plot(rms_q); plot(rms_mc);
legend("SARSA", "Q-Learning","Monte Carlo"); box on;
