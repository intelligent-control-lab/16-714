% 16-714 Advanced Control for Robotics
% Script for Lecture 9: MPC
% - Uses resolve_dynamics.m for dynamics/Jacobians:
%     dyn = resolve_dynamics(sys.name);
%     [A,B] = dyn.getAB('DT', struct('dim',dim,'dt',dt));
% - Uses the unified roll_out.m interface:
%     [tlist, xlist, ulist] = roll_out(name, ctrl, x0, simStruct)
% - Assumes synthesis('LQR'|'LQRn'|'nMPC', sys) returns controllers with field .u

clear; close all;

%% ========== Problem Spec ==========
problem = 2;  % 1: single-integrator, 2: double-integrator
sys = struct();

switch problem
    case 1
        sys.dt   = 0.5;
        sys.name = 'single_integrator';
        sys.Q    = 1;
        sys.R    = 10;
        sys.S    = 100;
        sys.x0   = [10];
        sys.N    = 10;

    case 2
        sys.dt   = 0.5;
        sys.name = 'double_integrator';
        sys.Q    = [1 0; 0 0];
        sys.R    = 1;
        sys.S    = 10*eye(2);
        sys.x0   = [10; 0];
        sys.N    = 10;

    otherwise
        error('Unknown problem id.');
end

% ---- Use resolve_dynamics to obtain (A,B) ----
dyn = resolve_dynamics(sys.name);              % must be on path
dim = length(sys.x0);                          % infer state dimension from x0
[sys.A, sys.B] = dyn.getAB('DT', struct('dim', dim, 'dt', sys.dt));
% (Optional) If your resolve_dynamics supports CT as well:
% [A_ct, B_ct] = dyn.getAB('CT', struct('dim', dim));

%% ========== Controllers ==========
cLQR  = synthesis('LQR',  sys);   % infinite-horizon LQR
cLQRn = synthesis('LQRn', sys);   % finite-horizon LQR (length = sys.N)
cMPC  = synthesis('nMPC', sys);   % (nominal) MPC with horizon sys.N

%% ========== Compare the P matrices ==========
figure('Name','P-matrix comparison'); hold on; box on; grid on;
leg = {};
% MPC (time-varying P_k)
for i = 1:size(cMPC.P,1)
    for j = 1:size(cMPC.P,2)
        plot(0:sys.N, permute(cMPC.P(i,j,1:sys.N+1), [3,1,2]), ...
            '-*', 'DisplayName', sprintf('P_{%d%d}(MPC)', i, j));
        leg{end+1} = sprintf('P_{%d%d}(MPC)', i, j); 
    end
end
% Infinite-horizon LQR (constant P)
for i = 1:size(cLQR.P,1)
    for j = 1:size(cLQR.P,2)
        plot(0:sys.N, cLQR.P(i,j)*ones(1,sys.N+1), '--', ...
            'DisplayName', sprintf('P_{%d%d}(LQR∞)', i, j));
        leg{end+1} = sprintf('P_{%d%d}(LQR∞)', i, j); 
    end
end
xlabel('k'); ylabel('P entries'); legend(leg, 'Location','best');

%% ========== Compare the gain K ==========
figure('Name','K-gain comparison'); hold on; box on; grid on;
leg = {};
% MPC (time-varying K_k)
for i = 1:size(cMPC.K,1)
    for j = 1:size(cMPC.K,2)
        plot(0:sys.N-1, permute(cMPC.K(i,j,1:sys.N), [3,1,2]), ...
            '-*', 'DisplayName', sprintf('K_{%d%d}(MPC)', i, j));
        leg{end+1} = sprintf('K_{%d%d}(MPC)', i, j); 
    end
end
% Infinite-horizon LQR (constant K)
for i = 1:size(cLQR.K,1)
    for j = 1:size(cLQR.K,2)
        plot(0:sys.N-1, cLQR.K(i,j)*ones(1,sys.N), '--', ...
            'DisplayName', sprintf('K_{%d%d}(LQR∞)', i, j));
        leg{end+1} = sprintf('K_{%d%d}(LQR∞)', i, j); 
    end
end
xlabel('k'); ylabel('K entries'); legend(leg, 'Location','best');

%% ========== Visualize trajectories (finite LQR vs infinite LQR vs MPC) ==========
simDT = struct('type','DT', 'K', sys.N, 'dt', sys.dt, 'integrator', 'ZOH');

[~, x_f,   ~] = roll_out(sys.name, cLQRn.u, sys.x0, simDT);
[~, x_inf, ~] = roll_out(sys.name, cLQR.u,  sys.x0, simDT);
[~, x_mpc, ~] = roll_out(sys.name, cMPC.u,  sys.x0, simDT);

figure('Name','State trajectory comparison'); hold on; box on; grid on;
k = 0:sys.N;
plot(k, x_f(1,:),   'DisplayName','Finite LQR');
plot(k, x_inf(1,:), 'DisplayName','Infinite LQR');
plot(k, x_mpc(1,:), 'DisplayName','MPC');
yline(0, '--k', 'DisplayName','reference');
xlabel('k'); ylabel('x_1'); legend('Location','best');

%% ========== MPC with Different Preview Horizons ==========
N0 = sys.N;  % keep original horizon
controllers = cell(1, N0);
for i = 1:N0
    sys_i = sys; sys_i.N = i;
    controllers{i} = synthesis('nMPC', sys_i);
end

% simulate over 2*N0 horizon
simDT_long = struct('type','DT', 'K', 2*N0, 'dt', sys.dt, 'integrator', 'ZOH');
[~, x_inf_long, ~] = roll_out(sys.name, cLQR.u, sys.x0, simDT_long);

x_mpc_h = zeros(length(sys.x0), 2*N0+1, N0);
for i = 1:N0
    [~, x_mpc_h(:,:,i), ~] = roll_out(sys.name, controllers{i}.u, sys.x0, simDT_long);
end

figure('Name','MPC horizon sweep'); hold on; box on; grid on;
plot(0:N0,   x_f(1,:),          'DisplayName','Finite LQR');
plot(0:2*N0, x_inf_long(1,:),   'DisplayName','Infinite LQR');
for i = 1:N0
    plot(0:2*N0, x_mpc_h(1,:,i), 'DisplayName', sprintf('MPC-%d', i),"color", [i/N0, 1-i/N0, i/N0]);
end
yline(0, '--k', 'DisplayName','reference');
xlabel('k'); ylabel('x_1'); legend('Location','best');

%% ========== Predicted Output at Each MPC Step ==========
% At every time step, a finite-horizon LQR of length sys.N is effectively solved.
x_pred_steps = zeros(length(sys.x0), N0+1, N0+1);

% First prediction from the initial state
[~, x_pred_steps(:,:,1), ~] = roll_out(sys.name, cLQRn.u, sys.x0, simDT);

% Receding-horizon predictions seeded by previous step's next state
for t = 1:N0
    x_seed = x_pred_steps(:, 2, t);  % one-step-ahead state from previous plan
    [~, x_pred_steps(:,:,t+1), ~] = roll_out(sys.name, cLQRn.u, x_seed, simDT);
end

figure('Name','Predicted trajectories per MPC step'); hold on; box on; grid on;
plot(0:N0,   x_f(1,:),        'DisplayName','Finite LQR');
plot(0:2*N0, x_inf_long(1,:), 'DisplayName','Infinite LQR');
for t = 1:N0
    plot((t-1):(t+N0-1), x_pred_steps(1,:,t), ...
        'DisplayName', sprintf('MPC Time %d', t-1), ...
        "color", [1-t/N0, t/N0, 1-t/N0]);
end
xlabel('k'); ylabel('x_1'); legend('Location','best');
