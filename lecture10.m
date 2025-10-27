% 16-714 Advanced Control for Robotics
% Script for Lecture 10: MPC Feasibility (strict bounds)
% Requires: resolve_dynamics.m, roll_out.m, synthesis.m

clear; close all; 

%% ========= Problem spec (define bounds in canonical vector form!) =========
problem = 2;  % 1: single_integrator, 2: double_integrator
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

        % Bounds — MUST match nx=1, nu=1 exactly (column vectors)
        sys.umax = [ 1];
        sys.umin = [-1];
        sys.xmax = [10];
        sys.xmin = [-1];

    case 2
        sys.dt   = 0.5;
        sys.name = 'double_integrator';
        sys.Q    = [1 0; 0 0];
        sys.R    = 1;
        sys.S    = 10*eye(2);
        sys.x0   = [10; 0];
        sys.N    = 10;

        % Bounds — MUST match nx=2, nu=1 exactly (column vectors)
        sys.umax = [ 1];
        sys.umin = [-1];
        sys.xmax = [ 10;  5];
        sys.xmin = [-10; -2];

    otherwise
        error('Unknown problem id.');
end

%% ========= Dynamics via resolve_dynamics (DT, ZOH) =========
dyn  = resolve_dynamics(sys.name);
nx   = length(sys.x0);
[sys.A, sys.B] = dyn.getAB('DT', struct('dim', nx, 'dt', sys.dt, 'dmode', 'ZOH'));
nu   = size(sys.B, 2);

%% ========= Controllers =========
cMPCn = synthesis('nMPC', sys);   % unconstrained nominal (reference)
cMPCu = synthesis('uMPC', sys);   % control-constrained
cMPCq = synthesis('qMPC', sys);   % state+control constrained

%% ========= Simulation =========
simhorizon = 3*sys.N;
simDT = struct('type','DT','K',simhorizon,'dt',sys.dt,'integrator','ZOH');

[tk, x_n, u_n] = roll_out(sys.name, cMPCn.u, sys.x0, simDT); 
[~,  x_u, u_u] = roll_out(sys.name, cMPCu.u, sys.x0, simDT);
[~,  x_q, u_q] = roll_out(sys.name, cMPCq.u, sys.x0, simDT);

tx = 0:simhorizon;       % states have K+1 samples
tu = 0:simhorizon-1;     % inputs have K samples

%% ========= Plot: state trajectories with bounds (no broadcasting) =========
figure('Name','MPC Feasibility - States'); clf;
tl = tiledlayout(nx,1,'TileSpacing','compact','Padding','compact');

for i = 1:nx
    ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(tx, x_n(i,:), 'DisplayName','Unconstrained');
    plot(tx, x_u(i,:), 'DisplayName','Control Constrained');
    plot(tx, x_q(i,:), 'DisplayName','State & Control Constrained');
    yline(0,'--k','HandleVisibility','off');

    % state bounds (as provided, exact size nx×1)
    plot(tx, sys.xmax(i)*ones(size(tx)), ':k', 'DisplayName','x_{max}');
    plot(tx, sys.xmin(i)*ones(size(tx)), ':k', 'DisplayName','x_{min}');

    ylabel(sprintf('x_%d',i));
    if i == 1, title(tl, 'State Trajectories with Constraints'); end
    if i == nx, xlabel('k'); legend('Location','best'); end
end

%% ========= Plot: control trajectories with bounds (no broadcasting) =========
figure('Name','MPC Feasibility - Inputs'); clf;
tl2 = tiledlayout(nu,1,'TileSpacing','compact','Padding','compact');

for j = 1:nu
    ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(tu, u_n(j,:), 'DisplayName','Unconstrained');
    plot(tu, u_u(j,:), 'DisplayName','Control Constrained');
    plot(tu, u_q(j,:), 'DisplayName','State & Control Constrained');

    % input bounds (as provided, exact size nu×1)
    plot(tu, sys.umax(j)*ones(size(tu)), ':k', 'DisplayName','u_{max}');
    plot(tu, sys.umin(j)*ones(size(tu)), ':k', 'DisplayName','u_{min}');

    ylabel(sprintf('u_%d',j));
    if j == 1, title(tl2, 'Control Trajectories with Constraints'); end
    if j == nu, xlabel('k'); legend('Location','best'); end
end

%% ========= Feasibility visualization (optional) =========
feasible = visualize_mpc(sys, cMPCq, simhorizon); 
