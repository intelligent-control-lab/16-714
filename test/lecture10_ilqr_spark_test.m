function result = lecture10_ilqr_spark_test(save_dir)
%LECTURE10_ILQR_SPARK_TEST Export MATLAB iLQR traces for SPARK validation.
%
% This file intentionally keeps the SPARK validation wrapper outside /lib.
% The local iLQR routine below mirrors the lecture iLQR algorithm but returns
% machine-readable traces and metadata for comparison with spark_version.

test_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(test_dir);
addpath(fullfile(repo_root, 'lib'));

if nargin < 1 || isempty(save_dir)
    save_dir = fullfile(repo_root, 'data');
end

spec = make_reference_spec(save_dir);
result = run_reference(spec);

fprintf('Reference artifact directory: %s\n', result.artifacts.directory);
fprintf('Final cost: %.12g\n', result.cost);
fprintf('Goal error norm: %.12g\n', norm(result.goal_error));

end

function spec = make_reference_spec(save_dir)
opts = struct();
opts.N = 150;
opts.dt = 0.1;
opts.integrator = 'Euler';
opts.initial_nominal = 'rollout';
opts.Q = diag([1.0, 1.0, 0.1]);
opts.R = diag([2.0, 2.5]);
opts.S = diag([100.0, 100.0, 50.0]);
opts.u_ref = zeros(2, 1);
opts.alphas = [1.0, 0.5, 0.25, 0.1, 0.05];
opts.maxIter = 120;
opts.tolJ = 1e-8;
opts.lambda0 = 1e-8;
opts.lambdaMax = 1e8;
opts.plot = false;
opts.verbose = true;

comparison = struct();
comparison.purpose = 'SPARK Agibot mobile-base iLQR parity check';
comparison.matlab_state_names = {'x', 'y', 'theta'};
comparison.matlab_control_names = {'v', 'omega'};
comparison.spark_robot_cfg = 'AgiBotG1MobileBaseDynamic1Config';
comparison.spark_agent_cfg = 'AgiBotG1MobileBaseAgent';
comparison.spark_use_sim_dynamics = false;
comparison.spark_state_names = {'LinearX', 'LinearY', 'RotYaw'};
comparison.spark_control_names = {'vLinearX', 'vRotYaw'};
comparison.spark_zero_control_names = {'vLinearY'};
comparison.spark_dt = opts.dt;
comparison.spark_horizon_steps = opts.N;
comparison.spark_control_limit_guard = 0.3;

spec = struct();
spec.name = 'ilqr_unicycle3_agibot_base_small';
spec.description = ['Unicycle3 iLQR reference for the first SPARK Agibot ', ...
                    'mobile-base implementation.'];
spec.model = 'unicycle3';
spec.x0 = [0.0; 0.0; 0.0];
spec.xg = [0.5; 0.25; 0.3];
spec.opts = opts;
spec.save_dir = save_dir;
spec.artifact_dir = save_dir;
spec.comparison = comparison;
end

function result = run_reference(spec)
result = local_ilqr(spec.model, spec.x0, spec.xg, spec.opts);
result.format_version = 1;
result.reference_name = spec.name;
result.reference_module = 'ilqr';
result.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
result.description = spec.description;
result.comparison = spec.comparison;
result.source = struct( ...
    'repo_root', fileparts(mfilename('fullpath')), ...
    'pipeline', 'lecture10_ilqr_spark_test', ...
    'matlab_version', version);
result.artifacts = struct();
result = save_reference_artifacts(result, spec);
end

function result = local_ilqr(model, x0, xg, opts)
N = get_opt(opts, 'N', 100);
dt = get_opt(opts, 'dt', 0.1);
Q = get_opt(opts, 'Q', diag([1, 1, 0.1]));
R = get_opt(opts, 'R', diag([1, 0.5]));
S = get_opt(opts, 'S', diag([50, 50, 10]));
u_ref = get_opt(opts, 'u_ref', zeros(2, 1));
alphas = get_opt(opts, 'alphas', [1.0 0.5 0.25 0.1 0.05]);
maxIter = get_opt(opts, 'maxIter', 100);
tolJ = get_opt(opts, 'tolJ', 1e-1);
lambda0 = get_opt(opts, 'lambda0', 1e-8);
lambdaMax = get_opt(opts, 'lambdaMax', 1e8);
integr = get_opt(opts, 'integrator', 'ZOH');
initial_nominal = get_opt(opts, 'initial_nominal', 'linear');
plotFlag = get_opt(opts, 'plot', true);
verbose = get_opt(opts, 'verbose', true);
improve_tol = 1e-12;

dyn = resolve_dynamics(model);
nu = numel(u_ref);
ubar = zeros(nu, N);
nx = numel(x0);
xbar = zeros(nx, N + 1);
for k = 0:N
    tau = k / N;
    xbar(:, k + 1) = (1 - tau) * x0 + tau * xg;
end

if strcmpi(initial_nominal, 'rollout') || strcmpi(initial_nominal, 'dynamic')
    ctrl_replay = @(x, t) replay_u(t, ubar, dt);
    simDT = struct('type', 'DT', 'K', N, 'dt', dt, 'integrator', integr);
    [~, xbar, ubar_] = roll_out(model, ctrl_replay, x0, simDT);
    ubar = ubar_;
elseif ~strcmpi(initial_nominal, 'linear')
    error('Unknown initial_nominal "%s". Use linear or rollout.', initial_nominal);
end

J = total_cost(Q, R, S, xbar, ubar, xg, u_ref);
if verbose
    fprintf('Iter %2d: J = %.6f\n', 0, J);
end
prevJ = J;

lambda = lambda0;
Jhist = J;
iter_hist = 0;
alpha_hist = [];
lambda_hist = lambda;
accepted_iter = 0;
backward_failures = 0;
line_search_rejections = 0;
last_it = 0;
status = 'max_iter';

for it = 1:maxIter
    last_it = it;
    A = cell(1, N);
    B = cell(1, N);
    d = cell(1, N);
    q = cell(1, N);
    r = cell(1, N);
    for k = 1:N
        xk = xbar(:, k);
        uk = ubar(:, k);
        [A{k}, B{k}] = dyn.getAB('DT', struct('dim', nx, 'dt', dt, 'dmode', 'ZOH', 'x', xk, 'u', uk));
        xk_plus = step(xk, uk, dt, dyn, 'Euler');
        d{k} = xk_plus - xbar(:, k + 1);
        q{k} = Q * (xbar(:, k) - xg);
        r{k} = R * (uk - u_ref);
    end

    P = cell(1, N + 1);
    s = cell(1, N + 1);
    K = cell(1, N);
    kfeed = cell(1, N);
    P{N + 1} = S;
    s{N + 1} = zeros(nx, 1);
    diverged = false;

    for k = N:-1:1
        Ak = A{k};
        Bk = B{k};
        dk = d{k};
        Pn = P{k + 1};
        sn = s{k + 1};
        qk = q{k};
        rk = r{k};

        Quu = R + Bk.' * Pn * Bk + lambda * eye(size(R));
        Qux = Bk.' * Pn * Ak;
        gu = rk + Bk.' * (Pn * dk + sn);

        [L, p] = chol(Quu, 'lower');
        if p > 0
            diverged = true;
            break;
        end
        invQuu = L' \ (L \ eye(size(Quu)));

        K{k} = -invQuu * Qux;
        kfeed{k} = -invQuu * gu;

        P{k} = Q + Ak.' * Pn * Ak - Qux.' * invQuu * Qux;
        s{k} = qk + Ak.' * (Pn * dk + sn) - Qux.' * invQuu * gu;
        P{k} = 0.5 * (P{k} + P{k}.');
    end

    if diverged
        backward_failures = backward_failures + 1;
        lambda = min(lambda * 10, lambdaMax);
        if lambda >= lambdaMax
            status = 'backward_pass_failed';
            if verbose
                warning('Backward pass failed (Quu not PD). Stopping.');
            end
            break;
        end
        continue;
    end

    simDT = struct('type', 'DT', 'K', N, 'dt', dt, 'integrator', integr);
    best = struct('J', inf, 'alpha', NaN, 'x', [], 'u', []);
    for a = alphas
        ctrl_affine = @(x, t) affine_ctrl(t, x, xbar, ubar, kfeed, K, dt, a);
        [~, xnew, unew] = roll_out(model, ctrl_affine, x0, simDT);
        Jcand = total_cost(Q, R, S, xnew, unew, xg, u_ref);
        if Jcand < best.J
            best.J = Jcand;
            best.alpha = a;
            best.x = xnew;
            best.u = unew;
        end
    end

    if (J - best.J) > max(improve_tol, 1e-12 * abs(J))
        xbar = best.x;
        ubar = best.u;
        J = best.J;
        lambda = max(lambda / 5, 1e-12);
        accepted_iter = accepted_iter + 1;
        Jhist(end + 1) = J; %#ok<AGROW>
        iter_hist(end + 1) = it; %#ok<AGROW>
        alpha_hist(end + 1) = best.alpha; %#ok<AGROW>
        lambda_hist(end + 1) = lambda; %#ok<AGROW>
    else
        line_search_rejections = line_search_rejections + 1;
        lambda = min(lambda * 10, lambdaMax);
        if lambda >= lambdaMax
            status = 'line_search_failed';
            if verbose
                warning('No improvement from line search. Stopping.');
            end
            break;
        end
        continue;
    end

    if verbose
        fprintf('Iter %2d: J = %.6f  (alpha=%.2f, lambda=%.1e)\n', it, J, best.alpha, lambda);
    end
    if abs(prevJ - J) < tolJ
        status = 'tolJ';
        break;
    end
    prevJ = J;
end

if verbose
    fprintf('Final state: [%s]^T\n', num2str(xbar(:, end).', '%g '));
    fprintf('Goal error:  [%s]^T\n', num2str((xbar(:, end) - xg).', '%g '));
end

if plotFlag
    figure;
    plot(xbar(1, :), xbar(2, :), 'LineWidth', 2);
    hold on;
    plot(x0(1), x0(2), 'ko', 'MarkerFaceColor', 'k');
    plot(xg(1), xg(2), 'r*', 'MarkerSize', 10);
    axis equal;
    grid on;
    xlabel('x');
    ylabel('y');
    title(sprintf('%s iLQR trajectory', model));
end

result = struct();
result.model = char(model);
result.x0 = x0(:);
result.xg = xg(:);
result.N = N;
result.dt = dt;
result.integrator = char(integr);
result.initial_nominal = char(initial_nominal);
result.Q = Q;
result.R = R;
result.S = S;
result.u_ref = u_ref(:);
result.x = xbar;
result.u = ubar;
result.t = 0:dt:(N * dt);
result.cost = J;
result.cost_history = Jhist;
result.cost_iteration = iter_hist;
result.alpha_history = alpha_hist;
result.lambda_history = lambda_hist;
result.iterations = last_it;
result.accepted_iterations = accepted_iter;
result.backward_failures = backward_failures;
result.line_search_rejections = line_search_rejections;
result.status = status;
result.final_state = xbar(:, end);
result.goal_error = xbar(:, end) - xg(:);
result.opts = opts;
end

function result = save_reference_artifacts(result, spec)
artifact_dir = spec.artifact_dir;
if ~exist(artifact_dir, 'dir')
    mkdir(artifact_dir);
end

artifacts = struct();
artifacts.directory = artifact_dir;
artifacts.mat = fullfile(artifact_dir, 'reference.mat');
artifacts.states_csv = fullfile(artifact_dir, 'states.csv');
artifacts.controls_csv = fullfile(artifact_dir, 'controls.csv');
artifacts.state_plot_png = fullfile(artifact_dir, 'reference_ilqr_state_plot.png');
artifacts.summary_json = fullfile(artifact_dir, 'summary.json');
artifacts.reference_json = fullfile(artifact_dir, 'reference.json');
result.artifacts = artifacts;

save(artifacts.mat, 'result', '-v7');
write_trajectory_csv(artifacts.states_csv, result.t(:), result.x.', {'x', 'y', 'theta'});
write_trajectory_csv(artifacts.controls_csv, result.t(1:end-1).', result.u.', {'v', 'omega'});
save_state_plot(result, artifacts.state_plot_png);
write_json(artifacts.summary_json, make_summary(result));

try
    write_json(artifacts.reference_json, result);
catch exc
    warning('lecture10_ilqr_spark_test:json', ...
        'Could not write full reference.json: %s', exc.message);
    result.artifacts.reference_json = '';
    save(artifacts.mat, 'result', '-v7');
end
end

function write_trajectory_csv(path, t, values, names)
table_values = array2table([t, values], 'VariableNames', [{'t'}, names]);
writetable(table_values, path);
end

function save_state_plot(result, output_path)
x = result.x;
u = result.u;
t = result.t(:).';
t_u = t(1:size(u, 2));

fig = figure('Visible', 'off', 'Position', [100 100 760 460]);

ax_path = subplot(2, 2, [1 3], 'Parent', fig);
hold(ax_path, 'on');
grid(ax_path, 'on');
box(ax_path, 'on');
plot(ax_path, x(1, :), x(2, :), 'Color', [0.0000 0.4470 0.7410], 'LineWidth', 2, ...
    'DisplayName', 'trajectory');
plot(ax_path, x(1, 1), x(2, 1), 'ko', 'MarkerSize', 4, 'DisplayName', 'start');
plot(ax_path, x(1, end), x(2, end), 'o', 'Color', [0.8500 0.3250 0.0980], ...
    'MarkerSize', 4, 'DisplayName', 'final');
plot(ax_path, result.xg(1), result.xg(2), 'r*', 'MarkerSize', 8, 'DisplayName', 'goal');
axis(ax_path, 'equal');
title(ax_path, 'Base Path');
xlabel(ax_path, 'x');
ylabel(ax_path, 'y');
legend(ax_path, 'Location', 'best');

ax_state = subplot(2, 2, 2, 'Parent', fig);
hold(ax_state, 'on');
grid(ax_state, 'on');
box(ax_state, 'on');
state_names = {'x', 'y', 'theta'};
state_colors = lines(3);
for i = 1:3
    plot(ax_state, t, x(i, :), 'LineWidth', 1.5, 'Color', state_colors(i, :), ...
        'DisplayName', state_names{i});
end
title(ax_state, 'Base State');
xlabel(ax_state, 't [s]');
ylabel(ax_state, 'state');
legend(ax_state, 'Location', 'best');

ax_control = subplot(2, 2, 4, 'Parent', fig);
hold(ax_control, 'on');
grid(ax_control, 'on');
box(ax_control, 'on');
control_names = {'v', 'omega'};
control_colors = lines(2);
for i = 1:2
    plot(ax_control, t_u, u(i, :), 'LineWidth', 1.5, 'Color', control_colors(i, :), ...
        'DisplayName', control_names{i});
end
title(ax_control, 'Base Control');
xlabel(ax_control, 't [s]');
ylabel(ax_control, 'control');
legend(ax_control, 'Location', 'best');

sgtitle(fig, 'iLQR Reference');

try
    exportgraphics(fig, output_path, 'Resolution', 160);
catch
    print(fig, output_path, '-dpng', '-r160');
end
close(fig);
end

function summary = make_summary(result)
summary = struct();
summary.format_version = result.format_version;
summary.reference_name = result.reference_name;
summary.reference_module = result.reference_module;
summary.created_at = result.created_at;
summary.model = result.model;
summary.N = result.N;
summary.dt = result.dt;
summary.integrator = result.integrator;
summary.initial_nominal = result.initial_nominal;
summary.status = result.status;
summary.iterations = result.iterations;
summary.accepted_iterations = result.accepted_iterations;
summary.cost = result.cost;
summary.cost_history = result.cost_history;
summary.final_state = result.final_state;
summary.goal_error = result.goal_error;
summary.max_abs_u = max(abs(result.u), [], 2);
summary.source = result.source;
summary.comparison = result.comparison;
summary.artifacts = result.artifacts;
end

function write_json(path, payload)
try
    encoded = jsonencode(payload, 'PrettyPrint', true);
catch
    encoded = jsonencode(payload);
end
fid = fopen(path, 'w');
if fid < 0
    error('lecture10_ilqr_spark_test:write_json', 'Could not open %s for writing.', path);
end
cleaner = onCleanup(@() fclose(fid));
fwrite(fid, encoded, 'char');
fwrite(fid, newline, 'char');
clear cleaner;
end

function val = get_opt(s, field, default)
if isfield(s, field)
    val = s.(field);
else
    val = default;
end
end

function u = replay_u(t, U, dt)
k = 1 + floor(t / dt);
k = max(1, min(size(U, 2), k));
u = U(:, k);
end

function u = affine_ctrl(t, x, xbar, ubar, kfeed, K, dt, alpha)
k = 1 + floor(t / dt);
k = max(1, min(size(ubar, 2), k));
u = ubar(:, k) + alpha * kfeed{k} + K{k} * (x - xbar(:, k));
end

function J = total_cost(Q, R, S, x, u, xg, uref)
Nloc = size(u, 2);
J = 0;
for kk = 1:Nloc
    dx = x(:, kk) - xg;
    du = u(:, kk) - uref;
    J = J + dx.' * Q * dx + du.' * R * du;
end
dT = x(:, Nloc + 1) - xg;
J = J + dT.' * S * dT;
end
