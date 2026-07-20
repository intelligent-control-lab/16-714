function result = lecture7_lqr_spark_test(save_dir)
%LECTURE7_LQR_SPARK_TEST Export lecture 7 LQR traces for policy validation.

test_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(test_dir);
addpath(fullfile(repo_root, 'lib'));

if exist('OCTAVE_VERSION', 'builtin')
    pkg load control;
end

if nargin < 1 || isempty(save_dir)
    save_dir = fullfile(repo_root, 'data', 'lecture7');
end

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

dyn_name = 'double_integrator';
x0 = [10; 10; 1; 5];
dim = numel(x0);
m = dim / 2;
dt = 0.5;
tmax = 10;
nsteps = round(tmax / dt);
q = [eye(m) zeros(m); zeros(m) zeros(m)];
r = eye(m);
model_integrator = 'ZOH';
sim_integrator = 'ZOH';

dyn = resolve_dynamics(dyn_name);
[A, B] = dyn.getAB('DT', struct('dim', dim, 'dt', dt, 'dmode', model_integrator));
[K, P] = dlqr(A, B, q, r);

u_dt = @(x, t) -K * x;
sim = struct('type', 'DT', 'K', nsteps, 'dt', dt, 'integrator', sim_integrator);
[tlist, xlist, ulist] = roll_out(dyn, u_dt, x0, sim);

state_csv = fullfile(save_dir, 'states.csv');
control_csv = fullfile(save_dir, 'controls.csv');
write_csv_with_header(state_csv, {'t', 'x1', 'x2', 'v1', 'v2'}, [tlist(:), xlist.']);
write_csv_with_header(control_csv, {'t', 'u1', 'u2'}, [(0:nsteps-1).' * dt, ulist.']);

fig_path = fullfile(save_dir, 'lqr_state_plot.png');
save_lqr_plot(fig_path, tlist, xlist, ulist);

result = struct();
result.format_version = 1;
result.reference_name = 'lecture7_lqr_double_integrator';
result.model = dyn_name;
result.dt = dt;
result.horizon_steps = nsteps;
result.x0 = x0;
result.Q = q;
result.R = r;
result.K = K;
result.P = P;
result.final_state = xlist(:, end);
result.artifacts = struct( ...
    'directory', save_dir, ...
    'states_csv', state_csv, ...
    'controls_csv', control_csv, ...
    'plot_png', fig_path);

save(fullfile(save_dir, 'reference.mat'), 'result');
write_summary_json(fullfile(save_dir, 'summary.json'), result);

fprintf('Reference artifact directory: %s\n', save_dir);
fprintf('Final state norm: %.12g\n', norm(result.final_state));

end

function write_csv_with_header(path, header, data)
fid = fopen(path, 'w');
if fid < 0
    error('Could not open %s for writing.', path);
end
fprintf(fid, '%s\n', strjoin(header, ','));
fclose(fid);
dlmwrite(path, data, '-append', 'delimiter', ',', 'precision', '%.16g');
end

function save_lqr_plot(path, tlist, xlist, ulist)
fig = figure('Visible', 'off', 'Position', [100, 100, 760, 460]);

subplot(2, 2, [1 3]);
plot(xlist(1, :), xlist(2, :), 'LineWidth', 2);
hold on;
plot(xlist(1, 1), xlist(2, 1), 'ko', 'MarkerSize', 4);
plot(xlist(1, end), xlist(2, end), 'ro', 'MarkerSize', 4);
axis equal;
grid on;
xlabel('x1');
ylabel('x2');
title('Position Trajectory');

subplot(2, 2, 2);
plot(tlist, xlist.', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('state');
title('State');
legend({'x1', 'x2', 'v1', 'v2'}, 'Location', 'best');

subplot(2, 2, 4);
plot(tlist(1:end-1), ulist.', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('control');
title('Control');
legend({'u1', 'u2'}, 'Location', 'best');

try
    exportgraphics(fig, path, 'Resolution', 160);
catch
    saveas(fig, path);
end
close(fig);
end

function write_summary_json(path, result)
payload = struct();
payload.format_version = result.format_version;
payload.reference_name = result.reference_name;
payload.model = result.model;
payload.dt = result.dt;
payload.horizon_steps = result.horizon_steps;
payload.final_state_norm = norm(result.final_state);
payload.created_at = datestr(now, 31);
payload.version = version;
payload.artifacts = result.artifacts;

try
    encoded = jsonencode(payload, 'PrettyPrint', true);
catch
    encoded = jsonencode(payload);
end

fid = fopen(path, 'w');
if fid < 0
    error('Could not open %s for writing.', path);
end
fprintf(fid, '%s\n', encoded);
fclose(fid);
end
