function result = lecture10_mpc_spark_test(save_dir)
%LECTURE10_MPC_SPARK_TEST Export lecture 10 MPC traces for policy validation.

test_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(test_dir);
addpath(fullfile(repo_root, 'lib'));

if exist('OCTAVE_VERSION', 'builtin')
    pkg load control;
    pkg load optim;
end

if nargin < 1 || isempty(save_dir)
    save_dir = fullfile(repo_root, 'data', 'lecture10_mpc');
end

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

sys = struct();
sys.dt = 0.5;
sys.name = 'double_integrator';
sys.Q = [1 0; 0 0];
sys.R = 1;
sys.S = 10 * eye(2);
sys.x0 = [10; 0];
sys.N = 10;
sys.umax = [1];
sys.umin = [-1];
sys.xmax = [10; 5];
sys.xmin = [-10; -2];

dyn = resolve_dynamics(sys.name);
[sys.A, sys.B] = dyn.getAB('DT', struct('dim', numel(sys.x0), 'dt', sys.dt, 'dmode', 'ZOH'));
simhorizon = 3 * sys.N;
sim = struct('type', 'DT', 'K', simhorizon, 'dt', sys.dt, 'integrator', 'ZOH');

cases = {
    'linear_mpc', 'nMPC';
    'input_constrained_mpc', 'uMPC';
    'constrained_mpc', 'qMPC';
};

result = struct();
result.format_version = 1;
result.reference_name = 'lecture10_mpc_double_integrator';
result.model = sys.name;
result.dt = sys.dt;
result.horizon_steps = sys.N;
result.execution_steps = simhorizon;
result.cases = struct();

for i = 1:size(cases, 1)
    case_name = cases{i, 1};
    synthesis_name = cases{i, 2};
    controller = synthesis(synthesis_name, sys);
    [tlist, xlist, ulist] = roll_out(sys.name, controller.u, sys.x0, sim);

    case_dir = fullfile(save_dir, case_name);
    if ~exist(case_dir, 'dir')
        mkdir(case_dir);
    end
    write_csv_with_header(fullfile(case_dir, 'states.csv'), {'t', 'x1', 'v1'}, [tlist(:), xlist.']);
    write_csv_with_header(fullfile(case_dir, 'controls.csv'), {'t', 'u1'}, [(0:simhorizon-1).' * sys.dt, ulist.']);

    case_result = struct();
    case_result.synthesis_name = synthesis_name;
    case_result.final_state = xlist(:, end);
    case_result.artifacts = struct( ...
        'directory', case_dir, ...
        'states_csv', fullfile(case_dir, 'states.csv'), ...
        'controls_csv', fullfile(case_dir, 'controls.csv'));
    result.cases.(case_name) = case_result;
end

save_summary_plot(fullfile(save_dir, 'lecture10_mpc_state_plot.png'), save_dir, cases, sys);
save(fullfile(save_dir, 'reference.mat'), 'result');
write_summary_json(fullfile(save_dir, 'summary.json'), result, cases);

fprintf('Reference artifact directory: %s\n', save_dir);
for i = 1:size(cases, 1)
    case_name = cases{i, 1};
    fprintf('%s final state norm: %.12g\n', case_name, norm(result.cases.(case_name).final_state));
end

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

function save_summary_plot(path, save_dir, cases, sys)
fig = figure('Visible', 'off', 'Position', [100, 100, 760, 460]);

subplot(2, 1, 1);
hold on;
grid on;
box on;
for i = 1:size(cases, 1)
    case_name = cases{i, 1};
    states = csvread(fullfile(save_dir, case_name, 'states.csv'), 1, 0);
    plot(states(:, 1), states(:, 2), 'LineWidth', 1.5, 'DisplayName', case_name);
end
plot([0, 3 * sys.N * sys.dt], [sys.xmax(1), sys.xmax(1)], ':k', 'DisplayName', 'x max');
plot([0, 3 * sys.N * sys.dt], [sys.xmin(1), sys.xmin(1)], ':k', 'DisplayName', 'x min');
xlabel('t [s]');
ylabel('x1');
title('Lecture 10 MPC State');
legend('Location', 'best');

subplot(2, 1, 2);
hold on;
grid on;
box on;
for i = 1:size(cases, 1)
    case_name = cases{i, 1};
    controls = csvread(fullfile(save_dir, case_name, 'controls.csv'), 1, 0);
    plot(controls(:, 1), controls(:, 2), 'LineWidth', 1.5, 'DisplayName', case_name);
end
plot([0, (3 * sys.N - 1) * sys.dt], [sys.umax(1), sys.umax(1)], ':k', 'DisplayName', 'u max');
plot([0, (3 * sys.N - 1) * sys.dt], [sys.umin(1), sys.umin(1)], ':k', 'DisplayName', 'u min');
xlabel('t [s]');
ylabel('u1');
title('Lecture 10 MPC Control');
legend('Location', 'best');

try
    exportgraphics(fig, path, 'Resolution', 160);
catch
    saveas(fig, path);
end
close(fig);
end

function write_summary_json(path, result, cases)
payload = struct();
payload.format_version = result.format_version;
payload.reference_name = result.reference_name;
payload.model = result.model;
payload.dt = result.dt;
payload.horizon_steps = result.horizon_steps;
payload.execution_steps = result.execution_steps;
payload.created_at = datestr(now, 31);
payload.version = version;

for i = 1:size(cases, 1)
    case_name = cases{i, 1};
    payload.cases.(case_name).final_state_norm = norm(result.cases.(case_name).final_state);
    payload.cases.(case_name).artifacts = result.cases.(case_name).artifacts;
end

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
