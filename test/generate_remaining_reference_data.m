function summary = generate_remaining_reference_data(save_root)
%GENERATE_REMAINING_REFERENCE_DATA Export reference artifacts for migrated examples.
if nargin < 1 || isempty(save_root)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    save_root = fullfile(repo_root, 'data');
end
if ~exist(save_root, 'dir'), mkdir(save_root); end

summary = struct();
summary.lecture2 = lecture2_reference(fullfile(save_root, 'lecture2'));
summary.lecture3 = lecture3_reference(fullfile(save_root, 'lecture3'));
summary.lecture12 = lecture12_reference(fullfile(save_root, 'lecture12'));
summary.lecture13 = lecture13_reference(fullfile(save_root, 'lecture13'));
summary.lecture14 = lecture14_reference(fullfile(save_root, 'lecture14'));
summary.lecture15 = lecture15_reference(fullfile(save_root, 'lecture15'));
summary.lecture16 = lecture16_reference(fullfile(save_root, 'lecture16'));
summary.lecture17 = lecture17_reference(fullfile(save_root, 'lecture17'));
summary.lecture18 = lecture18_reference(fullfile(save_root, 'lecture18'));
summary.lecture19 = lecture19_reference(fullfile(save_root, 'lecture19'));
summary.lecture21 = lecture21_reference(fullfile(save_root, 'lecture21'));
summary.lecture22 = lecture22_reference(fullfile(save_root, 'lecture22'));
summary.lecture23 = lecture23_reference(fullfile(save_root, 'lecture23'));
generate_reference_figures(save_root);
write_json(fullfile(save_root, 'reference_summary.json'), summary);
end

function summary = lecture2_reference(save_dir)
ensure_dir(save_dir);

dt = 0.1; tmax = 5; steps = round(tmax / dt); x0 = 0;
ctrl = @(x,t) 5 - x(1);
[t_ct, x_ct, u_ct] = rollout_ct(@(x,u) u, @(x,t) ctrl(x,t), x0, tmax, dt);
[t_zoh, x_zoh, u_zoh] = rollout_dt(@(x,u) u, @(x,t,k) ctrl(x,t), x0, steps, dt, 'zoh');
[t_euler, x_euler, u_euler] = rollout_dt(@(x,u) u, @(x,t,k) ctrl(x,t), x0, steps, dt, 'euler');
[t_rk4, x_rk4, u_rk4] = rollout_dt(@(x,u) u, @(x,t,k) ctrl(x,t), x0, steps, dt, 'rk4');
single_integrator_final = x_zoh(end,:);
write_csv(fullfile(save_dir, 'single_integrator_states.csv'), ...
    ["t","ct","zoh","euler","rk4"], [t_zoh(:), x_ct(:), x_zoh(:), x_euler(:), x_rk4(:)]);
write_csv(fullfile(save_dir, 'single_integrator_controls.csv'), ...
    ["t","ct","zoh","euler","rk4"], [t_zoh(1:end-1), u_ct(:), u_zoh(:), u_euler(:), u_rk4(:)]);

dt = 1; tmax = 5; steps = round(tmax / dt); x0 = [0; 0];
ctrl = @(x,t) 5 - x(1) - 2 * x(2);
[t_ct, x_ct, u_ct] = rollout_ct(@double_integrator_f, @(x,t) ctrl(x,t), x0, tmax, dt);
[t_zoh, x_zoh, u_zoh] = rollout_dt(@double_integrator_f, @(x,t,k) ctrl(x,t), x0, steps, dt, 'zoh');
[t_euler, x_euler, u_euler] = rollout_dt(@double_integrator_f, @(x,t,k) ctrl(x,t), x0, steps, dt, 'euler');
[t_rk4, x_rk4, u_rk4] = rollout_dt(@double_integrator_f, @(x,t,k) ctrl(x,t), x0, steps, dt, 'rk4');
write_csv(fullfile(save_dir, 'double_integrator_states.csv'), ...
    ["t","ct_x","ct_v","zoh_x","zoh_v","euler_x","euler_v","rk4_x","rk4_v"], ...
    [t_zoh(:), x_ct, x_zoh, x_euler, x_rk4]);
write_csv(fullfile(save_dir, 'double_integrator_controls.csv'), ...
    ["t","ct","zoh","euler","rk4"], [t_zoh(1:end-1), u_ct(:), u_zoh(:), u_euler(:), u_rk4(:)]);

goal = [5; 5]; x0 = [0; 0; 0]; dt = 0.1; steps = 50;
controller = @(x,t,k) [dot(goal - x(1:2), [cos(x(3)); sin(x(3))]); atan2(goal(2)-x(2), goal(1)-x(1)) - x(3)];
[t_uni, x_uni, u_uni] = rollout_dt(@unicycle3_f, controller, x0, steps, dt, 'euler');
write_csv(fullfile(save_dir, 'unicycle3_states.csv'), ["t","x","y","theta"], [t_uni(:), x_uni]);
write_csv(fullfile(save_dir, 'unicycle3_controls.csv'), ["t","v","omega"], [t_uni(1:end-1), u_uni]);

goal = [5; 5]; x0 = zeros(6,1); dt = 0.01; steps = 500;
controller = @(x,t,k) bicycle_controller(x, goal);
[t_bike, x_bike, u_bike] = rollout_dt(@bicycle_f, controller, x0, steps, dt, 'rk4', @(x,t) norm(x(1:2) - goal) <= 0.1);
write_csv(fullfile(save_dir, 'bicycle_dynamic_states.csv'), ["t","X","Y","vx","vy","r","psi"], [t_bike(:), x_bike]);
write_csv(fullfile(save_dir, 'bicycle_dynamic_controls.csv'), ["t","Fx","delta"], [t_bike(1:end-1), u_bike]);

summary.single_integrator_final = single_integrator_final;
summary.unicycle_goal_error = norm(x_uni(end,1:2)' - [5; 5]);
summary.bicycle_goal_error = norm(x_bike(end,1:2)' - [5; 5]);
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture3_reference(save_dir)
ensure_dir(save_dir);
rng(714);
dt = 0.1; steps = 50; q = zeros(7,1);
states = zeros(steps + 1, 7); controls = zeros(steps, 7); states(1,:) = q';
for k = 1:steps
    u = randn(7,1);
    q = q + dt * u;
    controls(k,:) = u';
    states(k+1,:) = q';
end
t = (0:steps)' * dt;
write_csv(fullfile(save_dir, 'random_joint_states.csv'), ["t","q1","q2","q3","q4","q5","q6","q7"], [t, states]);
write_csv(fullfile(save_dir, 'random_joint_controls.csv'), ["t","dq1","dq2","dq3","dq4","dq5","dq6","dq7"], [t(1:end-1), controls]);

q0 = [0.2; -0.5; 0.4]; q_goal = [0.6; -0.4; 0.8]; xy_goal = [1.6; 0.4]; lengths = [0.8; 0.6; 0.4];
q_joint = q0; q_cart = q0;
joint_states = zeros(steps+1,3); cart_states = zeros(steps+1,3);
joint_xy = zeros(steps+1,2); cart_xy = zeros(steps+1,2);
joint_controls = zeros(steps,3); cart_controls = zeros(steps,3);
joint_states(1,:) = q_joint'; cart_states(1,:) = q_cart';
joint_xy(1,:) = planar_fk(q_joint, lengths)'; cart_xy(1,:) = planar_fk(q_cart, lengths)';
for k = 1:steps
    u_joint = 0.7 * (q_goal - q_joint);
    q_joint = q_joint + dt * u_joint;
    joint_controls(k,:) = u_joint';
    joint_states(k+1,:) = q_joint';
    joint_xy(k+1,:) = planar_fk(q_joint, lengths)';

    qdot = 0.5 * pinv(planar_jacobian(q_cart, lengths)) * (xy_goal - planar_fk(q_cart, lengths));
    q_cart = q_cart + dt * qdot;
    cart_controls(k,:) = qdot';
    cart_states(k+1,:) = q_cart';
    cart_xy(k+1,:) = planar_fk(q_cart, lengths)';
end
write_csv(fullfile(save_dir, 'joint_servo_states.csv'), ["t","q1","q2","q3","ee_x","ee_y"], [t, joint_states, joint_xy]);
write_csv(fullfile(save_dir, 'joint_servo_controls.csv'), ["t","dq1","dq2","dq3"], [t(1:end-1), joint_controls]);
write_csv(fullfile(save_dir, 'cartesian_servo_states.csv'), ["t","q1","q2","q3","ee_x","ee_y"], [t, cart_states, cart_xy]);
write_csv(fullfile(save_dir, 'cartesian_servo_controls.csv'), ["t","dq1","dq2","dq3"], [t(1:end-1), cart_controls]);
summary.random_final_norm = norm(states(end,:));
summary.joint_final_error = norm(q_joint - q_goal);
summary.cartesian_final_error = norm(planar_fk(q_cart, lengths) - xy_goal);
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture12_reference(save_dir)
ensure_dir(save_dir);
result = run_frequency_ilc();
write_ilc_common(save_dir, result, true);
summary.initial_error_norm = result.error_norm(1);
summary.final_error_norm = result.error_norm(end);
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture13_reference(save_dir)
ensure_dir(save_dir);
result = run_time_ilc();
write_ilc_common(save_dir, result, false);
write_csv(fullfile(save_dir, 'learning_gain.csv'), ["row", compose("col%d", 0:size(result.L,2)-1)], [(0:size(result.L,1)-1)', result.L]);
summary.initial_error_norm = result.error_norm(1);
summary.final_error_norm = result.error_norm(end);
summary.gain_shape = size(result.L);
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture14_reference(save_dir)
ensure_dir(save_dir);
rng(714); n = 100; x = randn(2,n); AY = [1 0]; AZ = [1 1];
y = AY * x; z = AZ * x;
x_hat_y = AY' * y / (AY * AY');
z_hat_y = AZ * AY' * y / (AY * AY');
z_var_y = AZ * AZ' - AZ * AY' * AY * AZ' / (AY * AY');
x_hat_yz = x_hat_y + (AZ' - AY' * AY * AZ' / (AY * AY')) * (z - z_hat_y) / z_var_y;
write_csv(fullfile(save_dir, 'least_squares_samples.csv'), ...
    ["sample","x1","x2","y","z","xhat_y1","xhat_y2","xhat_yz1","xhat_yz2"], ...
    [(0:n-1)', x', y', z', x_hat_y', x_hat_yz']);
summary.max_error_y = max(abs(x - x_hat_y), [], 2);
summary.max_error_yz = max(abs(x - x_hat_yz), [], 2);
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture15_reference(save_dir)
ensure_dir(save_dir);
rng(714);
A = 0.5; B = 1; Bw = 1; C = 1; W = 0.01; V = 0.01; N = 101; dt = 1; x0 = 0; X0 = 1;
x0hat = x0 + randn * sqrt(X0);
x = zeros(N+1,1); u = zeros(N,1); y = zeros(N,1); x(1) = x0;
for k = 1:N
    t = (k-1) * dt;
    u(k) = sin(t / (N / 2 / pi));
    y(k) = C * x(k) + sqrt(V) * randn;
    x(k+1) = A * x(k) + B * u(k) + Bw * sqrt(W) * randn;
end
[xhat, z] = run_kf(A, B, C, Bw, W, V, x0hat, X0, u, y, false);
[xhat_ss, z_ss] = run_kf(A, B, C, Bw, W, V, x0hat, X0, u, y, true);
t = (0:N)' * dt;
write_csv(fullfile(save_dir, 'states.csv'), ["t","x"], [t, x]);
write_csv(fullfile(save_dir, 'controls.csv'), ["t","u"], [t(1:end-1), u]);
write_csv(fullfile(save_dir, 'estimates.csv'), ["t","x","y","kf","kfss","kf_var","kfss_var"], [t(1:N), x(1:N), y, xhat, xhat_ss, z, z_ss]);
summary.measurement_error = norm(y - x(1:N));
summary.kf_error = norm(xhat - x(1:N));
summary.kfss_error = norm(xhat_ss - x(1:N));
summary.seed = 714;
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture16_reference(save_dir)
ensure_dir(save_dir);
rng(714);
x0 = [4; 1]; X0 = 0.01 * [8 2; 2 3]; x0hat = x0 + sqrtm(X0) * randn(2,1);
W = 0.001; V = 0.001 * ones(2); Bw = ones(2,1); N = 10; dt = 1;
f_nom = @(x,u) x + u;
h_nom = @(x) [x(1)^4 * x(2); x(1) + x(2)^5];
Cfun = @(x) [4*x(1)^3*x(2), x(1)^4; 1, 5*x(2)^4];
x = zeros(N+1,2); y = zeros(N+1,2); u = zeros(N,2); x(1,:) = x0';
for k = 1:N
    y(k,:) = (h_nom(x(k,:)') + sqrt(0.001) * randn())';
    x(k+1,:) = (f_nom(x(k,:)', u(k,:)') + Bw * sqrt(W) * randn())';
end
y(N+1,:) = (h_nom(x(N+1,:)') + sqrt(0.001) * randn())';
[ekf_x, ekf_z] = run_ekf(f_nom, h_nom, Cfun, Bw, W, V, x0hat, X0, u, y);
[ukf_x, ukf_z] = run_ukf(f_nom, h_nom, Bw, W, V, x0hat, X0, u, y);
t = (0:N)' * dt;
write_csv(fullfile(save_dir, 'states.csv'), ["t","x1","x2"], [t, x]);
write_csv(fullfile(save_dir, 'measurements.csv'), ["t","y1","y2"], [t, y]);
write_csv(fullfile(save_dir, 'ekf_estimates.csv'), ["t","xhat1","xhat2"], [t(1:N), ekf_x]);
write_csv(fullfile(save_dir, 'ukf_estimates.csv'), ["t","xhat1","xhat2"], [t(1:N), ukf_x]);
write_csv(fullfile(save_dir, 'ekf_covariance.csv'), ["t","z11","z12","z21","z22"], [t(1:N), reshape_cov(ekf_z)]);
write_csv(fullfile(save_dir, 'ukf_covariance.csv'), ["t","z11","z12","z21","z22"], [t(1:N), reshape_cov(ukf_z)]);
summary.ekf_error = norm(ekf_x - x(1:N,:));
summary.ukf_error = norm(ukf_x - x(1:N,:));
summary.seed = 714;
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture17_reference(save_dir)
ensure_dir(save_dir);
N = 100; dt = 1; x = zeros(N+1,1); p_true = zeros(N,1); x(1) = 1;
for k = 1:N
    t = (k-1) * dt;
    p_true(k) = 1 + 0.05 * sin(t / 10);
    x(k+1) = p_true(k) * sin(x(k));
end
rls_p = zeros(N+1,1); sgd_p = zeros(N+1,1); rls_H = zeros(N+1,1);
rls_p(1) = 1; sgd_p(1) = 1; rls_H(1) = 0.5; lambda = 0.5; H0 = 0.5;
for k = 1:N
    g = sin(x(k));
    rls_p(k+1) = rls_p(k) + g * (x(k+1) - rls_p(k) * sin(x(k))) / (lambda * rls_H(k) + g * g);
    rls_H(k+1) = lambda * rls_H(k) + g * g;
    sgd_p(k+1) = sgd_p(k) + g * (x(k+1) - sgd_p(k) * sin(x(k))) / H0;
end
t = (0:N)' * dt;
write_csv(fullfile(save_dir, 'states.csv'), ["t","x"], [t, x]);
write_csv(fullfile(save_dir, 'parameter_estimates.csv'), ["t","p_true","rls","sgd","rls_learning_rate","sgd_learning_rate"], ...
    [t, [p_true; NaN], rls_p, sgd_p, 1 ./ rls_H, ones(N+1,1) / H0]);
summary.rls_final_error = abs(rls_p(end) - p_true(end));
summary.sgd_final_error = abs(sgd_p(end) - p_true(end));
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture18_reference(save_dir)
ensure_dir(save_dir);
rng(714);
dim = 4; dt = 0.2; N = 50; [A, B] = double_integrator_AB(dim, dt, 'zoh');
C = [eye(2), zeros(2)]; W = 0.1 * eye(2); V = 0.01 * eye(2); Bw = B;
x0 = [10; 5; 5; 0]; X0 = 0.001 * eye(dim); x0hat = x0 + sqrtm(X0) * randn(dim,1);
Q = [eye(2), zeros(2); zeros(2), zeros(2)]; R = eye(2);
[K, ~] = dlqr(A, B, Q, R);
x = zeros(N+1,dim); xhat = zeros(N+1,dim); y = zeros(N+1,2); u = zeros(N,2);
x(1,:) = x0'; xhat(1,:) = x0hat'; y(1,:) = (C * x0 + sqrtm(V) * randn(2,1))';
est = kfss_init(A, B, C, Bw, W, V, x0hat, X0);
for k = 1:N
    u(k,:) = (-K * xhat(k,:)')';
    x(k+1,:) = (A * x(k,:)' + B * u(k,:)' + Bw * sqrtm(W) * randn(2,1))';
    y(k+1,:) = (C * x(k+1,:)' + sqrtm(V) * randn(2,1))';
    [est, xhat(k+1,:)] = kfss_step(est, u(k,:)', y(k+1,:)');
end
t = (0:N)' * dt;
write_csv(fullfile(save_dir, 'states.csv'), ["t","x1","x2","v1","v2"], [t, x]);
write_csv(fullfile(save_dir, 'controls.csv'), ["t","u1","u2"], [t(1:end-1), u]);
write_csv(fullfile(save_dir, 'measurements.csv'), ["t","y1","y2"], [t, y]);
write_csv(fullfile(save_dir, 'estimates.csv'), ["t","xhat1","xhat2","vhat1","vhat2"], [t, xhat]);
summary.state_estimate_error = norm(xhat - x);
summary.final_state_norm = norm(x(end,:));
summary.seed = 714;
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture19_reference(save_dir)
ensure_dir(save_dir);
rng(714);
n = 2; m = 2; A = rand(n); B = rand(n,m); Astar = zeros(n); Ahat = eye(n); Bhat = eye(n,m);
x0 = [1; -1]; F = 10 * eye(n+m) + ones(n+m); N = 29; dt = 1;
[t, x, u, values, posterior] = run_mrac(A, B, Astar, Ahat, Bhat, F, x0, N, dt);
write_csv(fullfile(save_dir, 'states.csv'), ["t","x1","x2"], [t(:), x]);
write_csv(fullfile(save_dir, 'controls.csv'), ["t","u1","u2"], [t(1:end-1)', u]);
write_csv(fullfile(save_dir, 'analysis.csv'), ["t","value","posterior_error_1","posterior_error_2"], [(0:numel(values)-1)', values(:), posterior]);
summary.final_state_norm = norm(x(end,:));
summary.final_value = values(end);
summary.seed = 714;
summary.A = A; summary.B = B;
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture21_reference(save_dir)
ensure_dir(save_dir);
rng(714);
A = 1; B = 0.5; Q = 1; R = 1; tolerance = 0.01; kmax = 100; dt = 1; x0 = 1; n_ep = 100; alpha = 1;
[K, P] = dlqr(A, B, Q, R);
W_gt = [Q + A' * P * A, A' * P * B; B' * P * A, B' * P * B + R];
[~, x_lqr] = scalar_rollout(A, B, @(x,k) -K * x, x0, kmax, dt, tolerance);
[rms_mc, x_mc] = run_q_learning_family('mc', A, B, Q, R, [4 1; 1 4], W_gt, x0, kmax, dt, tolerance, alpha, 0.1, n_ep);
[rms_sarsa, x_sarsa] = run_q_learning_family('sarsa', A, B, Q, R, [4 1; 1 4], W_gt, x0, kmax, dt, tolerance, alpha, 0.1, n_ep);
[rms_q, x_q] = run_q_learning_family('qlearning', A, B, Q, R, [4 1; 1 4], W_gt, x0, kmax, dt, tolerance, alpha, 0.1, n_ep);
write_csv(fullfile(save_dir, 'rms_errors.csv'), ["episode","sarsa","qlearning","mc"], [(0:n_ep-1)', rms_sarsa(:), rms_q(:), rms_mc(:)]);
write_csv(fullfile(save_dir, 'lqr_states.csv'), ["k","x"], [(0:numel(x_lqr)-1)', x_lqr(:)]);
write_csv(fullfile(save_dir, 'mc_final_states.csv'), ["k","x"], [(0:numel(x_mc)-1)', x_mc(:)]);
write_csv(fullfile(save_dir, 'sarsa_final_states.csv'), ["k","x"], [(0:numel(x_sarsa)-1)', x_sarsa(:)]);
write_csv(fullfile(save_dir, 'qlearning_final_states.csv'), ["k","x"], [(0:numel(x_q)-1)', x_q(:)]);
summary.final_rms_mc = rms_mc(end);
summary.final_rms_sarsa = rms_sarsa(end);
summary.final_rms_qlearning = rms_q(end);
summary.seed = 714;
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture22_reference(save_dir)
ensure_dir(save_dir);
rng(714);
A = 0.5; B = 1; Q = 1; R = 1; x0 = 0.1; kmax = 50; dt = 1; n_ep = 20; seeds = randi(1e6, n_ep, 1);
[K, ~] = dlqr(A, B, Q, R);
[~, x_lqr] = scalar_rollout(A, B, @(x,k) -K*x, x0, kmax, dt, -1);
[pg_metrics] = run_pg_set(A, B, Q, R, K, x0, kmax, dt, n_ep, seeds, true);
[ilc_error, ilc_state] = run_scalar_ilc(A, B, x0, kmax, dt, n_ep, seeds, x_lqr);
write_csv(fullfile(save_dir, 'policy_gradient_metrics.csv'), ...
    ["episode","reinforce","reinforce_bl","actor_critic","reinforce_std","reinforce_bl_std","actor_critic_std"], ...
    [(0:n_ep-1)', pg_metrics]);
write_csv(fullfile(save_dir, 'ilc_error_norm.csv'), ["iteration","error_norm"], [(0:n_ep-1)', ilc_error(:)]);
write_csv(fullfile(save_dir, 'lqr_states.csv'), ["k","x"], [(0:numel(x_lqr)-1)', x_lqr(:)]);
write_csv(fullfile(save_dir, 'ilc_final_states.csv'), ["k","x"], [(0:numel(ilc_state)-1)', ilc_state(:)]);
summary.final_rms_reinforce = pg_metrics(end,1);
summary.final_rms_reinforce_bl = pg_metrics(end,2);
summary.final_rms_actor_critic = pg_metrics(end,3);
summary.final_ilc_error = ilc_error(end);
summary.seed = 714;
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function summary = lecture23_reference(save_dir)
ensure_dir(save_dir);
rng(714);
A = 1; B = 0.5; Q = 1; R = 1; x0 = 1; kmax = 9; dt = 1; n_ep = 200; tolerance = 0.01;
[K, ~] = dlqr(A, B, Q, R);
[~, x_lqr] = scalar_rollout(A, B, @(x,k) -K*x, x0, kmax, dt, tolerance);
pg_metrics = run_pg_set(A, B, Q, R, K, x0, kmax, dt, n_ep, [], false);
write_csv(fullfile(save_dir, 'policy_gradient_metrics.csv'), ...
    ["episode","reinforce","reinforce_bl","actor_critic","reinforce_std","reinforce_bl_std","actor_critic_std"], ...
    [(0:n_ep-1)', pg_metrics]);
write_csv(fullfile(save_dir, 'lqr_states.csv'), ["k","x"], [(0:numel(x_lqr)-1)', x_lqr(:)]);
summary.final_rms_reinforce = pg_metrics(end,1);
summary.final_rms_reinforce_bl = pg_metrics(end,2);
summary.final_rms_actor_critic = pg_metrics(end,3);
summary.seed = 714;
write_json(fullfile(save_dir, 'summary.json'), summary);
end

function write_ilc_common(save_dir, result, write_feedforward)
idx = size(result.states, 3);
write_csv(fullfile(save_dir, 'states.csv'), ["t","x","v","reference"], [result.t(:), result.states(:,:,idx), result.reference(:)]);
write_csv(fullfile(save_dir, 'controls.csv'), ["t","u"], [result.control_t(:), result.controls(:,:,idx)]);
write_csv(fullfile(save_dir, 'error_norm.csv'), ["iteration","error_norm"], [(0:numel(result.error_norm)-1)', result.error_norm(:)]);
if write_feedforward
    write_csv(fullfile(save_dir, 'feedforward.csv'), ["iteration", compose("u%d", 0:size(result.feedforward,2)-1)], [(0:size(result.feedforward,1)-1)', result.feedforward]);
end
end

function result = run_frequency_ilc()
dt = 1; horizon = 100; niter = 10; kp = 0.1; kd = 0.5;
t = (0:horizon)' * dt; ref = sin(0.2 * t);
states = zeros(horizon+1,2,niter+1); controls = zeros(horizon,1,niter+1);
[~, states(:,:,1), controls(:,:,1)] = ilc_rollout(zeros(1,horizon+1), true);
L_b = poly_add([1 -2 1], conv(dt, -[-kd/dt, kd/dt-kp])); L_a = dt;
errors = zeros(niter,horizon+1); ff = zeros(niter,horizon+1); e = zeros(niter,1);
for i = 1:niter
    errors(i,:) = ref' - states(:,1,i)';
    prev = zeros(1,horizon+1); if i > 1, prev = ff(i-1,:); end
    learned = ilc_filter(L_b, L_a, errors(i,:));
    ff(i,:) = filter(1, 1, prev + learned);
    [~, states(:,:,i+1), controls(:,:,i+1)] = ilc_rollout(ff(i,:), true);
    e(i) = norm(errors(i,:));
end
result = struct('t', t, 'control_t', t(1:end-1), 'reference', ref, 'states', states, 'controls', controls, ...
    'errors', errors, 'feedforward', ff, 'error_norm', e, 'L_b', L_b, 'L_a', L_a);

    function [tout, xout, uout] = ilc_rollout(feedforward, full_length)
        x = zeros(horizon+1,2); u = zeros(horizon,1); tout = t;
        for k = 1:horizon
            tk = (k-1) * dt;
            ufb = kp * (sin(0.2*tk) - x(k,1)) + kd * ((sin(0.2*(tk+dt)) - sin(0.2*tk))/dt - x(k,2));
            u(k) = feedforward(k) + ufb;
            x(k+1,:) = x(k,:) + dt * [x(k,2), u(k)];
        end
        xout = x; uout = u;
    end
end

function result = run_time_ilc()
dt = 1; horizon = 100; niter = 10; kp = 0.1; kd = 0.5;
t = (0:horizon)' * dt; ref = sin(0.2 * t);
[A, B] = double_integrator_AB(2, dt, 'euler'); K = [kp kd]; C = [1 0];
policy_L = time_ilc_gain(A, B, K, C, horizon);
states = zeros(horizon+1,2,niter+1); controls = zeros(horizon,1,niter+1);
[~, states(:,:,1), controls(:,:,1)] = ilc_rollout(zeros(1,horizon));
errors = zeros(niter,horizon+1); ff = zeros(niter,horizon); e = zeros(niter,1);
for i = 1:niter
    errors(i,:) = ref' - states(:,1,i)';
    prev = zeros(1,horizon); if i > 1, prev = ff(i-1,:); end
    ff(i,:) = prev + errors(i,:) * policy_L';
    [~, states(:,:,i+1), controls(:,:,i+1)] = ilc_rollout(ff(i,:));
    e(i) = norm(errors(i,:));
end
result = struct('t', t, 'control_t', t(1:end-1), 'reference', ref, 'states', states, 'controls', controls, ...
    'errors', errors, 'feedforward', ff, 'error_norm', e, 'L', policy_L);

    function [tout, xout, uout] = ilc_rollout(feedforward)
        x = zeros(horizon+1,2); u = zeros(horizon,1); tout = t;
        for k = 1:horizon
            tk = (k-1) * dt;
            ufb = kp * (sin(0.2*tk) - x(k,1)) + kd * ((sin(0.2*(tk+dt)) - sin(0.2*tk))/dt - x(k,2));
            u(k) = feedforward(k) + ufb;
            x(k+1,:) = x(k,:) + dt * [x(k,2), u(k)];
        end
        xout = x; uout = u;
    end
end

function [t, x, u] = rollout_dt(f, controller, x0, steps, dt, integrator, stop)
if nargin < 7, stop = []; end
x = zeros(steps+1, numel(x0)); u = []; t = zeros(steps+1,1); x(1,:) = x0(:)';
last = steps;
for k = 1:steps
    tk = (k-1) * dt;
    if ~isempty(stop) && stop(x(k,:)', tk), last = k-1; break; end
    uk = controller(x(k,:)', tk, k-1);
    if isempty(u), u = zeros(steps, numel(uk)); end
    u(k,:) = uk(:)';
    x(k+1,:) = step_state(f, x(k,:)', uk(:), dt, integrator)';
    t(k+1) = k * dt;
end
x = x(1:last+1,:); t = t(1:last+1); u = u(1:last,:);
end

function [t, x, u] = rollout_ct(f, controller, x0, tmax, sample_dt)
t = (0:sample_dt:tmax)';
ode = @(tt, xx) f(xx, controller(xx, tt));
[~, x] = ode45(ode, t, x0(:));
u = zeros(numel(t)-1, numel(controller(x(1,:)', t(1))));
for k = 1:numel(t)-1
    u(k,:) = controller(x(k,:)', t(k))';
end
end

function xnext = step_state(f, x, u, dt, integrator)
switch lower(integrator)
    case {'euler','zoh'}
        xnext = x + dt * f(x, u);
    case 'rk4'
        k1 = f(x, u); k2 = f(x + 0.5*dt*k1, u); k3 = f(x + 0.5*dt*k2, u); k4 = f(x + dt*k3, u);
        xnext = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
    otherwise
        error('Unknown integrator %s', integrator);
end
end

function dx = double_integrator_f(x, u)
n = numel(x) / 2; dx = [x(n+1:end); u(:)];
end

function dx = unicycle3_f(x, u)
dx = [u(1)*cos(x(3)); u(1)*sin(x(3)); u(2)];
end

function u = bicycle_controller(x, goal)
p = bicycle_params();
X = x(1); Y = x(2); vx = x(3); psi = x(6);
psi_des = atan2(goal(2)-Y, goal(1)-X);
e_psi = wrap_to_pi(psi_des - psi);
Rwb = [cos(psi), sin(psi); -sin(psi), cos(psi)];
e_vec_b = Rwb * (goal - [X; Y]);
v_long = max(abs(vx), p.eps_vx) * sign(vx + (vx == 0));
delta = e_psi + atan2(p.k_cte * e_vec_b(2), v_long);
delta = min(max(delta, -p.delta_max), p.delta_max);
dist = norm(goal - [X; Y]);
v_ref = min(p.v_ref_max, p.k_speed * dist);
fx = p.m * p.k_v * (v_ref - vx);
fx = min(max(fx, -p.Fmax), p.Fmax);
u = [fx; delta];
end

function dx = bicycle_f(x, u)
p = bicycle_params();
vx = x(3); vy = x(4); r = x(5); psi = x(6); fx = u(1); delta = u(2);
vx_eff = vx; if abs(vx_eff) < p.eps_vx, vx_eff = sign(vx_eff + (vx_eff == 0)) * p.eps_vx; end
alpha_f = atan2(vy + p.lf*r, vx_eff) - delta;
alpha_r = atan2(vy - p.lr*r, vx_eff);
fyf = -p.Cf * alpha_f; fyr = -p.Cr * alpha_r;
fxf = p.beta * fx; fxr = (1 - p.beta) * fx;
dvx = (fxf*cos(delta) - fyf*sin(delta) + fxr) / p.m + r * vy;
dvy = (fxf*sin(delta) + fyf*cos(delta) + fyr) / p.m - r * vx;
dr = (p.lf*(fxf*sin(delta) + fyf*cos(delta)) - p.lr*fyr) / p.Iz;
dx = [vx*cos(psi) - vy*sin(psi); vx*sin(psi) + vy*cos(psi); dvx; dvy; dr; r];
end

function p = bicycle_params()
p = struct('m',1500,'Iz',2250,'lf',1.2,'lr',1.6,'Cf',8e4,'Cr',9e4,'beta',0.5, ...
    'eps_vx',0.1,'k_cte',2.0,'k_v',1.2,'k_speed',0.8,'v_ref_max',8,'Fmax',6000,'delta_max',0.6);
end

function xy = planar_fk(q, lengths)
angles = cumsum(q(:));
xy = [sum(lengths(:) .* cos(angles)); sum(lengths(:) .* sin(angles))];
end

function J = planar_jacobian(q, lengths)
angles = cumsum(q(:)); J = zeros(2,3);
for j = 1:3
    J(1,j) = -sum(lengths(j:end) .* sin(angles(j:end)));
    J(2,j) =  sum(lengths(j:end) .* cos(angles(j:end)));
end
end

function angle = wrap_to_pi(angle)
angle = mod(angle + pi, 2*pi) - pi;
end

function [A, B] = double_integrator_AB(dim, dt, mode)
n = dim / 2;
A = [eye(n), dt * eye(n); zeros(n), eye(n)];
if strcmpi(mode, 'euler')
    B = [zeros(n); dt * eye(n)];
else
    B = [0.5 * dt^2 * eye(n); dt * eye(n)];
end
end

function y = poly_add(a, b)
n = max(numel(a), numel(b)); y = zeros(1,n);
y(end-numel(a)+1:end) = y(end-numel(a)+1:end) + a;
y(end-numel(b)+1:end) = y(end-numel(b)+1:end) + b;
end

function y = ilc_filter(b, a, error)
nshift = numel(b) - numel(a);
y = filter(b, a, error);
if nshift > 0, y(1:end-nshift) = y(nshift+1:end); end
end

function L = time_ilc_gain(A, B, K, C, horizon)
Acl = A - B * K;
[~, bar_B] = lift_dynamics(Acl, B, horizon);
c_pinv = pinv(C);
bar_inv_C = zeros((horizon+1) * size(c_pinv,1), horizon+1);
for i = 1:horizon+1
    rows = (i-1)*size(c_pinv,1)+1:i*size(c_pinv,1);
    bar_inv_C(rows,i) = c_pinv;
end
L = pinv(bar_B) * bar_inv_C;
end

function [bar_A, bar_B] = lift_dynamics(A, B, horizon)
nx = size(A,1); nu = size(B,2);
bar_A = zeros((horizon+1)*nx, nx); bar_B = zeros((horizon+1)*nx, horizon*nu);
bar_A(1:nx,:) = eye(nx);
for k = 1:horizon
    prev = (k-1)*nx+1:k*nx; cur = k*nx+1:(k+1)*nx; col = (k-1)*nu+1:k*nu;
    bar_A(cur,:) = A * bar_A(prev,:);
    bar_B(cur,:) = A * bar_B(prev,:);
    bar_B(cur,col) = B;
end
end

function [xhat, z] = run_kf(A, B, C, Bw, W, V, x0hat, X0, u, y, steady)
N = numel(u); xhat = zeros(N,1); z = zeros(N,1); xhat(1) = x0hat; z(1) = X0;
if steady
    Ms = dare(A', C', Bw*W*Bw', V);
end
for k = 1:N-1
    if steady
        M = Ms;
    else
        M = A * z(k) * A' + Bw * W * Bw';
    end
    gain = M * C' / (V + C * M * C');
    xprior = A * xhat(k) + B * u(k);
    xhat(k+1) = xprior + gain * (y(k+1) - C * xprior);
    z(k+1) = M - M * C' / (V + C * M * C') * C * M;
end
end

function [ekf_x, ekf_z] = run_ekf(f_nom, h_nom, Cfun, Bw, W, V, x0hat, X0, u, y)
N = size(u,1); ekf_x = zeros(N,2); ekf_z = zeros(2,2,N); ekf_x(1,:) = x0hat'; ekf_z(:,:,1) = X0;
for k = 1:N-1
    xcur = ekf_x(k,:)'; C = Cfun(xcur); A = eye(2);
    M = A * ekf_z(:,:,k) * A' + Bw * W * Bw';
    gain = M * C' / (V + C * M * C');
    xprior = f_nom(xcur, u(k,:)');
    ekf_x(k+1,:) = (xprior + gain * (y(k+1,:)' - h_nom(xcur)))';
    ekf_z(:,:,k+1) = M - M * C' / (C * M * C' + V) * C * M;
end
end

function [ukf_x, ukf_z] = run_ukf(f_nom, h_nom, Bw, W, V, x0hat, X0, u, y)
N = size(u,1); ukf_x = zeros(N,2); ukf_z = zeros(2,2,N); ukf_x(1,:) = x0hat'; ukf_z(:,:,1) = X0;
for k = 1:N-1
    [x_points, weights] = sigma_points(ukf_x(k,:)', ukf_z(:,:,k));
    for j = 1:size(x_points,2), x_points(:,j) = f_nom(x_points(:,j), u(k,:)'); end
    xprior = weighted_mean(x_points, weights);
    M = weighted_cov(x_points, weights) + Bw * W * Bw';
    [measure_points, weights2] = sigma_points(xprior, M);
    y_points = zeros(2, size(measure_points,2));
    for j = 1:size(measure_points,2), y_points(:,j) = h_nom(measure_points(:,j)); end
    yprior = weighted_mean(y_points, weights2);
    yvar = weighted_cov(y_points, weights2);
    xyvar = cross_cov(measure_points, xprior, y_points, yprior, weights2);
    gain = xyvar / (V + yvar);
    ukf_x(k+1,:) = (xprior + gain * (y(k+1,:)' - yprior))';
    ukf_z(:,:,k+1) = M - xyvar / (V + yvar) * xyvar';
end
end

function [points, weights] = sigma_points(mean, covar)
kappa = 2; n = numel(mean); points = zeros(n, 1 + 2*n); weights = zeros(1, 1 + 2*n);
points(:,1) = mean; weights(1) = kappa / (n + kappa);
root = real(sqrtm((n + kappa) * covar));
for i = 1:n
    points(:,2*i) = mean + root(:,i); points(:,2*i+1) = mean - root(:,i);
    weights(2*i) = 1 / (2 * (n + kappa)); weights(2*i+1) = weights(2*i);
end
end

function mean = weighted_mean(points, weights)
mean = points * weights(:) / sum(weights);
end

function covar = weighted_cov(points, weights)
mean = weighted_mean(points, weights); covar = zeros(size(points,1));
for i = 1:size(points,2)
    d = points(:,i) - mean; covar = covar + weights(i) * (d * d');
end
covar = covar / sum(weights);
end

function covar = cross_cov(x_points, x_mean, y_points, y_mean, weights)
covar = zeros(size(x_points,1), size(y_points,1));
for i = 1:size(x_points,2)
    covar = covar + weights(i) * ((x_points(:,i) - x_mean) * (y_points(:,i) - y_mean)');
end
covar = covar / sum(weights);
end

function data = reshape_cov(covariances)
N = size(covariances,3); data = zeros(N,4);
for k = 1:N
    z = covariances(:,:,k); data(k,:) = z(:)';
end
end

function est = kfss_init(A, B, C, Bw, W, V, xhat, Z)
Ms = dare(A', C', Bw * W * Bw', V);
est = struct('A',A,'B',B,'C',C,'V',V,'Ms',Ms,'xhat',xhat,'Z',Z);
est.gain = Ms * C' / (V + C * Ms * C');
est.Zss = Ms - Ms * C' / (V + C * Ms * C') * C * Ms;
end

function [est, xhat] = kfss_step(est, u, y)
xprior = est.A * est.xhat + est.B * u;
est.xhat = xprior + est.gain * (y - est.C * xprior);
est.Z = est.Zss; xhat = est.xhat';
end

function [t, states, controls, values, posterior] = run_mrac(A, B, Astar, Ahat, Bhat, F, x0, N, dt)
n = size(A,1); x = x0; phi = []; ref = zeros(N+1,n); ref(1,:) = x0';
for k = 1:N, ref(k+1,:) = (Astar * ref(k,:)')'; end
states = zeros(N+1,n); controls = zeros(N,size(B,2)); states(1,:) = x';
A_hist = zeros(n,n,N+1); B_hist = zeros(n,size(B,2),N+1); F_hist = zeros(size(F,1),size(F,2),N+1); err_hist = zeros(N+1,n); phi_hist = zeros(N,n+size(B,2));
A_hist(:,:,1) = Ahat; B_hist(:,:,1) = Bhat; F_hist(:,:,1) = F; err_hist(1,:) = x';
for k = 1:N
    tk = (k-1) * dt;
    if tk > 0 && ~isempty(phi)
        error = x - ref(k,:)';
        F = inv(inv(F) + phi * phi');
        Ahat = Ahat + error * (phi' * F(:,1:n));
        Bhat = Bhat + error * (phi' * F(:,n+1:end));
        A_hist(:,:,k) = Ahat; B_hist(:,:,k) = Bhat; F_hist(:,:,k) = F; err_hist(k,:) = error'; phi_hist(k-1,:) = phi';
    else
        A_hist(:,:,k) = Ahat; B_hist(:,:,k) = Bhat; F_hist(:,:,k) = F;
    end
    K = pinv(Bhat) * (Astar - Ahat);
    u = K * x; controls(k,:) = u';
    phi = [x; u];
    x = A * x + B * u;
    states(k+1,:) = x';
end
A_hist(:,:,N+1) = Ahat; B_hist(:,:,N+1) = Bhat; F_hist(:,:,N+1) = F; err_hist(N+1,:) = (x - ref(N+1,:)')';
values = zeros(N+1,1); posterior = zeros(N+1,n);
for k = 1:N+1
    ABtilde = [A_hist(:,:,k) - A, B_hist(:,:,k) - B];
    if k > 1
        ph = phi_hist(k-1,:)';
        posterior(k,:) = ((1 - ph' * F_hist(:,:,k) * ph) * err_hist(k,:)')';
    else
        posterior(k,:) = err_hist(k,:);
    end
    values(k) = norm(posterior(k,:))^2 + trace(ABtilde / F_hist(:,:,k) * ABtilde');
end
t = 0:N;
end

function [t, x, u] = scalar_rollout(A, B, controller, x0, kmax, dt, tolerance)
x = x0; u = [];
for k = 1:kmax
    if tolerance >= 0 && norm(x(end)) < tolerance, break; end
    if (k-1) * dt > kmax * dt, break; end
    uk = controller(x(end), k-1);
    u(end+1,1) = uk;
    x(end+1,1) = A * x(end) + B * uk;
end
t = (0:numel(x)-1)' * dt;
end

function [rms_values, final_states] = run_q_learning_family(alg, A, B, Q, R, W, Wgt, x0, kmax, dt, tolerance, alpha, epsilon, n_ep)
rms_values = zeros(n_ep,1); final_states = x0;
for ep = 1:n_ep
    [W, states] = q_episode(alg, A, B, Q, R, W, x0, kmax, dt, tolerance, alpha, epsilon);
    final_states = states;
    rms_values(ep) = norm(W - Wgt);
end
end

function [W, states] = q_episode(alg, A, B, Q, R, W, x0, kmax, dt, tolerance, alpha, epsilon)
stage = @(x,u) (Q*x^2 + R*u^2)/2; qfun = @(x,u,Wm) [x;u]' * Wm * [x;u] / 2;
grad = @(x,u) ([x;u] * [x;u]') / 2; value = @(x,Wm) x^2 * (Wm(1,1) - Wm(1,2) / Wm(2,2) * Wm(2,1)) / 2;
greedy = @(x,Wm) q_greedy(x, Wm, epsilon);
states = x0; controls = []; costs = [];
switch lower(alg)
    case 'mc'
        x = x0;
        for k = 1:kmax
            if norm(x) < tolerance, break; end
            u = greedy(x, W); controls(end+1,1) = u; costs(end+1,1) = stage(x,u);
            x = A*x + B*u; states(end+1,1) = x;
        end
        dW = zeros(size(W));
        for k = 1:numel(costs)
            G = sum(costs(k:end));
            dW = dW + alpha * (G - qfun(states(k), controls(k), W)) * grad(states(k), controls(k));
        end
        W = W + dW;
    case 'sarsa'
        x = x0; u = greedy(x, W); k = 0;
        while norm(x) >= tolerance && k * dt <= kmax * dt
            xnew = A*x + B*u; l = stage(x,u); unew = greedy(xnew, W);
            W = W + alpha * (l + qfun(xnew, unew, W) - qfun(x, u, W)) * grad(x,u);
            states(end+1,1) = xnew; x = xnew; u = unew; k = k + 1;
        end
    case 'qlearning'
        x = x0; k = 0;
        while norm(x) >= tolerance && k * dt <= kmax * dt
            u = greedy(x, W); xnew = A*x + B*u; l = stage(x,u);
            W = W + alpha * (l + value(xnew, W) - qfun(x, u, W)) * grad(x,u);
            states(end+1,1) = xnew; x = xnew; k = k + 1;
        end
end
end

function u = q_greedy(x, W, epsilon)
u_star = -W(2,1) / W(2,2) * x;
if rand < epsilon
    u = u_star - rand * x;
else
    u = u_star;
end
end

function metrics = run_pg_set(A, B, Q, R, K, x0, kmax, dt, n_ep, seeds, noisy)
metrics = zeros(n_ep,6);
agents = [struct('mu',-0.1,'sigma',0.1,'W',3), struct('mu',-0.1,'sigma',0.1,'W',3), struct('mu',-0.1,'sigma',0.1,'W',3)];
for ep = 1:n_ep
    if ~isempty(seeds), rng(seeds(ep)); end
    agents(1) = pg_update('reinforce', agents(1), A, B, Q, R, x0, kmax, dt, 0.01 * (~noisy) + 1e-4 * noisy, 1, noisy);
    if ~isempty(seeds), rng(seeds(ep)); end
    agents(2) = pg_update('baseline', agents(2), A, B, Q, R, x0, kmax, dt, 0.01 * (~noisy) + 1e-4 * noisy, 1.0 * (~noisy) + 1e-2 * noisy, noisy);
    if ~isempty(seeds), rng(seeds(ep)); end
    agents(3) = pg_update('actorcritic', agents(3), A, B, Q, R, x0, kmax, dt, 0.01 * (~noisy) + 1e-4 * noisy, 1.0 * (~noisy) + 1e-2 * noisy, noisy);
    for i = 1:3
        metrics(ep,i) = abs(agents(i).mu + K);
        metrics(ep,i+3) = agents(i).sigma;
    end
end
end

function agent = pg_update(alg, agent, A, B, Q, R, x0, kmax, dt, alpha_p, alpha_v, noisy)
stage = @(x,u) (Q*x^2 + R*u^2)/2;
sample = @(x) agent.mu * x + randn * agent.sigma * x;
gradp = @(x,u) [safe_div((u-agent.mu*x)*x, (agent.sigma*x)^2); -x + safe_div((u-agent.mu*x)^2*x, (agent.sigma*x)^2)];
switch alg
    case {'reinforce','baseline'}
        x = x0; xs = x; us = []; costs = [];
        for k = 1:kmax
            if ~noisy && norm(x) < 0.01, break; end
            u = sample(x); us(end+1,1) = u; costs(end+1,1) = stage(x,u);
            d = 0; if noisy, d = randn; end
            x = A*x + B*u + d; xs(end+1,1) = x;
        end
        dtheta = [0; 0]; dW = 0;
        for k = 1:numel(costs)
            G = sum(costs(k:end)); base = 0;
            if strcmp(alg, 'baseline'), base = agent.W * xs(k)^2 / 2; end
            dtheta = dtheta - bounded_step(alpha_p * (G - base), gradp(xs(k), us(k)));
            if strcmp(alg, 'baseline'), dW = dW + bounded_step(alpha_v * (G - base), xs(k)^2 / 2); end
        end
        agent.mu = agent.mu + dtheta(1); agent.sigma = agent.sigma + dtheta(2); agent.W = agent.W + dW;
    case 'actorcritic'
        x = x0;
        for k = 1:kmax
            if ~noisy && norm(x) < 0.01, break; end
            u = sample(x); d = 0; if noisy, d = randn; end
            xnew = A*x + B*u + d; l = stage(x,u);
            td = l + agent.W * xnew^2 / 2 - agent.W * x^2 / 2;
            step = bounded_step(alpha_p * td, gradp(x,u));
            agent.mu = agent.mu - step(1); agent.sigma = agent.sigma - step(2);
            agent.W = agent.W + bounded_step(alpha_v * td, x^2 / 2);
            x = xnew;
        end
end
agent.mu = min(max(agent.mu, -5), 5); agent.sigma = min(max(agent.sigma, 1e-4), 5);
if ~isfinite(agent.mu), agent.mu = 0; end
if ~isfinite(agent.sigma), agent.sigma = 1e-4; end
end

function out = bounded_step(scale, direction)
out = scale .* direction;
out(~isfinite(out)) = 0;
out = min(max(out, -1e6), 1e6);
end

function y = safe_div(a, b)
if abs(b) < 1e-14 || ~isfinite(a) || ~isfinite(b)
    y = 0;
else
    y = a / b;
end
end

function [errors, final_state] = run_scalar_ilc(A, B, x0, horizon, dt, niter, seeds, reference)
L = time_ilc_gain(A, B, 0, 1, horizon);
ff = zeros(niter,horizon); errors = zeros(niter,1);
rng(seeds(1)); x = scalar_noisy_rollout(A, B, zeros(1,horizon), x0, horizon);
for i = 1:niter
    e = reference(:)' - x(:)';
    prev = zeros(1,horizon); if i > 1, prev = ff(i-1,:); end
    ff(i,:) = prev + e * L';
    rng(seeds(i)); x = scalar_noisy_rollout(A, B, ff(i,:), x0, horizon);
    errors(i) = norm(e);
end
final_state = x;
end

function x = scalar_noisy_rollout(A, B, controls, x0, horizon)
x = zeros(horizon+1,1); x(1) = x0;
for k = 1:horizon
    x(k+1) = A*x(k) + B*controls(k) + randn;
end
end

function generate_reference_figures(save_root)
plot_legacy_case_figures(save_root);
plot_lecture2_figures(fullfile(save_root, 'lecture2'));
plot_lecture3_figures(fullfile(save_root, 'lecture3'));
plot_ilc_figures(fullfile(save_root, 'lecture12'), 'frequency_ilc_trajectory.png', 'frequency_ilc_error.png');
plot_ilc_figures(fullfile(save_root, 'lecture13'), 'time_domain_ilc_trajectory.png', 'time_domain_ilc_summary.png');
plot_kalman_figure(fullfile(save_root, 'lecture15'));
plot_ekf_ukf_figure(fullfile(save_root, 'lecture16'));
plot_rls_sgd_figure(fullfile(save_root, 'lecture17'));
plot_separation_figure(fullfile(save_root, 'lecture18'));
plot_mrac_figure(fullfile(save_root, 'lecture19'));
plot_value_learning_figure(fullfile(save_root, 'lecture21'));
plot_pg_ilc_figure(fullfile(save_root, 'lecture22'));
plot_actor_critic_figure(fullfile(save_root, 'lecture23'));
remove_legacy_extra_figures(save_root);
end

function plot_legacy_case_figures(save_root)
src = fullfile(save_root, 'lecture10', 'reference_ilqr_state_plot.png');
dst = fullfile(save_root, 'lecture10', 'ilqr_state_plot.png');
if isfile(src), copyfile(src, dst); end

for case_name = ["lqr", "finite_horizon_lqr", "linear_mpc"]
    case_dir = fullfile(save_root, 'lecture9', char(case_name));
    plot_state_control_case(case_dir, fullfile(case_dir, 'state_plot.png'), "Lecture 9 " + case_name);
end
for case_name = ["linear_mpc", "input_constrained_mpc", "constrained_mpc"]
    case_dir = fullfile(save_root, 'lecture10_mpc', char(case_name));
    plot_state_control_case(case_dir, fullfile(case_dir, 'state_plot.png'), "Lecture 10 MPC " + case_name);
end
end

function plot_state_control_case(case_dir, output_path, title_text)
states_path = fullfile(case_dir, 'states.csv');
controls_path = fullfile(case_dir, 'controls.csv');
if ~isfile(states_path) || ~isfile(controls_path), return; end
states = readmatrix(states_path);
controls = readmatrix(controls_path);
fig = figure('Visible','off','Position',[100 100 760 460]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile([2 1]); hold on; grid on;
if size(states,2) >= 3
    plot(states(:,2), states(:,3), 'LineWidth', 1.5);
    xlabel('x1'); ylabel('x2');
else
    plot(states(:,1), states(:,2), 'LineWidth', 1.5);
    xlabel('t'); ylabel('x');
end
title('Trajectory');
nexttile; hold on; grid on;
plot(states(:,1), states(:,2:end), 'LineWidth', 1.2);
xlabel('t'); ylabel('state'); title('State');
nexttile; hold on; grid on;
plot(controls(:,1), controls(:,2:end), 'LineWidth', 1.2);
xlabel('t'); ylabel('control'); title('Control');
sgtitle(title_text, 'Interpreter','none');
save_png(fig, output_path);
end

function plot_lecture2_figures(save_dir)
states = readmatrix(fullfile(save_dir, 'single_integrator_states.csv'));
controls = readmatrix(fullfile(save_dir, 'single_integrator_controls.csv'));
fig = figure('Visible','off','Position',[100 100 700 500]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; plot(states(:,1), states(:,2:end), 'LineWidth',1.2); grid on; title('Single Integrator State'); legend('ct','zoh','euler','rk4','Location','best');
nexttile; plot(controls(:,1), controls(:,2:end), 'LineWidth',1.2); grid on; title('Single Integrator Control'); legend('ct','zoh','euler','rk4','Location','best');
save_png(fig, fullfile(save_dir, 'single_integrator_plot.png'));

states = readmatrix(fullfile(save_dir, 'double_integrator_states.csv'));
controls = readmatrix(fullfile(save_dir, 'double_integrator_controls.csv'));
fig = figure('Visible','off','Position',[100 100 700 620]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile; plot(states(:,1), states(:,[2 4 6 8]), 'LineWidth',1.2); grid on; title('Double Integrator Position'); legend('ct','zoh','euler','rk4','Location','best');
nexttile; plot(states(:,1), states(:,[3 5 7 9]), 'LineWidth',1.2); grid on; title('Double Integrator Velocity'); legend('ct','zoh','euler','rk4','Location','best');
nexttile; plot(controls(:,1), controls(:,2:end), 'LineWidth',1.2); grid on; title('Double Integrator Control'); legend('ct','zoh','euler','rk4','Location','best');
save_png(fig, fullfile(save_dir, 'double_integrator_plot.png'));

states = readmatrix(fullfile(save_dir, 'unicycle3_states.csv'));
fig = figure('Visible','off','Position',[100 100 520 500]); hold on; grid on; axis equal;
plot(states(:,2), states(:,3), 'k', 'LineWidth',1.5); plot(5,5,'r*','MarkerSize',8);
xlabel('x'); ylabel('y'); title('Unicycle3 Path');
save_png(fig, fullfile(save_dir, 'unicycle3_path.png'));

states = readmatrix(fullfile(save_dir, 'bicycle_dynamic_states.csv'));
fig = figure('Visible','off','Position',[100 100 560 480]); hold on; grid on; axis equal;
plot(states(:,2), states(:,3), 'k', 'LineWidth',1.5); plot(5,5,'r*','MarkerSize',8);
xlabel('X'); ylabel('Y'); title('Dynamic Bicycle Path');
save_png(fig, fullfile(save_dir, 'bicycle_dynamic_path.png'));
end

function plot_lecture3_figures(save_dir)
states = readmatrix(fullfile(save_dir, 'random_joint_states.csv'));
fig = figure('Visible','off','Position',[100 100 700 420]); hold on; grid on;
plot(states(:,1), states(:,2:end), 'LineWidth',1.1);
xlabel('t [s]'); ylabel('joint position'); title('Joint-Space Single Integrator');
save_png(fig, fullfile(save_dir, 'random_joint_motion.png'));

joint = readmatrix(fullfile(save_dir, 'joint_servo_states.csv'));
cart = readmatrix(fullfile(save_dir, 'cartesian_servo_states.csv'));
fig = figure('Visible','off','Position',[100 100 840 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on; axis equal;
plot(joint(:,5), joint(:,6), 'k', 'LineWidth',1.5); plot(joint(end,5), joint(end,6), 'r*');
title('Joint-Space Servo'); xlabel('x'); ylabel('y');
nexttile; hold on; grid on; axis equal;
plot(cart(:,5), cart(:,6), 'k', 'LineWidth',1.5); plot(1.6, 0.4, 'r*');
title('Cartesian Servo'); xlabel('x'); ylabel('y');
save_png(fig, fullfile(save_dir, 'joint_vs_cartesian.png'));
end

function plot_ilc_figures(save_dir, trajectory_name, summary_name)
states = readmatrix(fullfile(save_dir, 'states.csv'));
controls = readmatrix(fullfile(save_dir, 'controls.csv'));
errors = readmatrix(fullfile(save_dir, 'error_norm.csv'));
fig = figure('Visible','off','Position',[100 100 700 520]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on;
plot(states(:,1), states(:,4), '--', 'Color', [0.4 0.4 0.4], 'LineWidth',1.2);
plot(states(:,1), states(:,2), 'r', 'LineWidth',1.5);
title('Output'); legend('reference','final','Location','best');
nexttile; hold on; grid on;
plot(controls(:,1), controls(:,2), 'r', 'LineWidth',1.5);
title('Final Input'); xlabel('t');
save_png(fig, fullfile(save_dir, trajectory_name));

fig = figure('Visible','off','Position',[100 100 540 340]); hold on; grid on;
plot(errors(:,1), errors(:,2), 'o-', 'LineWidth',1.3);
xlabel('iteration'); ylabel('error norm'); title('ILC Error');
if strcmp(summary_name, 'time_domain_ilc_summary.png') && isfile(fullfile(save_dir, 'learning_gain.csv'))
    gain = readmatrix(fullfile(save_dir, 'learning_gain.csv'));
    clf(fig); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile; plot(errors(:,1), errors(:,2), 'o-', 'LineWidth',1.3); grid on; title('ILC Error');
    nexttile; imagesc(gain(:,2:end)); colorbar; title('Learning Gain');
end
save_png(fig, fullfile(save_dir, summary_name));
end

function plot_kalman_figure(save_dir)
data = readmatrix(fullfile(save_dir, 'estimates.csv'));
fig = figure('Visible','off','Position',[100 100 700 540]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on;
plot(data(:,1), data(:,2), 'k'); plot(data(:,1), data(:,3), 'r'); plot(data(:,1), data(:,4), 'b'); plot(data(:,1), data(:,5), '--b');
title('Kalman Estimates'); legend('state','measurement','kf','kfss','Location','best');
nexttile; hold on; grid on;
plot(data(:,1), data(:,6)); plot(data(:,1), data(:,7)); title('Posterior Variance'); legend('kf','kfss','Location','best');
save_png(fig, fullfile(save_dir, 'kalman_filter.png'));
end

function plot_ekf_ukf_figure(save_dir)
states = readmatrix(fullfile(save_dir, 'states.csv'));
ekf = readmatrix(fullfile(save_dir, 'ekf_estimates.csv'));
ukf = readmatrix(fullfile(save_dir, 'ukf_estimates.csv'));
ekfz = readmatrix(fullfile(save_dir, 'ekf_covariance.csv'));
ukfz = readmatrix(fullfile(save_dir, 'ukf_covariance.csv'));
fig = figure('Visible','off','Position',[100 100 700 620]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
for idx = 1:2
    nexttile; hold on; grid on;
    plot(states(:,1), states(:,idx+1), 'k'); plot(ekf(:,1), ekf(:,idx+1), '--b'); plot(ukf(:,1), ukf(:,idx+1), 'b');
    title("State " + idx); legend('state','ekf','ukf','Location','best');
end
nexttile; hold on; grid on;
plot(ekfz(:,1), ekfz(:,2), '--b'); plot(ukfz(:,1), ukfz(:,2), 'b'); title('Covariance Entry'); legend('ekf z11','ukf z11','Location','best');
save_png(fig, fullfile(save_dir, 'ekf_ukf.png'));
end

function plot_rls_sgd_figure(save_dir)
states = readmatrix(fullfile(save_dir, 'states.csv'));
params = readmatrix(fullfile(save_dir, 'parameter_estimates.csv'));
fig = figure('Visible','off','Position',[100 100 700 600]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile; plot(states(:,1), states(:,2), 'k'); grid on; title('State');
nexttile; hold on; grid on; plot(params(:,1), params(:,2), 'k'); plot(params(:,1), params(:,3), 'r'); plot(params(:,1), params(:,4), 'b'); title('Parameter'); legend('true','rls','sgd','Location','best');
nexttile; hold on; grid on; plot(params(:,1), params(:,5), 'r'); plot(params(:,1), params(:,6), '--b'); title('Learning Rate'); legend('rls','sgd','Location','best');
save_png(fig, fullfile(save_dir, 'rls_sgd.png'));
end

function plot_separation_figure(save_dir)
states = readmatrix(fullfile(save_dir, 'states.csv'));
measurements = readmatrix(fullfile(save_dir, 'measurements.csv'));
estimates = readmatrix(fullfile(save_dir, 'estimates.csv'));
fig = figure('Visible','off','Position',[100 100 820 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on;
plot(states(:,2), states(:,3), 'k'); plot(measurements(:,2), measurements(:,3), 'r.'); plot(estimates(:,2), estimates(:,3), 'b');
title('Position Phase Plot'); legend('state','measurement','estimate','Location','best');
nexttile; hold on; grid on;
plot(states(:,2:end), 'Color', [0 0 0 0.45]); plot(estimates(:,2:end), 'b');
title('States And Estimates');
save_png(fig, fullfile(save_dir, 'separation_principle.png'));
end

function plot_mrac_figure(save_dir)
states = readmatrix(fullfile(save_dir, 'states.csv'));
analysis = readmatrix(fullfile(save_dir, 'analysis.csv'));
fig = figure('Visible','off','Position',[100 100 700 400]); hold on; grid on;
plot(states(:,1), states(:,2), 'k'); plot(states(:,1), states(:,3), 'Color',[0.35 0.35 0.35]);
plot(analysis(:,1), analysis(:,2), 'r'); plot(analysis(:,1), analysis(:,3), '-*b'); plot(analysis(:,1), analysis(:,4), '-*c');
title('Model Reference Adaptive Control'); legend('x1','x2','value','posterior error 1','posterior error 2','Location','best');
save_png(fig, fullfile(save_dir, 'mrac.png'));
end

function plot_value_learning_figure(save_dir)
rms = readmatrix(fullfile(save_dir, 'rms_errors.csv'));
lqr = readmatrix(fullfile(save_dir, 'lqr_states.csv'));
mc = readmatrix(fullfile(save_dir, 'mc_final_states.csv'));
sarsa = readmatrix(fullfile(save_dir, 'sarsa_final_states.csv'));
qlearning = readmatrix(fullfile(save_dir, 'qlearning_final_states.csv'));
fig = figure('Visible','off','Position',[100 100 700 560]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on;
plot(sarsa(:,1), sarsa(:,2)); plot(qlearning(:,1), qlearning(:,2)); plot(mc(:,1), mc(:,2)); plot(lqr(:,1), lqr(:,2), '--k');
title('Final Episode States'); legend('sarsa','qlearning','mc','lqr','Location','best');
nexttile; hold on; grid on;
plot(rms(:,1), rms(:,2)); plot(rms(:,1), rms(:,3)); plot(rms(:,1), rms(:,4));
title('Weight Error'); legend('sarsa','qlearning','mc','Location','best');
save_png(fig, fullfile(save_dir, 'value_learning.png'));
end

function plot_pg_ilc_figure(save_dir)
metrics = readmatrix(fullfile(save_dir, 'policy_gradient_metrics.csv'));
errors = readmatrix(fullfile(save_dir, 'ilc_error_norm.csv'));
lqr = readmatrix(fullfile(save_dir, 'lqr_states.csv'));
ilc = readmatrix(fullfile(save_dir, 'ilc_final_states.csv'));
fig = figure('Visible','off','Position',[100 100 700 560]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on;
plot(metrics(:,1), metrics(:,2), 'b'); plot(metrics(:,1), metrics(:,3), 'r'); plot(metrics(:,1), metrics(:,4), 'k');
title('Policy Parameter Error'); legend('reinforce','reinforce baseline','actor critic','Location','best');
nexttile; plot(errors(:,1), errors(:,2), 'o-'); grid on; title('ILC Error');
save_png(fig, fullfile(save_dir, 'pg_ilc_comparison.png'));

fig = figure('Visible','off','Position',[100 100 600 380]); hold on; grid on;
plot(lqr(:,1), lqr(:,2), '--k'); plot(ilc(:,1), ilc(:,2), 'r');
title('ILC Final State'); legend('lqr','ilc final','Location','best');
save_png(fig, fullfile(save_dir, 'ilc_final_state.png'));
end

function plot_actor_critic_figure(save_dir)
metrics = readmatrix(fullfile(save_dir, 'policy_gradient_metrics.csv'));
lqr = readmatrix(fullfile(save_dir, 'lqr_states.csv'));
fig = figure('Visible','off','Position',[100 100 700 560]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on;
plot(metrics(:,1), metrics(:,2), 'b'); plot(metrics(:,1), metrics(:,3), 'r'); plot(metrics(:,1), metrics(:,4), 'k');
title('Policy Parameter Error'); legend('reinforce','reinforce baseline','actor critic','Location','best');
nexttile; hold on; grid on;
plot(lqr(:,1), lqr(:,2), '--k'); title('LQR State Reference'); legend('lqr','Location','best');
save_png(fig, fullfile(save_dir, 'actor_critic.png'));
end

function save_png(fig, output_path)
ensure_dir(fileparts(output_path));
try
    exportgraphics(fig, output_path, 'Resolution', 160);
catch
    saveas(fig, output_path);
end
close(fig);
end

function remove_legacy_extra_figures(save_root)
paths = [
    string(fullfile(save_root, 'lecture9', 'lecture9_state_plot.png'))
    string(fullfile(save_root, 'lecture10', 'reference_ilqr_state_plot.png'))
    string(fullfile(save_root, 'lecture10_mpc', 'lecture10_mpc_state_plot.png'))
];
for i = 1:numel(paths)
    if isfile(paths(i)), delete(paths(i)); end
end
end

function ensure_dir(path)
if ~exist(path, 'dir'), mkdir(path); end
end

function write_csv(path, header, data)
fid = fopen(path, 'w');
if fid < 0, error('Could not open %s for writing.', path); end
fprintf(fid, '%s\n', strjoin(string(header), ','));
fclose(fid);
writematrix(data, path, 'WriteMode', 'append');
end

function write_json(path, payload)
try
    encoded = jsonencode(payload, 'PrettyPrint', true);
catch
    encoded = jsonencode(payload);
end
fid = fopen(path, 'w');
if fid < 0, error('Could not open %s for writing.', path); end
fprintf(fid, '%s\n', encoded);
fclose(fid);
end
