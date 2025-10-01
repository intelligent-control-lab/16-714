% 16-714 Advanced Control for Robotics
% Script for Lecture 13: ILC Time Domain 

%% Problem Specification
sys.dt   = 1;                       % Sampling time
sys.name = 'double_integrator';
sys.N    = 100;

% Use resolve_dynamics to obtain A,B (and keep for possible f_dt use)
dyn = resolve_dynamics(sys.name);
[sys.A, sys.B] = dyn.getAB('DT', struct('dim',2,'dt',sys.dt,'integrator','Euler'));

% Plant TF (for ILC synthesis)
sys.a = [1 -2 1];
sys.b = [sys.dt];
sys.C = [1 0];                      % measuring position
sys.x0 = [0;0];

% Reference signal
ref        = @(t) sin(0.2*t);
simhorizon = sys.N;
niter      = 10;
sigma      = 0;                  % control disturbance; 0 if no disturbance

%% Specify feedback control
kp = 0.1; kd = 0.5;
sys.K = [kp kd];                    % (used by time-domain ILC synthesis if needed)
sys.c = [-kd/sys.dt, kd/sys.dt - kp];  % (for freq-domain analysis)

ufb = @(x,t) kp*(ref(t) - sys.C*x) + ...
             kd*((ref(t+sys.dt)-ref(t))/sys.dt - x(2)) + ...
             sigma * randn;

%% New roll_out call (DT with struct)
sim = struct('type','DT','K',simhorizon,'dt',sys.dt,'integrator','Euler');

% First trajectory (Iteration 0)
[tlist, x_fb, u_fb] = roll_out(sys.name, ufb, sys.x0, sim);

% Plot iteration 0
figure(1); clf;
subplot(2,1,1); hold on; grid on;
plot(tlist, ref(tlist), '--');
plot(tlist, x_fb(1,:), 'k');
xlabel('time (s)'); ylabel('position');
title('Output'); xlim([tlist(1), tlist(end)]);

subplot(2,1,2); hold on; grid on;
plot(tlist(1:end-1), u_fb);
xlabel('time (s)'); ylabel('input');
title('Input'); xlim([tlist(1), tlist(end-1)]);

%% Specify the learning gain (Time-domain ILC)
cILCt = synthesis('ILCt', sys);

%% ILC loop
error  = zeros(niter, simhorizon+1);
ff     = zeros(niter, simhorizon);      % feedforward is defined on k = 0..K-1
e      = zeros(1, niter);
legstr = ["Reference", "Iter 0"];

for i = 1:niter
    % tracking error from previous rollout
    error(i,:) = ref(tlist) - x_fb(1,:);

    % update FF using time-domain learning filter
    if i > 1
        ff(i,1:simhorizon) = cILCt.u(error(i,:), ff(i-1,1:simhorizon));
    else
        ff(i,1:simhorizon) = cILCt.u(error(i,:), zeros(1, simhorizon));
    end

    % combined controller: FF(k) + FB(x,t); index FF by k = t/dt
    uff = @(x,t) ff(i, round(t/sys.dt)+1) + ufb(x,t);

    % rollout with updated controller
    [tlist, x_fb, u_fb] = roll_out(sys.name, uff, sys.x0, sim);

    % plots for iteration i
    subplot(2,1,1);
    plot(tlist, x_fb(1,:), 'Color', [i/niter, 0, 0]);
    subplot(2,1,2);
    plot(tlist(1:end-1), u_fb, 'Color', [i/niter, 0, 0]);

    % error norm
    e(i) = norm(error(i,:));
    legstr = [legstr, "Iter " + num2str(i)];
end

subplot(2,1,1); legend(legstr, 'Location','best');

figure(2); clf; grid on;
plot(0:niter-1, e, 'o-');
xlabel('ILC iteration'); ylabel('||e||_2');
title('Norm Error across ILC Iterations');

%% Compare time-domain ILC vs frequency-domain ILC (matrix form of L)
cILCw = synthesis('ILCw', sys);

% Build Toeplitz-like L from FIR b-coeffs of learning filter
L_matrix_form = zeros(simhorizon, simhorizon+1);
l = length(cILCw.L.b);
for i = 1:simhorizon
    j = 0;
    while j < l && i + j <= simhorizon + 1
        L_matrix_form(i, i+j) = cILCw.L.b(l-j);
        j = j + 1;
    end
end
% Example check:
% ff(1,1:simhorizon)' - L_matrix_form * error(1,:)'

figure(10); clf;
subplot(1,2,1);
imagesc(cILCt.L); axis image; colorbar;
title('L matrix (time-domain ILC)');

subplot(1,2,2);
imagesc(L_matrix_form); axis image; colorbar;
title('L matrix (freq-domain ILC as matrix)');
