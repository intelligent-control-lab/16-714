% 16-714 Advanced Control for Robotics
% Script for Lecture 12 ILC Frequency Domain  

%% Problem Specification
sys.dt   = 1;                          % Sampling time
sys.name = 'double_integrator';

% --- Use resolve_dynamics to obtain getAB ---
dyn = resolve_dynamics(sys.name);      % returns struct with dyn.getAB, f_ct/f_dt, etc.
[sys.A, sys.B] = dyn.getAB('DT', struct('dim',2,'dt',sys.dt,'integrator','Euler'));

% Plant TF and controller TF parts for ILC synthesis
sys.a  = [1 -2 1];
sys.b  = [sys.dt];
sys.C  = [1 0];                        % only measuring the position
sys.x0 = [0;0];

% Reference signal
ref        = @(t) sin(0.2*t);
simhorizon = 100;                      % number of steps K
niter      = 10;
sigma      = 0;                        % control disturbance; 0 if no disturbance

%% Specify feedback control
kp = 0.1; kd = 0.5;
ufb = @(x,t) kp*(ref(t) - sys.C*x) + ...
             kd*((ref(t+sys.dt)-ref(t))/sys.dt - x(2)) + ...
             sigma * randn;           % proportional + derivative on error

% Controller TF coefficients for ILC design (note: non-causal form used in analysis)
sys.c = [-kd/sys.dt, kd/sys.dt - kp];

%% New roll_out call (DT with struct)
sim = struct('type','DT', 'K', simhorizon, 'dt', sys.dt, 'integrator','Euler');

% First trajectory (Iteration 0)
[tlist, x_fb, u_fb] = roll_out(sys.name, ufb, sys.x0, sim);

%% Plot iteration 0
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

%% ILC learning gain (frequency-domain)
cILC = synthesis('ILCw', sys);

%% ILC loop
error  = zeros(niter, simhorizon+1);
ff     = zeros(niter, simhorizon+1);
e      = zeros(1, niter);
legstr = ["Reference", "Iter 0"];

for i = 1:niter
    error(i,:) = ref(tlist) - x_fb(1,:);

    if i > 1
        ff(i,:) = cILC.u(error(i,:), ff(i-1,:));
    else
        ff(i,:) = cILC.u(error(i,:), zeros(1, simhorizon+1));
    end

    % feedforward + feedback (index via t/dt + 1)
    uff = @(x,t) ff(i, round(t/sys.dt)+1) + ufb(x,t);

    [tlist, x_fb, u_fb] = roll_out(sys.name, uff, sys.x0, sim);

    subplot(2,1,1); plot(tlist, x_fb(1,:), 'Color', [i/niter, 0, 0]);
    subplot(2,1,2); plot(tlist(1:end-1), u_fb, 'Color', [i/niter, 0, 0]);

    e(i) = norm(error(i,:));
    legstr = [legstr, "Iter " + num2str(i)];
end

subplot(2,1,1); legend(legstr, 'Location','best');
%% Error Across ILC Iterations
figure(2); clf; grid on;
plot(0:niter-1, e, 'o-');
xlabel('ILC iteration'); ylabel('||e||_2');
title('Norm Error across ILC Iterations');
