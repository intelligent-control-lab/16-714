% 16-714 Advanced Control for Robotics
% Script for Lecture 18: Separation Principle

%% Problem Parameters
dim = 4;
dyn_name = "double_integrator";
sys.x0    = [10;5;5;0];
sys.X0    = 0.001*eye(dim);
sys.x0hat = sys.x0 + sqrt(sys.X0)*randn(dim,1); % a posteriori initial estimate
sys.dt    = 0.2;

% Get A,B via resolve_dynamics.getAB
dynspec   = resolve_dynamics(dyn_name);
[sys.A, sys.B] = dynspec.getAB('DT', struct('dim', dim, 'dt', sys.dt, 'dmode','ZOH'));

% Noise, plant, and measurement
sys.W  = 0.1  * eye(dim/2);     % process noise on acceleration channels
sys.V  = 0.01 * eye(dim/2);    % measurement noise on position channels
sys.Bw = sys.B;                 % process noise input matrix

sys.fnominal = @(x,u) sys.A * x + sys.B * u;
sys.f        = @(x,u) sys.fnominal(x,u) + sys.Bw * (sqrtm(sys.W) * randn(dim/2,1));
sys.hnominal = @(x) x(1:2);
sys.h        = @(x) sys.hnominal(x) + sqrtm(sys.V) * randn(dim/2,1);
sys.C        = [eye(2) zeros(2)];

% LQR cost
sys.Q = [eye(dim/2) zeros(dim/2); zeros(dim/2) zeros(dim/2)];
sys.R = eye(dim/2);

Tmax  = 50;
sys.N = Tmax;

%% Synthesis (controller) and Estimator (steady-state KF)
c   = synthesis('LQR', sys);
est = estimator('KFss', sys);  % provides update_x(xhat,u,y,Z) and update_Z(Z); est.Z preallocated

%% (Option 1) Simulate with roll_out only

% Model: DT map via Direct (we include process noise in fd)
model = struct();
model.name = 'double_int_sep_dt';
model.f_dt   = @(x,u) sys.f(x,u);   % used when sim.integrator='Direct'

% Controller: estimate feedback using current xhat (provided by opts.observer)
ctrl = struct('mode','estimate', 'pi', @(xhat,t) c.u(xhat,t));

% Observer: wrap the KF with persistent Z so roll_out can call update(xhat,u,y,t)
observer.xhat0 = sys.x0hat;
observer.update = estimator_update_handle(est);

% Sensor: measurement function y = h(x) (u unused but included to match signature)
sensor = @(x,t) sys.h(x);

% Simulation settings
sim = struct('type','DT', 'K', sys.N, 'dt', sys.dt, 'integrator','Direct');

% opts: plug in sensor/observer and log what we need
opts = struct();
opts.sensor      = sensor;
opts.observer    = observer;
opts.log_fields  = {'y','xhat'};
opts.verbose     = false;

[tlist, xlist, ulist, log] = roll_out(model, ctrl, sys.x0, sim, opts);

% Pull logs
ylist     = log.y;        % size 2 x N
xhat_hist = log.xhat;     % size 4 x N

% %% (Option 2) Simulate with explicit for hoop
% % Model: DT map via Direct (we include process noise in fd)
% model = struct();
% model.name = 'double_int_sep_dt';
% model.f_dt   = @(x,u) sys.f(x,u);   % used when sim.integrator='Direct'
% 
% xlist = zeros(dim,Tmax); xlist(:,1) = sys.x0; 
% ylist = zeros(dim/2,Tmax); ylist(:,1) = sys.x0hat(1:2); 
% for k = 1:Tmax 
%     ulist(:,k) = c.u(est.xhat(:,k),(k-1)*sys.dt); 
%     xlist(:,k+1) = step(xlist(:,k), ulist(:,k), sys.dt, model, 'Direct'); 
%     ylist(:,k+1) = sys.h(xlist(:,k+1)); 
%     est.xhat(:,k+1) = est.update_x(est.xhat(:,k),ulist(:,k),ylist(:,k+1),est.Z(:,:,k)); 
%     est.Z(:,:,k+1) = est.update_Z(est.Z(:,:,k)); 
% end
% 
% xhat_hist = est.xhat;
%% Plots
figure(1); clf; hold on;
plot(xlist(1,:), xlist(2,:), 'k');
plot(ylist(1,:), ylist(2,:), 'r');
plot(xhat_hist(1,:), xhat_hist(2,:), 'b');
legend("Ground truth x", "Measurement y", "A posteriori estimate xhat (KF)");
title('Separation Principle (position phase plot)');
xlabel('x_1'); ylabel('x_2');

figure(2); clf;
for i = 1:dim
    subplot(4,1,i); hold on;
    plot(xlist(i,:), 'k');
    if i <= 2, plot(ylist(i,:), 'r'); end
    plot(xhat_hist(i,:), 'b');
    if i == 1
        legend("Ground truth x", "Measurement y", "A posteriori estimate xhat (KF)", ...
               'Location','best');
        title('States, measurements, and KF estimates');
    end
    xlabel('k');
end


