% 16-714 Advanced Control for Robotics
% Script for Lecture 15: Kalman Filter
%% Problem Setup (system + noise)
% DT model:
%   x_{k+1} = A x_k + B u_k + w_k
%   y_k     = C x_k + v_k
sys.A  = 0.5;
sys.B  = 1;
sys.Bw = 1;         % process noise gain
sys.C  = 1;

sys.N  = 101;       % K steps -> x has K+1 samples, u has K samples, y has K samples
sys.dt = 1;

% initial condition and covariances
sys.x0    = 0;
sys.X0    = 1;                          % prior covariance on x0
sys.x0hat = sys.x0 + randn * sqrt(sys.X0);

% noise (scalar for this 1D example)
sys.W = 0.01;       % process noise variance
sys.V = 0.1;        % measurement noise variance

%% roll_out model (DT "Direct" map via step.m)
% Provide a DT update for the 'direct' integrator. step() will call f_dt.
model.f_dt = @(x,u,dt) sys.A*x + sys.B*u + sys.Bw*sqrt(sys.W)*randn();
sys.fnominal = @(x,u)  sys.A*x + sys.B*u;

%% Controller: state feedback u = sin( t / (N/2/pi) )
% roll_out passes tk=k*dt, so we can use t directly.
ctrl = @(x,t) sin( t / (sys.N/2/pi) );

%% Measurement model (logged inside roll_out)
sensor = @(x,u,t) sys.C*x + sqrt(sys.V)*randn();

%% Simulation config (DT, direct)
sim_cfg = struct('type','DT', 'K', sys.N, 'dt', sys.dt, 'integrator','Direct');

%% Logging options
opts = struct();
opts.sensor     = sensor;
opts.log_fields = {'y'};     % ask roll_out to record measurements y_k

%% Simulate (5-arg call)
[tlist, xlist, ulist, log] = roll_out(model, ctrl, sys.x0, sim_cfg, opts);

% Shapes:
%   tlist: 1×(K+1) with tlist(k+1)=k*dt, k=0..K
%   xlist: 1×(K+1)
%   ulist: 1×K      (u_0..u_{K-1})
%   log.y: 1×K      (y_0..y_{K-1})

%% Plot measurements vs ground truth
figure(1); clf; hold on;
plot(tlist, xlist, 'k', 'DisplayName','Ground truth x');
plot(tlist(1:numel(ulist)), log.y, 'r', 'DisplayName','Measurement y');  % y_k aligns with t_k, k=0..K-1
xlabel('k'); ylabel('value'); grid on;

%% Estimators
% estimator('KF', sys) should read A,B,C,W,V,X0,x0hat from sys
est   = estimator('KF', sys);
estss = estimator('KFss', sys);   

%% Filtering loop
% roll_out gives y_k (not y_{k+1}), so we use:
%   given (xhat_k, Z_k) and input u_k, measurement y_k -> compute (xhat_{k+1}, Z_{k+1})
for k = 1:K
    est.xhat(:,k+1)   = est.update_x(  est.xhat(:,k),   uk(k), yk(k),   est.Z(:,:,k));
    est.Z(:,:,k+1)    = est.update_Z(  est.Z(:,:,k));

    estss.xhat(:,k+1) = estss.update_x(estss.xhat(:,k), uk(k), yk(k), estss.Z(:,:,k));
    estss.Z(:,:,k+1)  = estss.update_Z(estss.Z(:,:,k));
end

%% Error reports
% Compare on time indices that exist for both signals
fprintf('Measurement Error (||y - x||_2):        %g\n', norm(yk - xk(1:K)));
fprintf('Estimation Error (KF,   ||xhat - x||):   %g\n', norm(est.xhat(:)  - xk(:)));
fprintf('Estimation Error (KFss, ||xhat - x||):   %g\n', norm(estss.xhat(:) - xk(:)));

%% Plot estimates
figure(1);
plot(t, est.xhat,  'b',   'DisplayName','A posteriori xhat (KF)');
plot(t, estss.xhat,'--b', 'DisplayName','A posteriori xhat (KFss)');
legend('Location','best');

%% Plot posterior variances (1D so these are scalars)
figure(2); clf; hold on;
plot(0:K, squeeze(est.Z),   'DisplayName','A posterior variance (KF)');
plot(0:K, squeeze(estss.Z), 'DisplayName','A posterior variance (KFss)');
xlabel('k'); ylabel('variance'); grid on; legend('Location','best');
