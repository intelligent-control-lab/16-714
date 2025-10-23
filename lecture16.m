% 16-714 Advanced Control for Robotics
% Script for Lecture 16: EKF and UKF

%% Problem Setup
problem = 4;
switch problem
    case 1
        sys.x0 = [1;2];
        sys.X0 = 0.001*[1 0;0 2]; sys.x0hat = sys.x0 + sqrt(sys.X0)*randn(2,1); % A posteriori
        sys.W = 0.001; sys.V = 0.001; sys.Bw = 1;
        sys.fnominal = @(x,u) x + u;
        sys.f = @(x,u) sys.fnominal(x,u) + sys.Bw * sqrt(sys.W) * randn();
        sys.hnominal = @(x) [x(1)^2 + x(2)^2; x(1)*x(2)];
        sys.h = @(x) sys.hnominal(x) + sqrt(sys.V) * randn();
        sys.A = @(x,u) eye(2);
        sys.C = @(x) [2*x(1) 2*x(2); x(2) x(1)];

    case 2
        sys.x0 = [1;2];
        sys.X0 = 0.01*[1 0;0 2]; sys.x0hat = sys.x0 + sqrt(sys.X0)*randn(2,1);
        sys.W = 0.001; sys.V = 0.001; sys.Bw = 1;
        sys.fnominal = @(x,u) x + u;
        sys.f = @(x,u) sys.fnominal(x,u) + sys.Bw * sqrt(sys.W) * randn();
        sys.hnominal = @(x) [x(1)^2 + x(2)^2; x(1)*x(2)];
        sys.h = @(x) sys.hnominal(x) + sqrt(sys.V) * randn();
        sys.A = @(x,u) eye(2);
        sys.C = @(x) [2*x(1) 2*x(2); x(2) x(1)];

    case 3
        sys.x0 = [0.1; -1];
        sys.X0 = 0.01*[8 2; 2 3]; sys.x0hat = sys.x0 + sqrt(sys.X0)*randn(2,1);
        sys.W = 0.001; sys.V = 0.001; sys.Bw = 1;
        sys.fnominal = @(x,u) x + u;
        sys.f = @(x,u) sys.fnominal(x,u) + sys.Bw * sqrt(sys.W) * randn();
        sys.hnominal = @(x) [x(1)^2 + x(2)^2; x(1)*x(2)];
        sys.h = @(x) sys.hnominal(x) + sqrt(sys.V) * randn();
        sys.A = @(x,u) eye(2);
        sys.C = @(x) [2*x(1) 2*x(2); x(2) x(1)];

    case 4
        sys.x0 = [4; 1];
        sys.X0 = 0.01*[8 2; 2 3]; sys.x0hat = sys.x0 + sqrt(sys.X0)*randn(2,1);
        sys.W = 0.001; sys.V = 0.001; sys.Bw = 1;
        sys.fnominal = @(x,u) x + u;
        sys.f = @(x,u) sys.fnominal(x,u) + sys.Bw * sqrt(sys.W) * randn();
        sys.hnominal = @(x) [x(1)^4 * x(2); x(1) + x(2)^5];
        sys.h = @(x) sys.hnominal(x) + sqrt(sys.V) * randn();
        sys.A = @(x,u) eye(2);
        sys.C = @(x) [4*x(1)^3*x(2)  x(1)^4;  1  5*x(2)^4];
end

sys.N  = 10;
sys.dt = 1;

%% Controller (state feedback u(x,t))
u = @(x,t) [0;0];

%% Simulation spec for new roll_out
sim = struct();
sim.type       = 'DT';
sim.K          = sys.N;
sim.dt         = sys.dt;
sim.integrator = 'Direct';   % your DT map is given directly by sys.f

% (opts is optional; include it if you use any flags like time-varying)
opts = struct();  
opts.sensor     = @(x,u,t) sys.h(x);
opts.log_fields = {'y'};     % ask roll_out to record measurements y_k

model.f_dt = @(x,u,dt) sys.f(x,u);
%% Simulate
% Expected outputs: [tlist, xlist, ulist, log]
% log.y is assumed to store measurements y_k
[tlist, xlist, ulist, log] = roll_out(model, u, sys.x0, sim, opts);

% Get measurement sequence
if isstruct(log) && isfield(log,'y')
    ylist = log.y;
else
    % Fallback: synthesize measurements from sys.h (will include noise if sys.h adds it)
    ylist = zeros(length(sys.h(sys.x0)), size(xlist,2));
    for k = 1:size(xlist,2)
        ylist(:,k) = sys.h(xlist(:,k));
    end
end

%% Plot ground truth
figure(1); clf;
subplot(2,1,1); hold on;
plot(tlist, xlist(1,:),'k'); ylabel('x_1');
subplot(2,1,2); hold on;
plot(tlist, xlist(2,:),'k'); ylabel('x_2'); xlabel('k');

%% Estimators
ekf = estimator('EKF', sys);
ukf = estimator('UKF', sys);

%% Filtering
for k = 1:sys.N-1
    ekf.xhat(:,k+1) = ekf.update_x(ekf.xhat(:,k), ulist(:,k), ylist(:,k+1), ekf.Z(:,:,k));
    ekf.Z(:,:,k+1)  = ekf.update_Z(ekf.Z(:,:,k), ekf.xhat(:,k), ulist(:,k));

    ukf.xhat(:,k+1) = ukf.update_x(ukf.xhat(:,k), ulist(:,k), ylist(:,k+1), ukf.Z(:,:,k));
    ukf.Z(:,:,k+1)  = ukf.update_Z(ukf.Z(:,:,k), ukf.xhat(:,k), ulist(:,k));
end

disp('Estimation Error (EKF):')
disp(norm(ekf.xhat - xlist(1:sys.N)));
disp('Estimation Error (UKF):')
disp(norm(ukf.xhat - xlist(1:sys.N)));

% Plot estimates
figure(1);
subplot(2,1,1); hold on;
plot(tlist(1:sys.N), ekf.xhat(1,:),'--b');
plot(tlist(1:sys.N), ukf.xhat(1,:),'b');
legend('Ground truth x_1','EKF','UKF');

subplot(2,1,2); hold on;
plot(tlist(1:sys.N), ekf.xhat(2,:),'--b');
plot(tlist(1:sys.N), ukf.xhat(2,:),'b');
legend('Ground truth x_2','EKF','UKF');

% Plot covariance entries
figure(2); clf; hold on;
legstr = strings(0);
for k = 1:2
    for l = 1:2
        plot(tlist(1:sys.N), permute(ekf.Z(k,l,:),[3,2,1]));
        plot(tlist(1:sys.N), permute(ukf.Z(k,l,:),[3,2,1]));
        legstr = [legstr, "EKF Z_{"+num2str(k)+num2str(l)+"}", "UKF Z_{"+num2str(k)+num2str(l)+"}"]; %#ok<AGROW>
    end
end
legend(legstr);

