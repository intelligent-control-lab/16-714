% 16-714 Advanced Control for Robotics
% Script for Lecture 15: Kalman Filter

%% Problem Setup
% x' = 0.5*x + u + w;
% y = x + v
sys.A = 0.5;
sys.B = 1;
sys.Bw = 1;
sys.C = 1;
sys.x0 = 0;
sys.N = 101;
sys.dt = 1;
% Noise models
sys.X0 = 1; sys.x0hat = sys.x0 + randn * sqrt(sys.X0);%A posteriori
sys.W = 0.01;
sys.V = 0.1;
% Dynamics
sys.fnominal = @(x,u) sys.A*x+sys.B*u;
sys.f = @(x,u) sys.fnominal(x,u) + sys.Bw * sqrt(sys.W) * randn();
sys.hnominal = @(x) sys.C * x;
sys.h = @(x) sys.hnominal(x) + sqrt(sys.V) * randn();
sys.simmode = 'Direct'; % The dynamics are directly specified in discrete time
% Controller
u = @(x,t) sin(t/(sys.N/2/pi));
% Simulate
[~, xlist, ulist, ylist] = roll_out(sys, u, sys.x0, 'DT', sys.N, sys.dt, 'Direct');
% Plot
figure(1);clf;hold on;
plot(xlist,'k')
plot(ylist,'r')
%% Estimate States
est = estimator('KF',sys);
estss = estimator('KF',sys);

%% Filtering
for k = 1:sys.N
    est.xhat(:,k+1) = est.update_x(est.xhat(:,k),ulist(:,k),ylist(:,k+1),est.Z(:,:,k));
    est.Z(:,:,k+1) = est.update_Z(est.Z(:,:,k));

    estss.xhat(:,k+1) = estss.update_x(estss.xhat(:,k),ulist(:,k),ylist(:,k+1),estss.Z(:,:,k));
    estss.Z(:,:,k+1) = estss.update_Z(estss.Z(:,:,k));
end
display('Measurement Error:')
display(norm(ylist-xlist))
display('Estimation Error (KF):')
display(norm(est.xhat-xlist));
display('Estimation Error (KFss):')
display(norm(estss.xhat-xlist));

% Plot
figure(1);
plot(est.xhat,'b')
plot(estss.xhat,'--b')
legend("Ground truth x", "Measurement y", "A posteriori estimate xhat (KF)", "A posteriori estimate xhat (KFss)")
figure(2);clf;hold on;
plot(permute(est.Z,[3,2,1]))
plot(permute(estss.Z,[3,2,1]))
legend("A posterior variance (KF)", "A posterior variance (KFss)")