% 16-714 Advanced Control for Robotics
% Script for Lecture 18: Seperation Principle

%% Problem Parameters
dim = 4; dyn = "double_integrator";
sys.x0 = [10;5;5;0];
sys.X0 = 0.001*eye(4); sys.x0hat = sys.x0 + sqrt(sys.X0)*randn(dim,1);%A posteriori
sys.dt = 0.2;

% Dynamic parameters
[sys.A, sys.B] = getAB(dyn,dim,'DT',sys.dt );
sys.W = 0.1 * eye(2);
sys.V = 0.001 * eye(2);
sys.Bw = sys.B;
sys.fnominal = @(x,u) sys.A * x + sys.B * u;
sys.f = @(x,u) sys.fnominal(x,u) + sys.Bw * sqrt(sys.W) * randn(dim/2,1);
sys.hnominal = @(x) x(1:2);
sys.h = @(x) sys.hnominal(x) + sqrt(sys.V) * randn(dim/2,1);
sys.C = [eye(2) zeros(2)];

% Cost function
sys.Q = [eye(dim/2) zeros(dim/2); zeros(dim/2) zeros(dim/2)];
sys.R = eye(dim/2);

Tmax = 50;
sys.N = Tmax;

%% Get controller and estimator seperately
c = synthesis('LQR',sys);
est = estimator('KFss',sys);

xlist = zeros(dim,Tmax); xlist(:,1) = sys.x0;
ylist = zeros(dim/2,Tmax); ylist(:,1) = sys.x0hat(1:2);
for k = 1:Tmax
    ulist(:,k) = c.u(est.xhat(:,k),(k-1)*sys.dt);
    xlist(:,k+1) = step(xlist(:,k), ulist(:,k), sys.dt, dyn, 'Euler');
    ylist(:,k+1) = sys.h(xlist(:,k+1));
    est.xhat(:,k+1) = est.update_x(est.xhat(:,k),ulist(:,k),ylist(:,k+1),est.Z(:,:,k));
    est.Z(:,:,k+1) = est.update_Z(est.Z(:,:,k));
end

%%
figure(1); clf; hold on;
plot(xlist(1,:),xlist(2,:),'k');
plot(ylist(1,:),ylist(2,:),'r');
plot(est.xhat(1,:),est.xhat(2,:),'b');
legend("Ground truth x", "Measurement y", "A posteriori estimate xhat (KF)")

figure(2); clf; hold on;
for i = 1:dim
subplot(str2num(['41',num2str(i)])); hold on
plot(xlist(i,:),'k')
if i < 3 
    plot(ylist(i,:),'r')
end
plot(est.xhat(i,:),'b')
if i == 1
legend("Ground truth x", "Measurement y", "A posteriori estimate xhat (KF)")
end
end

