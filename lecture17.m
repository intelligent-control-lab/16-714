% 16-714 Advanced Control for Robotics
% Script for Lecture 17: RLS

%% Problem Setup
sys.x0 = 1;
sys.W = 0;        
sys.Bw = 1;

sys.fnominal = @(x, p) p*sin(x); 
sys.f = @(x, p) sys.fnominal(x,p) + sys.Bw * sqrt(sys.W) * randn();
sys.grad = @(x, p) sin(x);

p = @(x,t) 1+0.05*sin(t/10);%0.99+0.01*t;
sys.p0hat = 1;
sys.H0 = 0.5;
sys.lambda = 0.5;

sys.N = 100;
sys.dt = 1;

sys.simmode = 'Direct'; % The dynamics are directly specified in discrete time
% Simulate (note the time varying parameter is injected to the system as the control u)
[~, xlist,plist] = roll_out(sys, p, sys.x0, 'DT', sys.N, sys.dt, 'Direct');
% Plot
figure(1);clf;
subplot(311);hold on;
plot(xlist(1,:),'k')

%% Estimate Parameters
rls = estimator('RLS', sys);
sgd = estimator('SGD', sys);

%% Adaptation
for k = 1:sys.N
    rls.phat(:,k+1) = rls.update_p(rls.phat(:,k),xlist(:,k),xlist(:,k+1),rls.H(:,:,k));
    rls.H(:,:,k+1) = rls.update_H(rls.H(:,:,k), rls.phat(:,k), xlist(:,k));
    sgd.phat(:,k+1) = sgd.update_p(sgd.phat(:,k),xlist(:,k),xlist(:,k+1));
end

% Plot
figure(1);
subplot(312);hold on;
plot(plist,'k')
plot(rls.phat(1,:),'r')
plot(sgd.phat(1,:),'b')
legend("Ground truth x", "RLS", "SGD")

subplot(313);hold on;
plot(permute(1/rls.H(1,1,:),[3,2,1]),'r')
plot(1/sys.H0.*ones(1,sys.N),'b')
legend("RLS learning rate", "SGD learning rate")