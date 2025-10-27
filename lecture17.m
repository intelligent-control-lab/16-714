% 16-714 Advanced Control for Robotics
% Script for Lecture 17: RLS 

%% Problem Setup
sys.x0 = 1;
sys.W = 0;        
sys.Bw = 1;

sys.fnominal = @(x, p) p*sin(x); 
sys.f        = @(x, p) sys.fnominal(x,p) + sys.Bw * sqrt(sys.W) * randn();
sys.grad     = @(x, p) sin(x);

p = @(x,t) 1 + 0.05*sin(t/10);   % time-varying parameter (acts as control u)
sys.p0hat = 1;
sys.H0 = 0.5;
sys.lambda = 0.5;

sys.N  = 100;
sys.dt = 1;

%% New roll_out usage
% Model wrapper for DT 'Direct' integration: x_{k+1} = f(x_k, u_k)
model.name = 'scalar_sin_dt';
model.f_dt    = @(x,u,t) sys.f(x,u);   % DT one-step map

% Controller supplies the time-varying parameter as u = p(x,t)
ctrl = @(x,t) p(x,t);

% Simulation settings
sim = struct('type','DT', 'K', sys.N, 'dt', sys.dt, 'integrator','Direct');

% Optional flags (leave empty if not needed)
opts = struct();  % e.g., opts.timevarying = true;

% Simulate
[tlist, xlist, ulist, log] = roll_out(model, ctrl, sys.x0, sim, opts);

%% Plot states
figure(1); clf;
subplot(3,1,1); hold on;
plot(xlist(1,:),'k','LineWidth',1.5);
xlabel('k'); ylabel('x_k'); title('State trajectory');

%% Estimators
rls = estimator('RLS', sys);
sgd = estimator('SGD', sys);

%% Adaptation
for k = 1:sys.N
    rls.phat(:,k+1) = rls.update_p(rls.phat(:,k), xlist(:,k), xlist(:,k+1), rls.H(:,:,k));
    rls.H(:,:,k+1)  = rls.update_H(rls.H(:,:,k), rls.phat(:,k), xlist(:,k));
    sgd.phat(:,k+1) = sgd.update_p(sgd.phat(:,k), xlist(:,k), xlist(:,k+1));
end

%% Plot parameter (u = p) and estimates
subplot(3,1,2); hold on;
plot(ulist(1,:),'k','LineWidth',1.5);         % ground-truth p_k (as ulist)
plot(rls.phat(1,:),'r','LineWidth',1.2);
plot(sgd.phat(1,:),'b','LineWidth',1.2);
xlabel('k'); ylabel('p_k');
legend('Ground truth p','RLS','SGD','Location','Best');
title('Parameter estimates');

%% Plot learning rates
subplot(3,1,3); hold on;
plot(permute(1./rls.H(1,1,:),[3,2,1]),'r','LineWidth',1.2);
plot((1/sys.H0).*ones(1,sys.N+1),'b--','LineWidth',1.0);
xlabel('k'); ylabel('step size');
legend('RLS learning rate','SGD learning rate baseline','Location','Best');
title('Learning rate comparison');
