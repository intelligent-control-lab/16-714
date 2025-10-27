% 16-714 Advanced Control for Robotics
% Script for Lecture 23: Actor Critic

%% Problem Setup
sys.A = 1; sys.B = 0.5;
sys.f_dt = @(x,u,t) sys.A*x + sys.B*u;
sys.Q = 1; sys.R = 1;
sys.h = @(x,u) (x'*sys.Q*x + u'*sys.R*u)/2;
sys.tolerance = 0.01;
sys.kmax = 9; sys.dt = 1; sys.x0 = 1;
sys.tcondition = @(x,t) norm(x) < sys.tolerance || t > sys.kmax * sys.dt;

%% LQR solution
c = synthesis('LQR',sys);

simDT  = struct( ...
    'type','DT', ...
    'K',   sys.kmax, ...                 % nominal horizon steps
    'dt',  sys.dt, ...
    'integrator','Direct', ...
    'stop', @(xk,tk,log) sys.tcondition(xk,tk) );  % early stop condition
optsDT = struct(); 

[~, x_lqr] = roll_out(sys, c.u, sys.x0, simDT, optsDT);

%% Learning paramters
sys.W0 = 3; sys.delta = 1;

sys.policy.name = 'Gaussian';
sys.policy.init.mu = -0.1*eye(size(sys.B'));
sys.policy.init.sigma = 0.1*eye(size(sys.B'));

n_ep = 200;
alpha = 0.01; %policy
alpha_baseline = 1; %value
%% REINFORCE
reinforce = PGagent(sys);
rms_reinforce = zeros(1,n_ep);
std_reinforce = zeros(1,n_ep);
figure(1);clf;subplot(311);hold on;
for episode = 1:n_ep
    x_list = reinforce.update('REINFORCE',alpha,sys);
    plot(0:length(x_list)-1, x_list,'color',[0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep])
    rms_reinforce(episode) = norm(reinforce.theta.mu + c.K);
    std_reinforce(episode) = reinforce.theta.sigma;
end
plot(0:length(x_lqr)-1, x_lqr,'--k')
title("REINFORCE")
box on;
%% REINFORCE with Baseline
reinforce_baseline = PGagent(sys);
rms_reinforce_baseline = zeros(1,n_ep);
std_reinforce_baseline = zeros(1,n_ep);
subplot(312);hold on;
alphas.p = alpha;
alphas.v = alpha_baseline;
for episode = 1:n_ep
    x_list = reinforce_baseline.update('REINFORCEbl',alphas,sys);
    plot(0:length(x_list)-1, x_list,'color',[0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep])
    rms_reinforce_baseline(episode) = norm(reinforce_baseline.theta.mu + c.K);
    std_reinforce_baseline(episode) = reinforce_baseline.theta.sigma;
end
plot(0:length(x_lqr)-1, x_lqr,'--k')
title("REINFORCE with Baseline")
box on;
%% Actor Critic
actor_critic = PGagent(sys);
rms_ac = zeros(1,n_ep);
std_ac = zeros(1,n_ep);
subplot(313);hold on;
for episode = 1:n_ep
    x_list = actor_critic.update('ActorCritic',alphas,sys);
    plot(0:length(x_list)-1, x_list,'color',[0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep])
    rms_ac(episode) = norm(actor_critic.theta.mu + c.K);
    std_ac(episode) = actor_critic.theta.sigma;
end
plot(0:length(x_lqr)-1, x_lqr,'--k')
title("Actor Critic")
box on;
%%
figure(2);clf;hold on
plot(rms_reinforce,'b');
plot(rms_reinforce_baseline,'r');
plot(rms_ac,'k')
% plot(rms_reinforce + std_reinforce,'--b');
% plot(rms_reinforce - std_reinforce,'--b');
% plot(rms_reinforce_baseline + std_reinforce_baseline,'--r');
% plot(rms_reinforce_baseline - std_reinforce_baseline,'--r');
% plot(rms_ac + std_ac,'--k');
% plot(rms_ac - std_ac,'--k');
legend("REINFORCE", "REINFORCE with Baseline", "Actor Critic")
box on;
