% 16-714 Advanced Control for Robotics
% Script for Lecture 22: Compare ILC with PG

%% Problem Setup
sys.A = 0.5; sys.B = 1;

scenario = 'repetitive';%'non-repetitive', 'state-dependent'
switch scenario
    case 'repetitive'
        sys.f = @(x,u,t) sys.A*x + sys.B*u + sin(0.1*pi*t);
    case 'non-repetitive'
        d = randn;
        sys.f = @(x,u,t) sys.A*x + sys.B*u + d;
    case 'state-dependent'
        sys.f = @(x,u,t) sys.A*x + sys.B*u + sin(0.1*pi*t)*x;
end
sys.Q = 1; sys.R = 1;
sys.h = @(x,u) (x'*sys.Q*x + u'*sys.R*u)/2;
sys.kmax = 50; sys.dt = 1; sys.x0 = 0.1;
sys.tcondition = @(x,t) t > sys.kmax * sys.dt;

sysnominal = sys;
sysnominal.f = @(x,u) sys.A * x + sys.B * u;

%% LQR solution
c = synthesis('LQR',sys);
[~, x_fb] = roll_out(sys, c.u, sys.x0, 'DT', sys.tcondition, sys.dt, 'Direct');
[~, x_lqr] = roll_out(sysnominal, c.u, sys.x0, 'DT', sys.tcondition, sys.dt, 'Direct');


%% Learning paramters
sys.W0 = 3; sys.delta = 1;

sys.policy.name = 'Gaussian';
sys.policy.init.mu = -0.1*zeros(size(sys.B'));
sys.policy.init.sigma = 0.1*eye(size(sys.B'));

n_ep = 20;
alpha = 0.0001; %policy
alpha_baseline = 0.01; %value
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
legend("REINFORCE", "REINFORCE with Baseline", "Actor Critic")
box on;

%%
sys.K = 0; sys.N = sys.kmax; sys.C = 1;
cILCt = synthesis('ILCt', sys);

%% ILC
error = zeros(n_ep,sys.kmax+1);
ff = zeros(n_ep, sys.kmax);
e = zeros(1,n_ep);
string = ["Reference", "Iter 0"];
figure(3);clf;hold on;
for i = 1:n_ep
    error(i,:) = x_lqr(1:sys.kmax+1) - x_fb(1:sys.kmax+1);
    % filter the error signal
    if i > 1
        ff(i,1:sys.kmax) = cILCt.u(error(i,:), ff(i-1,1:sys.kmax));
    else
        ff(i,1:sys.kmax) = cILCt.u(error(i,:), zeros(1,sys.kmax));
    end
    % Set the feedforward controller
    uff = @(x,t) ff(i,t/sys.dt+1);
    % Simulate
    [tlist, x_fb, u_fb] = roll_out(sys, uff, sys.x0, 'DT', sys.tcondition, sys.dt, 'Direct');
    % plot
    subplot(211); hold on
    plot(0:sys.kmax, x_fb(1,1:sys.kmax+1), 'color',[i/n_ep,0,0])
    xlim([0,sys.kmax])
    subplot(212); hold on
    plot(0:sys.kmax-1, u_fb(1,1:sys.kmax),'color',[i/n_ep,0,0])
    % Compute error norm
    e(i) = norm(error(i,:));
    string = [string,"Iter "+num2str(i)];
end
subplot(211);
legend(string);title('Output')
subplot(212);
title('Input')
figure(4);clf;
plot(0:n_ep-1, e);
title('Norm Error in Different ILC Iterations')
