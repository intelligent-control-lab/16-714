% 16-714 Advanced Control for Robotics
% Script for Lecture 22: Compare ILC with PG  (time-varying via model only)

%% Problem Setup
sys.A = 0.5; sys.B = 1;

% Pre-generate seeds (for reproducibility)
n_ep = 20;
seeds = randi(1e6, n_ep, 1);

scenario = 'non-repetitive'; % 'non-repetitive', 'state-dependent'
switch scenario
    case 'repetitive'
        sys.f_dt = @(x,u,t) sys.A*x + sys.B*u + sin(0.1*pi*t);
    case 'non-repetitive'
        sys.f_dt = @(x,u,t) sys.A*x + sys.B*u + randn();
    case 'state-dependent'
        sys.f_dt = @(x,u,t) sys.A*x + sys.B*u + sin(0.1*pi*t)*x;
end
sys.Q = 1; sys.R = 1;
sys.h = @(x,u) (x'*sys.Q*x + u'*sys.R*u)/2;
sys.kmax = 50; sys.dt = 1; sys.x0 = 0.1;
sys.tcondition = @(x,t) t > sys.kmax * sys.dt;

% roll_out configuration (no special timevarying opts needed)
simDT = struct('type','DT','K',sys.kmax,'dt',sys.dt, ...
               'integrator','Direct', ...
               'stop', @(xk,tk,log) sys.tcondition(xk,tk));
optsDT = struct('integrator','Direct'); % keep minimal opts if your roll_out expects it

sysnominal = sys;
sysnominal.f = @(x,u,t) sys.A * x + sys.B * u;

%% LQR solution
c = synthesis('LQR',sys);
[~, x_fb]  = roll_out(sys,        c.u, sys.x0, simDT, optsDT);
[~, x_lqr] = roll_out(sysnominal, c.u, sys.x0, simDT, optsDT);

%% Learning parameters
sys.W0 = 3; sys.delta = 1;
sys.policy.name       = 'Gaussian';
sys.policy.init.mu    = -0.1*zeros(size(sys.B')); % scalar
sys.policy.init.sigma = 0.1*eye(size(sys.B'));    % scalar

alpha = 1.0e-4;         % policy step
alpha_baseline = 1.0e-2; % value step
alphas = struct('p',alpha,'v',alpha_baseline);

%% REINFORCE
reinforce = PGagent(sys);
rms_reinforce = zeros(1,n_ep);
std_reinforce = zeros(1,n_ep);
figure(1); clf; subplot(311); hold on;
for episode = 1:n_ep
    rng(seeds(episode));
    x_list = reinforce.update('REINFORCE', alpha, sys); % uses roll_out on time-varying sys
    plot(0:length(x_list)-1, x_list, 'color', [0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep])
    rms_reinforce(episode) = norm(reinforce.theta.mu + c.K);
    std_reinforce(episode) = reinforce.theta.sigma;
end
plot(0:length(x_lqr)-1, x_lqr, '--k'); title("REINFORCE"); box on;

%% REINFORCE with Baseline
reinforce_bl = PGagent(sys);
rms_reinforce_bl = zeros(1,n_ep);
std_reinforce_bl = zeros(1,n_ep);
subplot(312); hold on;
for episode = 1:n_ep
    rng(seeds(episode));
    x_list = reinforce_bl.update('REINFORCEbl', alphas, sys);
    plot(0:length(x_list)-1, x_list, 'color', [0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep])
    rms_reinforce_bl(episode) = norm(reinforce_bl.theta.mu + c.K);
    std_reinforce_bl(episode) = reinforce_bl.theta.sigma;
end
plot(0:length(x_lqr)-1, x_lqr, '--k'); title("REINFORCE with Baseline"); box on;

%% Actor-Critic
actor_critic = PGagent(sys);
rms_ac = zeros(1,n_ep);
std_ac = zeros(1,n_ep);
subplot(313); hold on;
for episode = 1:n_ep
    rng(seeds(episode));
    x_list = actor_critic.update('ActorCritic', alphas, sys);
    plot(0:length(x_list)-1, x_list, 'color', [0.2+0.8*episode/n_ep, 1-0.8*episode/n_ep, 1-0.9*episode/n_ep])
    rms_ac(episode) = norm(actor_critic.theta.mu + c.K);
    std_ac(episode) = actor_critic.theta.sigma;
end
plot(0:length(x_lqr)-1, x_lqr, '--k'); title("Actor Critic"); box on;

%% Theta convergence
figure(2); clf; hold on
plot(rms_reinforce,'b');
plot(rms_reinforce_bl,'r');
plot(rms_ac,'k');
legend("REINFORCE","REINFORCE+BL","Actor-Critic"); box on;

%% ILC setup and loop
sys.K = 0; sys.N = sys.kmax; sys.C = 1;
cILCt = synthesis('ILCt', sys);

error = zeros(n_ep,sys.kmax+1);
ff    = zeros(n_ep, sys.kmax);
e     = zeros(1,n_ep);
legstr = ["Reference", "Iter 0"];

figure(3); clf; hold on;
for i = 1:n_ep
    rng(seeds(i));
    error(i,:) = x_lqr(1:sys.kmax+1) - x_fb(1:sys.kmax+1);

    % feedforward update
    if i > 1
        ff(i,1:sys.kmax) = cILCt.u(error(i,:), ff(i-1,1:sys.kmax));
    else
        ff(i,1:sys.kmax) = cILCt.u(error(i,:), zeros(1,sys.kmax));
    end

    % feedforward controller
    uff = @(x,t) ff(i, t/sys.dt + 1);

    % simulate on time-varying sys (no special opts)
    [tlist, x_fb, u_fb] = roll_out(sys, uff, sys.x0, simDT, optsDT);

    % plots
    subplot(211); hold on
    plot(0:sys.kmax,   x_fb(1,1:sys.kmax+1), 'color', [i/n_ep,0,0]); xlim([0,sys.kmax])
    subplot(212); hold on
    plot(0:sys.kmax-1, u_fb(1,1:sys.kmax),   'color', [i/n_ep,0,0])

    e(i) = norm(error(i,:));
    legstr = [legstr, "Iter " + num2str(i)];
end
subplot(211); legend(legstr); title('Output'); subplot(212); title('Input');

figure(4); clf; plot(0:n_ep-1, e); title('Norm Error in Different ILC Iterations');
