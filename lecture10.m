% 16-714 Advanced Control for Robotics
% Script for Lecture 10: MPC feasibility

%% Specify problem
problem = 2;
switch problem
    case 1
        sys.dt = 0.5; % Sampling time
        sys.name = 'single_integrator';
        [sys.A,sys.B] = getAB('single_integrator',1,'DT',sys.dt);
        sys.Q = 1;
        sys.R = 10;
        sys.S = 100;
        sys.x0 = [10];
        sys.N = 10;
        sys.umax = 1;
        sys.umin = -1;
        sys.xmax = 10;
        sys.xmin = -1;
    case 2
        sys.dt = 0.5; % Sampling time
        sys.name = 'double_integrator';
        [sys.A,sys.B] = getAB('double_integrator',2,'DT',sys.dt);
        sys.Q = [1 0; 0 0];
        sys.R = 1;
        sys.S = 10*eye(2);
        sys.x0 = [10; 0];
        sys.N = 10;
        sys.umax = 1;
        sys.umin = -1;
        sys.xmax = [10; 5];
        sys.xmin = [-10; -2];
end

%% Controllers
cMPCn = synthesis('nMPC', sys);
cMPCu = synthesis('uMPC', sys);
cMPCq = synthesis('qMPC', sys);

%% Visualiza trajectory
simhorizon = 3*sys.N;
[~, x_mpcn, u_mpcn] = roll_out(sys.name, cMPCn.u, sys.x0, 'DT', simhorizon, sys.dt);
[~, x_mpcu, u_mpcu] = roll_out(sys.name, cMPCu.u, sys.x0, 'DT', simhorizon, sys.dt);
[~, x_mpcq, u_mpcq] = roll_out(sys.name, cMPCq.u, sys.x0, 'DT', simhorizon, sys.dt);

figure(1); clf; hold on
subplot(length(sys.xmax),1,1);
for i = 1:length(sys.xmax)
    subplot(length(sys.xmax),1,i);hold on;
    plot(0:simhorizon, x_mpcn(i,:));
    plot(0:simhorizon, x_mpcu(i,:));
    plot(0:simhorizon, x_mpcq(i,:));
    plot(0:simhorizon, zeros(1,simhorizon+1),".k")
    plot(0:simhorizon, sys.xmax(i).*ones(1,simhorizon+1),"--k")
    plot(0:simhorizon, sys.xmin(i).*ones(1,simhorizon+1),"--k")
    box on
    legend("Unconstrained", "Control Constrained", "State and Control Constrained")
end

figure(2); clf; hold on
plot(0:simhorizon-1, u_mpcn);
plot(0:simhorizon-1, u_mpcu);
plot(0:simhorizon-1, u_mpcq);
plot(0:simhorizon-1, sys.umax*ones(1,simhorizon),"--k")
plot(0:simhorizon-1, sys.umin*ones(1,simhorizon),"--k")
box on
legend("Unconstrained", "Control Constrained", "State and Control Constrained")

%%
feasibile = visualize_mpc(sys, cMPCq, simhorizon);