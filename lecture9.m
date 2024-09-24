% 16-714 Advanced Control for Robotics
% Script for Lecture 9: MPC

%% Specify problem
problem = 1;
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
    case 2
        sys.dt = 0.5; % Sampling time
        sys.name = 'double_integrator';
        [sys.A,sys.B] = getAB('double_integrator',2,'DT',sys.dt);
        sys.Q = [1 0; 0 0];
        sys.R = 1;
        sys.S = 10*eye(2);
        sys.x0 = [10; 0];
        sys.N = 10;
end

%% Controllers
cLQR = synthesis('LQR', sys);
cLQRn = synthesis('LQRn', sys);
cMPC = synthesis('MPC', sys);

%% Compare the P matrices
figure; hold on;
string = [];
for i = 1:size(cMPC.P,1)
    for j = 1:size(cMPC.P,2)
        plot(0:sys.N,permute(cMPC.P(i,j,1:sys.N+1),[3,1,2]),"LineStyle","-","Marker","*")
        string = [string, "P_{"+num2str(i)+num2str(j)+"}"];
    end
end
for i = 1:size(cLQR.P,1)
    for j = 1:size(cLQR.P,1)
        plot(0:sys.N,cLQR.P(i,j)*ones(1,sys.N+1),"LineStyle","--")
    end
end
box on
legend(string)

%% Compare the gain K
figure; hold on;
string = [];
for i = 1:size(cMPC.K,1)
    for j = 1:size(cMPC.K,2)
        plot(0:sys.N-1,permute(cMPC.K(i,j,1:sys.N),[3,1,2]),"LineStyle","-","Marker","*")
        string = [string, "K_{"+num2str(i)+num2str(j)+"}"];
    end
end
for i = 1:size(cLQR.K,1)
    for j = 1:size(cLQR.K,1)
        plot(0:sys.N-1,cLQR.K(i,j)*ones(1,sys.N),"LineStyle","--")
    end
end
box on
legend(string)
%% Visualiza trajectory

[~, x_lqr_f, ~] = roll_out(sys.name, cLQRn.u, sys.x0, 'DT', sys.N, sys.dt);
[~, x_lqr_inf, ~] = roll_out(sys.name, cLQR.u, sys.x0, 'DT', sys.N, sys.dt);
[~, x_lqr_mp, ~] = roll_out(sys.name, cMPC.u, sys.x0, 'DT', sys.N, sys.dt);

figure; hold on
plot(0:sys.N, x_lqr_f(1,:));
plot(0:sys.N, x_lqr_inf(1,:));
plot(0:sys.N, x_lqr_mp(1,:));
plot(0:sys.N, zeros(1,sys.N+1),"--k")
box on
legend("Finite LQR", "Infinite LQR", "MPC")

%% MPC with Different Preview Horizon
x_lqr_mps = zeros(length(sys.x0), 2*sys.N+1, sys.N);
N = sys.N;
% Set up N different controllers
for i = 1:N
    sys.N = i;
    c{i} = synthesis('MPC', sys);
end
% Simulate for 2*N horizon
[~, x_lqr_inf, ~] = roll_out(sys.name, cLQR.u, sys.x0, 'DT', 2*sys.N, sys.dt);
for i = 1:N
    [~, x_lqr_mps(:,:,i), ~] = roll_out(sys.name, c{i}.u, sys.x0, 'DT', 2*sys.N, sys.dt);
end
% Visualize
figure; hold on
plot(0:N, x_lqr_f(1,:));
plot(0:2*N, x_lqr_inf(1,:));
string = ["Finite LQR", "Infinite LQR"];
for i = 1:N
    plot(0:2*N, x_lqr_mps(1,:,i), "color", [i/N, 1-i/N, i/N]);
    string = [string, "MPC-"+num2str(i)];
end
plot(0:2*N, zeros(1,2*N+1),"--k")
box on
legend(string)
%% Predicted Output at Each MPC Step
% Essentially at every time step, a finite horizon LQR is solved.
x_lqr_mp_steps = zeros(length(sys.x0), sys.N+1, sys.N+1);
[~, x_lqr_mp_steps(:,:,1), ~] = roll_out(sys.name, cLQRn.u, sys.x0, 'DT', sys.N, sys.dt);
for t = 1:sys.N
    [~, x_lqr_mp_steps(:,:,t+1), ~] = roll_out(sys.name, cLQRn.u, x_lqr_mp_steps(:,2,t), 'DT', sys.N, sys.dt);
end

figure; hold on
plot(0:sys.N, x_lqr_f(1,:));
plot(0:2*sys.N, x_lqr_inf(1,:));
string = ["Finite LQR", "Infinite LQR"];
for t = 1:sys.N
    plot(t-1:t+sys.N-1, x_lqr_mp_steps(1,:,t), "color", [1-t/N, t/N, 1-t/N]);
    string = [string, "MPC Time "+num2str(t-1)];
end
box on
legend(string)
