% 16-714 Advanced Control for Robotics
% Script for Lecture 13: ICL Time Domain

%% Problem Specification
sys.dt = 1; % Sampling time
sys.name = 'double_integrator';
[sys.A,sys.B] = getAB('double_integrator',2,'DT',sys.dt,'Euler');
sys.a = [1 -2 1];
sys.b = [sys.dt];
sys.C = [1 0]; % only measuring the position
sys.x0 = [0;0];
sys.N = 100;

% Reference signal
ref = @(t) (sin(0.2*t));
simhorizon = sys.N;
niter = 10;
sigma = 0.02; % control disturbance; 0 if no disturbance

%% Specify feedback control
kp = 0.1; kd = 0.5;
sys.K = [kp kd]; % For time domain use
sys.c = [-kd/sys.dt, kd/sys.dt - kp]; % for frequency domain use
ufb = @(x,t) kp*(ref(t) - sys.C*x) + kd*((ref(t+sys.dt)-ref(t))/sys.dt - x(2)) + sigma * randn; % proportional error feedback

% First trajectory (Iteration 0)
[tlist, x_fb, u_fb] = roll_out(sys.name, ufb, sys.x0, 'DT', simhorizon, sys.dt, 'Euler');

% Plot
figure(1);clf;subplot(211); hold on
plot(0:simhorizon,ref(tlist),'--')
plot(0:simhorizon,x_fb(1,:),'k')
subplot(212)
plot(0:simhorizon-1,u_fb)

%% Specify the learning gain
cILCt = synthesis('ILCt', sys);

%% ILC
error = zeros(niter,simhorizon+1);
ff = zeros(niter, simhorizon);
e = zeros(1,niter);
string = ["Reference", "Iter 0"];
for i = 1:niter
    error(i,:) = ref(tlist) - x_fb(1,:);
    % filter the error signal
    if i > 1
        ff(i,1:simhorizon) = cILCt.u(error(i,:), ff(i-1,1:simhorizon));
    else
        ff(i,1:simhorizon) = cILCt.u(error(i,:), zeros(1,simhorizon));
    end
    % Set the feedforward controller
    uff = @(x,t) ff(i,t/sys.dt+1) + ufb(x,t);
    % Simulate
    [tlist, x_fb, u_fb] = roll_out(sys.name, uff, sys.x0, 'DT', simhorizon, sys.dt, 'Euler');
    % plot
    subplot(211); hold on
    plot(0:simhorizon, x_fb(1,:), 'color',[i/niter,0,0])
    xlim([0,simhorizon])
    subplot(212); hold on
    plot(0:simhorizon-1, u_fb,'color',[i/niter,0,0])
    % Compute error norm
    e(i) = norm(error(i,:));
    string = [string,"Iter "+num2str(i)];
end
subplot(211);
legend(string);title('Output')
subplot(212);
title('Input')
figure(2);clf;
plot(0:niter-1, e);
title('Norm Error in Different ILC Iterations')

%% Compare two different ILC versions by turning the learning transfer 
% function from the frequency domain ILC into a matrix
cILCw = synthesis('ILCw', sys);
L_matrix_form = zeros(simhorizon, simhorizon+1);
l = length(cILCw.L.b);
for i = 1:simhorizon
    j = 0;
    while j < l && i + j <= simhorizon + 1
        L_matrix_form(i,i+j) = cILCw.L.b(l-j);
        j = j+1;
    end
end

% To verify the L matrix is constructed correctly
% ff(1,1:simhorizon)' -  L_matrix_form * error(1,:)'

figure(10);subplot(121)
image((cILCt.L+1.5)*100)
colorbar('XTickLabel',{'-1','0','1'},'XTick',[50,150,250])
title('L matrix in the time domain ILC')
subplot(122)
image((L_matrix_form+1.5)*100)
colorbar('XTickLabel',{'-1','0','1'},'XTick',[50,150,250])
title('L matrix in the frequency domain ILC')