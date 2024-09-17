% 16-714 Advanced Control for Robotics
% Script for Lecture 7: LQR

%% Problem Parameters
dim = 4; dyn = 'double_integrator';
x0 = [10;10;1;5];
Tmax = 10;
dt = 0.1;
Q = [eye(dim/2) zeros(dim/2); zeros(dim/2) zeros(dim/2)]; 
R = eye(dim/2);

%% Discrete time LQR
[A,B] = getAB(dyn,dim,'DT',dt);

% Compute the Riccati Equation recursively
N = 100;
P = zeros(dim, dim, N); P(:,:,N) = eye(dim);
for i = N-1:-1:1
    P(:,:,i) = Q + A' * P(:,:,i+1) * A - A' * P(:,:,i+1) * B * inv(R + B' * P(:,:,i+1) * B) * B' * P(:,:,i+1) * A;
end
% Plot
figure(1); hold on; string = [];
for k = 1:dim
    for l = 1:dim
        plot(permute(P(k,l,:),[3,2,1]));
        string = [string "p_{"+num2str(k)+num2str(l)+"}"];
    end
end
legend(string)
% Compare with the dlqr solution
% K == control_gain
% Pd == P(:,:,1)
control_gain = inv(R + B' * P(:,:,1) * B) * B' * P(:,:,1) * A;
[K, Pd] = dlqr(A, B, Q, R);

%% Continuous Time LQR
[Ac,Bc] = getAB(dyn,dim,'CT');

[Kc, Pc] = lqr(Ac, Bc, Q, R);

u = @(x,t) -Kc*x;
[tlistc, xlistc, ulistc] = roll_out(dyn, u, x0, 'CT', Tmax);

figure(2);clf;hold on;
plot(xlistc(1,:), xlistc(2,:),'k')
%% Compare the Continuous Time and Discrete Time LQR
u = @(x,t) -K*x;
[tlist, xlist, ulist] = roll_out(dyn, u, x0, 'CT', Tmax, dt);

plot(xlist(1,:), xlist(2,:),'r')
legend("Continuous Time LQR Trajectory","Discrete Time LQR Trajectory")

%% Alternative implementation
sys.dt = dt;
sys.A = A;
sys.B = B;
sys.Q = Q;
sys.R = R;
c = synthesis('LQR',sys);