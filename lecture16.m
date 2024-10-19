% 16-714 Advanced Control for Robotics
% Script for Lecture 16: EKF and UKF

%% Problem Setup
problem = 3;
switch problem
    case 1
        sys.x0 = [1;2];
        sys.X0 = 0.001*[1 0;0 2]; sys.x0hat = sys.x0 + sqrt(sys.X0)*randn(2,1);%A posteriori
        sys.W = 0.001;
        sys.V = 0.001;
        sys.Bw = 1;
        sys.fnominal = @(x,u) x+u;
        sys.f = @(x,u) sys.fnominal(x,u) + sys.Bw * sqrt(sys.W) * randn();
        sys.hnominal = @(x) [x(1)^2+x(2)^2;x(1)*x(2)];
        sys.h = @(x) sys.hnominal(x) + sqrt(sys.V) * randn();
        sys.A = @(x,u) eye(2);
        sys.C = @(x) [2*x(1) 2*x(2); x(2) x(1)];
    case 2
        sys.x0 = [1;2];
        sys.X0 = 0.01*[1 0;0 2]; sys.x0hat = sys.x0 + sqrt(sys.X0)*randn(2,1);%A posteriori
        sys.W = 0.001;
        sys.V = 0.001;
        sys.Bw = 1;
        sys.fnominal = @(x,u) x+u;
        sys.f = @(x,u) sys.fnominal(x,u) + sys.Bw * sqrt(sys.W) * randn();
        sys.hnominal = @(x) [x(1)^2+x(2)^2;x(1)*x(2)];
        sys.h = @(x) sys.hnominal(x) + sqrt(sys.V) * randn();
        sys.A = @(x,u) eye(2);
        sys.C = @(x) [2*x(1) 2*x(2); x(2) x(1)];
    case 3
        sys.x0 = [0.1;-1];
        sys.X0 = 0.01*[8 2;2 3]; sys.x0hat = sys.x0 + sqrt(sys.X0)*randn(2,1);%A posteriori
        sys.W = 0.001;
        sys.V = 0.001;
        sys.Bw = 1;
        sys.fnominal = @(x,u) x+u;
        sys.f = @(x,u) sys.fnominal(x,u) + sys.Bw * sqrt(sys.W) * randn();
        sys.hnominal = @(x) [x(1)^2+x(2)^2;x(1)*x(2)];
        sys.h = @(x) sys.hnominal(x) + sqrt(sys.V) * randn();
        sys.A = @(x,u) eye(2);
        sys.C = @(x) [2*x(1) 2*x(2); x(2) x(1)];
    case 4
        sys.x0 = [4;1];
        sys.X0 = 0.01*[8 2;2 3]; sys.x0hat = sys.x0 + sqrt(sys.X0)*randn(2,1);%A posteriori
        sys.W = 0.001;
        sys.V = 0.001;
        sys.Bw = 1;
        sys.fnominal = @(x,u) x+u;
        sys.f = @(x,u) sys.fnominal(x,u) + sys.Bw * sqrt(sys.W) * randn();
        sys.hnominal = @(x) [x(1)^4*x(2);x(1)+x(2)^5];
        sys.h = @(x) sys.hnominal(x) + sqrt(sys.V) * randn();
        sys.A = @(x,u) eye(2);
        sys.C = @(x) [4*x(1)^3*x(2) x(1)^4; 1 5*x(2)^4];
end

sys.N = 10;
sys.dt = 1;

sys.simmode = 'Direct'; % The dynamics are directly specified in discrete time
% Controller
u = @(x,t) [0;0];
% Simulate
[~, xlist, ulist, ylist] = roll_out(sys, u, sys.x0, 'DT', sys.N, sys.dt, 'Direct');
% Plot
figure(1);clf;
subplot(211);hold on;
plot(xlist(1,:),'k')
subplot(212);hold on;
plot(xlist(2,:),'k')

%% Estimate States
ekf = estimator('EKF',sys);
ukf = estimator('UKF',sys);

%% Filtering
for k = 1:sys.N
    ekf.xhat(:,k+1) = ekf.update_x(ekf.xhat(:,k),ulist(:,k),ylist(:,k+1),ekf.Z(:,:,k));
    ekf.Z(:,:,k+1) = ekf.update_Z(ekf.Z(:,:,k), ekf.xhat(:,k),ulist(:,k));

    ukf.xhat(:,k+1) = ukf.update_x(ukf.xhat(:,k),ulist(:,k),ylist(:,k+1),ukf.Z(:,:,k));
    ukf.Z(:,:,k+1) = ukf.update_Z(ukf.Z(:,:,k), ukf.xhat(:,k), ulist(:,k));
end
display('Estimation Error (EKF):')
display(norm(ekf.xhat-xlist));
display('Estimation Error (UKF):')
display(norm(ukf.xhat-xlist));

% Plot
figure(1);
subplot(211);hold on;
plot(ekf.xhat(1,:),'--b')
plot(ukf.xhat(1,:),'b')
legend("Ground truth x", "EKF", "UKF")
subplot(212);hold on;
plot(ekf.xhat(2,:),'--b')
plot(ukf.xhat(2,:),'b')
legend("Ground truth x", "EKF", "UKF")

figure(2);clf;hold on;
string = [];
for k = 1:2
    for l = 1:2
        plot(permute(ekf.Z(k,l,:),[3,2,1]));
        plot(permute(ukf.Z(k,l,:),[3,2,1]));
        string = [string "EKF Z_{"+num2str(k)+num2str(l)+"}" "UKF Z_{"+num2str(k)+num2str(l)+"}"];
    end
end
legend(string)
