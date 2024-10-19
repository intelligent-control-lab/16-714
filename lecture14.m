% 16-714 Advanced Control for Robotics
% Script for Lecture 14: Least Square

%% Problem setup
% X ~ randn([0,0],eye(2))
% Y = AY * X, Z = AZ * X
% E[Y] = 0; E[Z] = 0;
% VAR[Y] = AY*AY'; VAR[Z] = AZ*AZ'; 
% COV[X,Y] = AY'; COV[X,Z] = AZ'; COV[Z,Y] = AZ*AY';
n = 100; % number of random samples
x = randn(2,n); % samples for X
AY = [1 0]; AZ = [1 1];
y = AY*x;
z = AZ*x;

%% Estimate X only using Y (Now we do not assume we know the exact sample values of x)
% E[X|Y] = E[X] + COV[X,Y]*inv(VAR[Y])[Y-E(Y)];
% VAR[X|Y] = VAR[X] - COV[X,Y]*inv(VAR[Y])*COV[Y,X];
x_hat_y = AY'* y / (AY*AY');

%% Estimate X using Y and Z
% E[Z|Y] = E[Z] + COV[Z,Y]*inv(VAR[Y])[Y-E(Y)];
% VAR[Z|Y] = VAR[Z] - COV[Z,Y]*inv(VAR[Y])*COV[Y,Z];
z_hat_y = AZ*AY' * y / (AY*AY');
z_var_y = AZ*AZ' - AZ*AY'*AY*AZ'/ (AY*AY');
z_tilde_y = z - z_hat_y;
% E[X|Y,Z] = E[X|Y] + E[tilde X|tilde Z]
% COV[tilde X, tilde Z] = COV(X,Z) - COV(X,Y)*COV(Y,Z)/VAR(Y);
x_hat_yz = x_hat_y + (AZ' - AY'*AY*AZ'/ (AY*AY')) * z_tilde_y / z_var_y;

%% Show Results
% Now let us check if we get a good estimation of x
display('Max estimation error for E[X|Y]:')
display(max(abs(x-x_hat_y)'));
display('Max estimation error for E[X|Y,Z]:')
display(max(abs(x-x_hat_yz)'));