% 16-714 Advanced Control for Robotics
% Script for the investiment problem with state constraints

r = [0,1,0,1,0,1];
p = [1,2,3,3,2,1];

% Backward propagation
lambda = zeros(2,7);
lambda(:,7) = [-1;0];
u = zeros(1,6);
for i = 6:-1:1
    if [-p(i),1]*lambda(:,i+1) > 0
        lambda(:,i) = [1,0;p(i)+r(i),0] * lambda(:,i+1);
        u(i) = -1;
    else
        lambda(:,i) = [0,1/p(i);0,r(i)/p(i)+1] * lambda(:,i+1);
        u(i) = 1;
    end
end

% Forward simualtion
x = zeros(2,7);
x(:,1) = [10,0];
for i = 1:6
    if u(i) == 1
        u(i) = (x(1,i)+r(i)*x(2,i))/p(i);
    else
        u(i) = -x(2,i);
    end
    x(:,i+1) = [1,r(i);0,1]*x(:,i) + [-p(i);1]*u(i);
end