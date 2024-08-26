%% Library of dynamic equations; All in Continuous Time
% Currently the following dynamics are supported
% Single Integrator
% Double Integrator
% Dubins Car
% Unicycle 3 state
% Unicycle 4 STATE
% First order arm J2C
% Second order arm J2C

function xdot = f(x, u, name)
switch name
    case 'single_integrator'
        xdot = u;
    case 'double_integrator'
        n = length(x)/2;
        xdot(1:n,1) = x(n+1:2*n);
        xdot(n+1:2*n,1) = u;
    case 'dubins'
        v = 1; % velocity is constant
        xdot(1,1) = v*cos(x(3)); % x position
        xdot(2,1) = v*sin(x(3)); % y position
        xdot(3,1) = u; % theta
    case 'unicycle3'
        xdot(1,1) = u(1)*cos(x(3)); % x position
        xdot(2,1) = u(1)*sin(x(3)); % y position
        xdot(3,1) = u(2); % theta
    case 'unicycle4'
        xdot(1,1) = x(3)*cos(x(4)); % x position
        xdot(2,1) = x(3)*sin(x(4)); % y position
        xdot(3,1) = u(1); % velocity
        xdot(4,1) = u(2); % theta
    case 'arm_first_order' % position only model
        global robot; ik = inverseKinematics("RigidBodyTree",robot);
        endEffector = "EndEffector_Link";
        theta = ik(endEffector,trvec2tform(x'),[1 1 1 1 1 1],homeConfiguration(robot));
        geoJacob = geometricJacobian(robot, theta, endEffector);
        xdot = geoJacob(1:length(x),:) * u;
    case 'arm_second_order' % position only model
        global robot; ik = inverseKinematics("RigidBodyTree",robot);
        endEffector = "EndEffector_Link";
        theta = ik(endEffector,trvec2tform(x'),[1 1 1 1 1 1],homeConfiguration(robot));
        geoJacob = geometricJacobian(robot, theta, endEffector);
        % a heuristic implementation for JacobDot
        eps = 0.1;
        geoJacobDot = (geometricJacobian(robot, theta+eps, endEffector)-geometricJacobian(robot, theta-eps, endEffector))/(2*eps);
        n = length(x)/2;
        xdot(1:n,1) = x(n+1:2*n);
        xdot(n+1:2*n,1) = geoJacob(1:n,:) * u + geoJacobDot(1:n,:);
    otherwise
        xdot = zeros(size(x));
end
end