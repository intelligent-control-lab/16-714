%% Library of dynamic equations; All in Continuous Time

function xdot = f(x, u, name)
%    Supported names:
%      'single_integrator', 'double_integrator',
%      'dubins', 'unicycle3', 'unicycle4',
%      'arm_first_order', 'arm_second_order',
%      'bicycle_dynamic' (planar single-track with linear tires)

xdot = zeros(size(x));

switch name

    % --------------------- simple point-mass models ---------------------
    case 'single_integrator'
        xdot = u;

    case 'double_integrator'
        n = length(x)/2;
        xdot(1:n,1)       = x(n+1:2*n);
        xdot(n+1:2*n,1)   = u;

    case 'dubins'
        p = par('dubins');              % constants in par.m
        v = p.v;
        xdot(1,1) = v * cos(x(3));      % x
        xdot(2,1) = v * sin(x(3));      % y
        xdot(3,1) = u(1);               % theta rate (turn input)

    case 'unicycle3'
        % x = [x; y; theta], u = [v; w]
        xdot(1,1) = u(1) * cos(x(3));
        xdot(2,1) = u(1) * sin(x(3));
        xdot(3,1) = u(2);

    case 'unicycle4'
        % x = [x; y; v; theta], u = [a; w]
        xdot(1,1) = x(3) * cos(x(4));
        xdot(2,1) = x(3) * sin(x(4));
        xdot(3,1) = u(1);
        xdot(4,1) = u(2);

    % ------------------------- arm ---------------------------
    case 'arm_first_order' % position-only model
        p = par('arm_first_order');
        global robot;
        ik = inverseKinematics("RigidBodyTree", robot);
        theta = ik(p.endEffector, trvec2tform(x'), [1 1 1 1 1 1], homeConfiguration(robot));
        geoJacob = geometricJacobian(robot, theta, p.endEffector);
        xdot = geoJacob(1:length(x),:) * u;

    case 'arm_second_order' % position + velocity states
        p = par('arm_second_order');
        global robot;
        ik = inverseKinematics("RigidBodyTree", robot);
        theta = ik(p.endEffector, trvec2tform(x'), [1 1 1 1 1 1], homeConfiguration(robot));
        geoJacob = geometricJacobian(robot, theta, p.endEffector);
        % heuristic Jacobian time-derivative
        epsh = p.eps;
        Jp = geometricJacobian(robot, theta + epsh, p.endEffector);
        Jm = geometricJacobian(robot, theta - epsh, p.endEffector);
        Jdot = (Jp - Jm) / (2*epsh);

        n = length(x)/2;
        xdot(1:n,1)       = x(n+1:2*n);
        xdot(n+1:2*n,1)   = geoJacob(1:n,:) * u + Jdot(1:n,:);

    % -------------------- bicycle dynamic (planar) ----------------------
    case {'bicycle_dynamic','bicycle'}  % alias
        % State:  x = [X; Y; v_x; v_y; r; psi]
        % Input:  u = [F_x; delta_f]
        p = par('bicycle_dynamic');

        X   = x(1); %#ok<NASGU>
        Y   = x(2); %#ok<NASGU>
        vx  = x(3); vy = x(4);
        r   = x(5); psi = x(6);
        Fx  = u(1); delta = u(2);

        % robust longitudinal speed for slip-angles
        vx_eff = vx;
        if abs(vx_eff) < p.eps_vx
            sgn = sign(vx_eff); if sgn == 0, sgn = 1; end
            vx_eff = sgn * p.eps_vx;
        end

        % slip angles (exact geometry)
        alpha_f = atan2(vy + p.lf*r, vx_eff) - delta;
        alpha_r = atan2(vy - p.lr*r, vx_eff);

        % lateral tire forces (linear)
        Fyf = -p.Cf * alpha_f;
        Fyr = -p.Cr * alpha_r;

        % split longitudinal force (FWD/RWD/AWD)
        Fxf = p.beta * Fx;
        Fxr = (1 - p.beta) * Fx;

        % body-frame dynamics
        dvx = ( Fxf*cos(delta) - Fyf*sin(delta) + Fxr )/p.m + r*vy;
        dvy = ( Fxf*sin(delta) + Fyf*cos(delta) + Fyr )/p.m - r*vx;
        dr  = ( p.lf*(Fxf*sin(delta) + Fyf*cos(delta)) - p.lr*Fyr )/p.Iz;

        % world-frame kinematics
        dpsi = r;
        dX = vx*cos(psi) - vy*sin(psi);
        dY = vx*sin(psi) + vy*cos(psi);

        xdot = [dX; dY; dvx; dvy; dr; dpsi];

    otherwise
        % unknown name -> zeros already set
end
end
