function p = par(name)
% par  Parameter registry for f(x,u,name)
%      Usage: p = par('bicycle_dynamic');

switch name

    case 'dubins'
        p.v = 1.0;                 % constant forward speed [m/s]

    case 'bicycle_dynamic'
        % Planar bicycle (single-track) with linear tires
        p.m      = 1500;           % mass [kg]
        p.Iz     = 2250;           % yaw inertia [kg m^2]
        p.lf     = 1.2;            % CG -> front axle [m]
        p.lr     = 1.6;            % CG -> rear axle [m]
        p.Cf     = 8e4;            % front cornering stiffness [N/rad]
        p.Cr     = 9e4;            % rear  cornering stiffness [N/rad]
        p.beta   = 0.5;            % drive split (0=RWD, 1=FWD, 0.5=AWD)
        p.eps_vx = 0.1;            % low-speed guard [m/s]
        p.viz_wheel_len = 0.8;   % purely for drawing
        p.viz_wheel_w   = 0.25;
        p.k_cte = 2.0; p.k_v = 1.2; p.k_speed = 0.8;
        p.v_ref_max = 8; p.Fmax = 6000; p.delta_max = 0.6;

    case 'arm_first_order'
        p.endEffector = "EndEffector_Link";

    case 'arm_second_order'
        p.endEffector = "EndEffector_Link";
        p.eps = 0.1;               % finite-difference step for Jdot

    otherwise
        p = struct();              % models that need no constants
end
end