function xnew = step(x, u, dt, dyn, mode)
% STEP  Advance one step/sample for CT/DT models given resolved dyn.
%
%   mode: 'Euler' | 'ZOH' | 'RK4' | 'Direct'
%
% Requirements:
%   - 'Euler', 'ZOH', 'RK4' require dyn.f_ct (continuous-time RHS)
%   - 'Direct' requires dyn.f_dt (discrete-time map)

    arguments
        x  double {mustBeVector}
        u  double {mustBeVector}
        dt double {mustBeNonnegative}
        dyn struct
        mode (1,:) char
    end
    x = x(:); u = u(:);
    md = lower(mode);

    switch md
        case 'euler'
            assert(~isempty(dyn.f_ct), 'Euler requires CT RHS dyn.f_ct.');
            xnew = x + dt * col(dyn.f_ct(x,u));

        case 'zoh'
            assert(~isempty(dyn.f_ct), 'ZOH requires CT RHS dyn.f_ct.');
            if dt == 0, xnew = x; return; end
            odefn = @(t,y) col(dyn.f_ct(y,u));  % u held constant over [0,dt]
            opts  = odeset('RelTol',1e-6,'AbsTol',1e-8);
            [~, Y] = ode45(odefn, [0, dt], x, opts);
            xnew   = Y(end,:).';

        case 'rk4'
            assert(~isempty(dyn.f_ct), 'RK4 requires CT RHS dyn.f_ct.');
            f  = @(x_) col(dyn.f_ct(x_, u));   % hold u constant over the step
            k1 = f(x);
            k2 = f(x + 0.5*dt*k1);
            k3 = f(x + 0.5*dt*k2);
            k4 = f(x + dt*k3);
            xnew = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

        case 'direct'
            assert(~isempty(dyn.f_dt), 'Direct mode requires DT map dyn.f_dt.');
            xnew = col(dyn.f_dt(x,u));

        otherwise
            error('Unknown mode "%s". Use Euler, ZOH, RK4, or Direct.', mode);
    end
end

function y = col(z), y = z(:); end
