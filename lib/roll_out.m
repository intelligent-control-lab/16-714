function [tlist, xlist, ulist, log] = roll_out(model, ctrl, x0, sim, opts)
% ROLL_OUT  Unified simulator for CT/DT models, state/measurement/estimate control.
%
% Usage (typical):
%   [t,x,u] = roll_out("unicycle3", @(x,t)[1;0.2], [0;0;0], ...
%                      struct('type','CT','T',10,'integrator','ode45'));
%
%   sim.type       : 'CT' | 'DT'
%   sim.T          : horizon (seconds) if type='CT'
%   sim.K          : number of steps if type='DT'
%   sim.dt         : step size (seconds) if type='DT'
%   sim.integrator : 'ode45'|'Euler'|'ZOH'|'RK4'|'Direct'
%   sim.stop       : optional @(xk,tk,log)->true to stop early (DT only)
%
%   ctrl           : controller; either function handle u(x,t)  (state feedback),
%                    or struct('mode','state'|'measurement'|'estimate','pi', @pi)
%
%   opts.sensor      : optional @(x,t)->y  (default: y=x if needed)
%   opts.observer    : optional struct with fields:
%                        xhat0   : initial estimate (vector)
%                        update  : @(xhat,u,y,t)->xhat_next
%   opts.log_fields  : cellstr of {'y','xhat'} etc. to record
%   opts.cost      : optional @(xk,uk,xnext,tk)->scalar cost (DT only);
%                      if provided, per-step costs are stored in log.cost(k+1)
%   opts.verbose     : logical (default false)
%
% Outputs:
%   tlist : 1×N time stamps
%   xlist : n×N states
%   ulist : m×(N-1) inputs
%   log   : struct with requested logs (e.g., log.y, log.xhat, log.cost)
%
% Notes:
%   • Time-varying dynamics are handled by letting model.f_dt/f_ct optionally
%     accept time:  f_dt(x,u,t) / f_ct(x,u,t). If they don't, we call f_dt(x,u)
%     / f_ct(x,u) transparently.
%   • CT simulation currently assumes state-feedback (u = pi(x,t)) for clarity.
%   • Dependencies: resolve_dynamics.m, step.m (expects dyn.f_dt(x,u) for DT).

    arguments
        model
        ctrl
        x0 double {mustBeVector}
        sim struct
        opts struct = struct()
    end

    % --- Normalize opts once: inject defaults for missing fields ---
    defaults = struct( ...
        'sensor',      [], ...
        'observer',    struct(), ...
        'log_fields',  {{}}, ...
        'cost',      [], ...     % per-step cost handle (DT only)
        'verbose',     false );
    opts = mergestruct(defaults, opts);

    % Normalize initial state shape
    x0 = x0(:);

    % Resolve model into a unified interface
    if isstring(model) || ischar(model)
        dyn = resolve_dynamics(model);
    else
        dyn = model;
    end

    % Normalize controller
    [ctrl_mode, pi_handle] = normalize_controller(ctrl);

    % Initialize logging struct
    log = struct();
    want_y     = any(strcmp(opts.log_fields, 'y'));
    want_xhat  = any(strcmp(opts.log_fields, 'xhat'));
    has_cost = isa(opts.cost, 'function_handle');

    switch upper(sim.type)
        % =========================
        % Continuous-time simulation
        % =========================
        case 'CT'
            T = require_field(sim, 'T', 'CT simulation requires sim.T (horizon in seconds).');

            if ~strcmp(ctrl_mode, 'state')
                error(['For tutorial clarity, CT simulation currently supports state-feedback only. ', ...
                       'Use DT with Euler/ZOH/RK4 for measurement/estimate feedback.']);
            end

            % ODE RHS: allow model.f_ct to be with/without time
            odefn = @(t, x) call_f_ct_with_optional_t(dyn, x, pi_handle(x,t), t);

            switch lower(default_if_missing(sim, 'integrator', 'ode45'))
                case {'','ode45'}
                    opts_ode = odeset('RelTol',1e-6,'AbsTol',1e-8);
                    soln = ode45(odefn, [0, T], x0, opts_ode);
                    tlist = soln.x;
                    xlist = soln.y;

                    % inputs aligned with samples except terminal
                    m = size(pi_handle(x0,0),1);
                    ulist = zeros(m, max(1,numel(tlist)-1));
                    for i = 1:numel(tlist)-1
                        ulist(:,i) = pi_handle(xlist(:,i), tlist(i));
                    end
                otherwise
                    error('CT simulation: integrator must be ''ode45'' (for now).');
            end

        % =======================
        % Discrete-time simulation
        % =======================
        case 'DT'
            K  = require_field(sim, 'K',  'DT simulation requires sim.K (#steps).');
            dt = require_field(sim, 'dt', 'DT simulation requires sim.dt (step size).');
            integrator = lower(default_if_missing(sim, 'integrator', 'Euler'));
            has_sens   = isfield(opts, 'sensor') && isa(opts.sensor, 'function_handle');

            % Initialize histories
            tlist      = zeros(1, K+1);
            xlist      = zeros(numel(x0), K+1);  xlist(:,1) = x0;

            % Determine input dimension from first call
            u0 = pi_handle(x0, 0);
            m     = numel(u0);
            ulist = zeros(m, K);   % u_0 ... u_{K-1}

            % Optional: observer
            if isfield(opts, 'observer') && isfield(opts.observer, 'xhat0')
                xhat  = opts.observer.xhat0(:);
                obs_f = opts.observer.update;
                if want_xhat, log.xhat = zeros(numel(xhat), K); end
            else
                xhat  = [];
                obs_f = [];
            end

            % Optional logs
            if want_y,    log.y        = []; end
            if has_cost, log.cost  = zeros(1, K); end

            % Main DT loop
            for k = 0:(K-1)
                tk = k*dt;
                tlist(k+1) = tk;
                xk = xlist(:, k+1);

                % Measurement
                if has_sens
                    yk = opts.sensor(xk, tk);
                else
                    yk = xk;
                end
                if want_y
                    if isempty(log.y), log.y = zeros(numel(yk), K); end
                    log.y(:,k+1) = yk;
                end

                % Control mode
                if strcmp(ctrl_mode,'state')
                    uk = pi_handle(xk, tk);
                elseif strcmp(ctrl_mode,'measurement')
                    uk = pi_handle(yk, tk);
                elseif strcmp(ctrl_mode,'estimate')
                    if isempty(obs_f)
                        error('Estimate feedback requested but opts.observer.update is missing.');
                    end
                    if k > 0
                        xhat = obs_f(xhat, ulist(:,k), yk, tk);
                    end
                    if want_xhat
                        log.xhat(:,k+1) = xhat;
                    end
                    uk = pi_handle(xhat, tk);
                else
                    error('Unknown ctrl mode "%s".', ctrl_mode);
                end

                % Store input
                ulist(:,k+1) = uk;

                % Advance one step
                % Build a time-baked dynamics wrapper so step() can keep calling dyn.f_dt(x,u)
                dyn_k = wrap_dyn_with_time(dyn, tk);  % f_dt(x,u) -> calls f_dt(x,u,tk) if available
                switch integrator
                    case {'euler','zoh','rk4','direct'}
                        x_next = step(xk, uk, dt, dyn_k, integrator);
                    otherwise
                        error('Unknown DT integrator "%s".', integrator);
                end

                % Store next state
                xlist(:,k+2) = x_next;

                % cost (after x and x_next are known)
                if has_cost
                    log.cost(1, k+1) = opts.cost(xk, uk, x_next, tk);
                end

                % Early stopping (optional)
                if isfield(sim, 'stop') && ~isempty(sim.stop) && isa(sim.stop, 'function_handle')
                    if sim.stop(x_next, tk+dt, log)
                        tlist(k+2) = tk+dt;
                        % Trim histories to the point of stopping
                        tlist = tlist(1:k+2);
                        xlist = xlist(:,1:k+2);
                        ulist = ulist(:,1:k+1);
                        if want_y,     log.y       = log.y(:,1:k+1);    end
                        if want_xhat,  log.xhat    = log.xhat(:,1:k+1); end
                        if has_cost, log.cost  = log.cost(1,1:k+1); end
                        return;
                    end
                end
            end

            % finalize time stamps
            tlist(end) = K*dt;

        otherwise
            error('sim.type must be ''CT'' or ''DT''.');
    end
end

% --------------------- helpers (local to this file) ----------------------

function dyn_t = wrap_dyn_with_time(dyn, t)
% Returns a struct with f_dt(x,u) and (optionally) f_ct(x,u) that internally
% call the original with time if available, else without time.
    dyn_t = dyn;

    if isfield(dyn, 'f_dt') && ~isempty(dyn.f_dt) && isa(dyn.f_dt,'function_handle')
        dyn_t.f_dt = @(x,u) call_f_dt_with_optional_t(dyn, x, u, t);
    end
    if isfield(dyn, 'f_ct') && ~isempty(dyn.f_ct) && isa(dyn.f_ct,'function_handle')
        dyn_t.f_ct = @(x,u) call_f_ct_with_optional_t(dyn, x, u, t);
    end
end

function fx = call_f_dt_with_optional_t(dyn, x, u, t)
% Try f_dt(x,u,t) and fall back to f_dt(x,u)
    try
        fx = dyn.f_dt(x, u, t);
    catch
        fx = dyn.f_dt(x, u);
    end
end

function fx = call_f_ct_with_optional_t(dyn, x, u, t)
% Try f_ct(x,u,t) and fall back to f_ct(x,u)
    try
        fx = dyn.f_ct(x, u, t);
    catch
        fx = dyn.f_ct(x, u);
    end
end

function val = require_field(s, fname, errmsg)
    if ~isfield(s, fname) || isempty(s.(fname))
        error(errmsg);
    end
    val = s.(fname);
end

function val = default_if_missing(s, fname, default_val)
    if ~isfield(s, fname) || isempty(s.(fname))
        val = default_val;
    else
        val = s.(fname);
    end
end

function [mode, pi_handle] = normalize_controller(ctrl)
    if isa(ctrl, 'function_handle')
        mode = 'state';
        pi_handle = ctrl;
    elseif isstruct(ctrl)
        if ~isfield(ctrl, 'mode') || isempty(ctrl.mode), ctrl.mode = 'state'; end
        if ~isfield(ctrl, 'pi') || isempty(ctrl.pi)
            error('Controller struct requires field .pi (the policy handle).');
        end
        mode = lower(ctrl.mode);
        pi_handle = ctrl.pi;
    else
        error('ctrl must be a function handle u(x,t) or a struct with fields mode, pi.');
    end
end

function s = mergestruct(a,b)
    s = a;
    if ~isstruct(b), return; end
    fn = fieldnames(b);
    for i = 1:numel(fn)
        s.(fn{i}) = b.(fn{i});
    end
end
