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
%   opts.sensor    : optional @(x,u,t)->y  (default: y=x if needed)
%   opts.observer  : optional struct with fields:
%                      xhat0   : initial estimate (vector)
%                      update  : @(xhat,u,y,t)->xhat_next
%   opts.timevarying : optional @(dyn,k,t)->dyn_k  (DT only)
%   opts.log_fields  : cellstr of {'y','xhat'} etc. to record
%   opts.verbose     : logical (default false)
%
% Outputs:
%   tlist : 1×N time stamps
%   xlist : n×N states
%   ulist : m×(N-1) inputs
%   log   : struct with requested logs (e.g., log.y, log.xhat)
%
% Notes:
%   • CT simulation currently assumes state-feedback (u = pi(x,t)) for clarity.
%     For measurement/estimate feedback, prefer DT with Euler/ZOH/RK4.
%   • Models can be provided by name (string/char) or as a struct of handles.
%
% Dependencies (below in this file set or separate files):
%   resolve_dynamics.m, step.m (updated), local helpers in this file.

    arguments
        model
        ctrl
        x0 double {mustBeVector}
        sim struct
        opts.sensor = []
        opts.observer = []
        opts.timevarying = []
        opts.log_fields = {}
        opts.verbose (1,1) logical = false
    end

    % Normalize initial state shape
    x0 = x0(:);

    % Resolve model into a unified interface
    dyn = resolve_dynamics(model);

    % Normalize controller
    [ctrl_mode, pi_handle] = normalize_controller(ctrl);

    % Initialize logging struct
    log = struct();
    want_y    = any(strcmp(opts.log_fields, 'y'));
    want_xhat = any(strcmp(opts.log_fields, 'xhat'));

    switch upper(sim.type)
        % =========================
        % Continuous-time simulation
        % =========================
        case 'CT'
            T = require_field(sim, 'T', 'CT simulation requires sim.T (horizon in seconds).');

            % Tutorial-friendly constraint for clarity:
            if ~strcmp(ctrl_mode, 'state')
                error(['For tutorial clarity, CT simulation currently supports state-feedback only. ', ...
                       'Use DT with Euler/ZOH/RK4 for measurement/estimate feedback.']);
            end
            % Build RHS for ODE45
            switch lower(sim.integrator)
                case {'','ode45'}
                    odefn = @(t, x) dyn.f_ct(x, pi_handle(x,t));  % u = pi(x,t)
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
            has_tv     = ~isempty(opts.timevarying) && isa(opts.timevarying, 'function_handle');
            has_sens   = ~isempty(opts.sensor)     && isa(opts.sensor, 'function_handle');

            % Initialize histories
            tlist      = zeros(1, K+1);
            xlist      = zeros(numel(x0), K+1);  xlist(:,1) = x0;

            % Determine input dimension from first call (robustly)
            try
                u0 = pi_handle(x0, 0);
            catch
                % If measurement/estimate, we still try with identity y=x / xhat=x
                u0 = pi_handle(x0, 0);
            end
            m     = numel(u0);
            ulist = zeros(m, K);   % inputs u_0 ... u_{K-1}

            % Optional: observer
            if ~isempty(opts.observer)
                xhat  = opts.observer.xhat0(:);
                obs_f = opts.observer.update;
                if want_xhat, log.xhat = zeros(numel(xhat), K); end
            else
                xhat  = [];
                obs_f = [];
            end

            % Optional: measurement log
            if want_y, log.y = []; end

            % Main DT loop
            for k = 0:(K-1)
                tk = k*dt;
                tlist(k+1) = tk;
                xk = xlist(:, k+1);

                % Time-varying dynamics (per step) if provided
                if has_tv
                    dyn_k = opts.timevarying(dyn, k, tk);
                else
                    dyn_k = dyn;
                end

                % Measurement
                if strcmp(ctrl_mode,'state')
                    uk = pi_handle(xk, tk);

                elseif strcmp(ctrl_mode,'measurement')
                    if has_sens
                        yk = opts.sensor(xk, (k>0)*ulist(:,k) + 0, tk); % pass u_{k-1} if needed
                    else
                        yk = xk; % default measurement: identity
                    end
                    if want_y
                        if isempty(log.y), log.y = zeros(numel(yk), K); end
                        log.y(:,k+1) = yk;
                    end
                    uk = pi_handle(yk, tk);

                elseif strcmp(ctrl_mode,'estimate')
                    % measurement
                    if has_sens
                        yk = opts.sensor(xk, (k>0)*ulist(:,k) + 0, tk);
                    else
                        yk = xk;
                    end
                    if want_y
                        if isempty(log.y), log.y = zeros(numel(yk), K); end
                        log.y(:,k+1) = yk;
                    end
                    % observer update
                    if isempty(obs_f)
                        error('Estimate feedback requested but opts.observer.update is missing.');
                    end
                    xhat = obs_f(xhat, (k>0)*ulist(:,k) + 0, yk, tk);
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
                switch integrator
                    case {'euler','zoh','rk4','direct'}
                        x_next = step(xk, uk, dt, dyn_k, integrator);
                    otherwise
                        error('Unknown DT integrator "%s".', integrator);
                end

                xlist(:,k+2) = x_next;

                % Early stopping (optional)
                if isfield(sim, 'stop') && ~isempty(sim.stop) && isa(sim.stop, 'function_handle')
                    if sim.stop(x_next, tk+dt, log)
                        tlist(k+2) = tk+dt;
                        % Trim histories to the point of stopping
                        tlist = tlist(1:k+2);
                        xlist = xlist(:,1:k+2);
                        ulist = ulist(:,1:k+1);
                        if want_y,    log.y    = log.y(:,1:k+1); end
                        if want_xhat, log.xhat = log.xhat(:,1:k+1); end
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