function dyn = resolve_dynamics(dynamics)
% RESOLVE_DYNAMICS  Normalize dynamics spec into a uniform interface.
% Output fields:
%   dyn.f_ct(x,u) : CT RHS (for CT sim, Euler/ZOH/RK4)
%   dyn.f_dt(x,u) : DT map (for 'Direct')
%   dyn.h         : optional output map
%   dyn.getAB     : optional linearization hook
%   dyn.name      : name or []

    dyn = struct('f_ct',[], 'f_dt',[], 'h',[], 'getAB',[], 'name',[]);

    if isstring(dynamics) || ischar(dynamics)
        name   = string(dynamics);
        dyn.name  = name;
        dyn.f_ct  = @(x,u) f(x,u,name);                 % your provided name-based f
        dyn.f_dt  = [];                                 % unless you define a name-based fd
        dyn.h     = [];                                 % default none
        dyn.getAB = @(mode, opts) getAB_local(name, mode, opts);

    elseif isstruct(dynamics)
        if isfield(dynamics,'f')  && isa(dynamics.f, 'function_handle')
            dyn.f_ct = @(x,u) dynamics.f(x,u);
        end
        if isfield(dynamics,'fd') && isa(dynamics.fd,'function_handle')
            dyn.f_dt = @(x,u) dynamics.fd(x,u);
        end
        if isfield(dynamics,'h')  && isa(dynamics.h, 'function_handle')
            dyn.h    = dynamics.h;
        end
        if isfield(dynamics,'getAB') && isa(dynamics.getAB,'function_handle')
            dyn.getAB = @(mode, opts) getAB_local(name, mode, opts);
        end
    else
        error('Unsupported dynamics type. Use name (string/char) or struct with handles.');
    end

    if isempty(dyn.f_ct) && isempty(dyn.f_dt)
        error('Dynamics must provide CT RHS (.f_ct) or DT map (.f_dt).');
    end
end

% ========================================================================
% Local helper: [A,B] = getAB_local(name, mode, opts)
% - name: 'single_integrator' | 'double_integrator'
% - mode: 'CT' | 'DT'  (char or string ok)
% - opts: struct with fields:
%       .dim   (required, positive integer)
%       .dt    (required for 'DT', positive)
%       .dmode (optional, default 'ZOH'; supports 'ZOH','EULER')
% ========================================================================
function [A,B] = getAB_local(name, mode, opts)
    % ---- Robust parsing (accept char or string) -------------------------
    if isstring(name), name = char(name); end
    if isstring(mode), mode = char(mode); end

    % ---- Validate opts struct ------------------------------------------
    if ~(isstruct(opts) && isfield(opts,'dim'))
        error('getAB:opts', 'opts must be a struct with at least opts.dim.');
    end
    dim = opts.dim;
    if ~(isscalar(dim) && dim>0 && mod(dim,1)==0)
        error('getAB:dim', 'opts.dim must be a positive integer.');
    end

    isCT = strcmpi(mode,'CT');
    isDT = strcmpi(mode,'DT');
    if ~(isCT || isDT)
        error('getAB:mode', 'mode must be ''CT'' or ''DT'' (got %s).', mode);
    end

    if isDT
        if ~isfield(opts,'dt') || ~(isscalar(opts.dt) && isfinite(opts.dt) && opts.dt>0)
            error('getAB:dt', 'opts.dt must be provided and positive for DT mode.');
        end
        dt = opts.dt;
        if isfield(opts,'dmode') && ~isempty(opts.dmode)
            dmode = upper(string(opts.dmode));
        else
            dmode = "ZOH";
        end
    end

    switch lower(name)
        case 'single_integrator'
            % xdot = u;  x in R^dim, u in R^dim
            if isCT
                A = zeros(dim, dim);
                B = eye(dim, dim);
            else
                % ZOH == Forward Euler here
                A = eye(dim, dim);
                B = dt * eye(dim, dim);
            end

        case 'double_integrator'
            % x = [p; v], dim = 2n, u in R^n
            n = dim/2;
            if mod(n,1) ~= 0
                error('getAB:dim', 'double_integrator expects even state dimension; got dim=%d.', dim);
            end
            if isCT
                A = [zeros(n,n) eye(n,n);
                     zeros(n,n) zeros(n,n)];
                B = [zeros(n,n);
                     eye(n,n)];
            else
                switch char(dmode)
                    case 'ZOH'
                        A = [eye(n,n) dt*eye(n,n);
                             zeros(n,n) eye(n,n)];
                        B = [0.5*dt^2*eye(n,n);
                             dt*eye(n,n)];
                    case {'EULER','FORWARDEULER','FE'}
                        A_ct = [zeros(n,n) eye(n,n);
                                zeros(n,n) zeros(n,n)];
                        B_ct = [zeros(n,n);
                                eye(n,n)];
                        A = eye(2*n) + dt*A_ct;
                        B = dt*B_ct;
                    otherwise
                        warning('getAB:dmode','Unknown dmode="%s"; falling back to ZOH.', char(dmode));
                        A = [eye(n,n) dt*eye(n,n);
                             zeros(n,n) eye(n,n)];
                        B = [0.5*dt^2*eye(n,n);
                             dt*eye(n,n)];
                end
            end

        otherwise
            warning('getAB:unknown','Unknown dynamics "%s"; returning zero A,B.', name);
            A = zeros(dim, dim);
            B = zeros(dim, dim);
    end
end