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
        dyn.getAB = @(varargin) getAB(name, varargin{:});

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
            dyn.getAB = @(varargin) dynamics.getAB(varargin{:});
        end
    else
        error('Unsupported dynamics type. Use name (string/char) or struct with handles.');
    end

    if isempty(dyn.f_ct) && isempty(dyn.f_dt)
        error('Dynamics must provide CT RHS (.f_ct) or DT map (.f_dt).');
    end
end
