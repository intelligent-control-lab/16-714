% when mode is CT
% varargin only specifies the maximum simulation time
%
% when mode is DT
% varargin specifies 1) maximum simulation step; 2) dt; 3) discretization
% mode (ZOH or Euler)

function [tlist, xlist, ulist] = roll_out(dynamics, u, x0, mode, varargin)
switch mode
    case 'CT'
        odefn = @(t,y) f(y,u(y,t),dynamics);
        opts = odeset('RelTol',1e-6,'AbsTol',1e-5);
        soln = ode45(odefn,[0,varargin{1}],x0,opts);
        tlist = soln.x;
        xlist = soln.y(:,:);
        for i = 1:length(tlist)-1 ulist(:,i) = u(xlist(:,i),tlist(i)); end
    case 'DT'
        dt = varargin{2};
        tlist = 0:dt:dt*varargin{1};
        xlist = zeros(size(x0,1),length(tlist)); xlist(:,1) = x0;
        try
            m = size(u(x0,1),1);
        catch
            error(['no solution for the initial state ',num2str(x0)]);
        end
        ulist = zeros(m,length(tlist)-1);
        if length(varargin)>2
            simmode = varargin{3};
        else
            simmode = 'Euler';
        end
        for k = 1:varargin{1}
            try
                ulist(:,k) = u(xlist(:,k),tlist(k));
            catch
                warning(['no solution at time step ',num2str(k)]);
                ulist(:,k) = zeros(m,1);
            end
            xlist(:,k+1) = step(xlist(:,k), ulist(:,k), dt, dynamics, simmode);
        end
        
end