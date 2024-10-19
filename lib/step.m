function xnew = step(x,u,dt,dynamics,mode)
switch mode
    case 'ZOH'
        if isa(dynamics, 'string')
            odefn = @(t,y) f(y,u,dynamics);
        else
            odefn = @(t,y) dynamics.f(y,u);
        end
        soln = ode45(odefn,[0,dt],x);
        xnew = soln.y(:,end);
    case 'Euler' 
        if isa(dynamics, 'string')
            xnew = x + f(x,u,dynamics).*dt;
        else
            xnew = x + dynamics.f(x,u).*dt;
        end
    case 'Direct' % the dynamics are directly specified in discrete time
        xnew = dynamics.f(x,u);
end
end