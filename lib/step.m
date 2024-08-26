function xnew = step(x,u,dt,name,mode)
switch mode
    case 'ZOH'
        odefn = @(t,y) f(y,u,name);
        soln = ode45(odefn,[0,dt],x);
        xnew = soln.y(:,end);
    case 'Euler' 
        xnew = x + f(x,u,name).*dt;
end
end