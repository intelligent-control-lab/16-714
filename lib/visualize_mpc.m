function feasible = visualize_mpc(sys,c,Kmax)
feasible = 1;
dt = sys.dt;
tlist = 0:dt:dt*Kmax;
n = length(sys.x0);
xlist = zeros(n,Kmax+1); xlist(:,1) = sys.x0;
xlist_steps = zeros(n, sys.N+1, Kmax+1);
try
    m = size(c.u(sys.x0,1),1);
catch
    error(['no solution for the initial state ',num2str(sys.x0')]);
end
ulist = zeros(m,Kmax);
for k = 1:Kmax
    try
        ulist(:,k) = c.u(xlist(:,k),tlist(k));
        xlist_steps(:,:,k) = reshape(c.xref(xlist(:,k)), n, sys.N+1);
    catch
        warning(['no solution at time step ',num2str(k)]);
        feasible = -1;
        ulist(:,k) = zeros(m,1);
    end
    xlist(:,k+1) = step(xlist(:,k), ulist(:,k), sys.dt, sys.name, 'Euler');
end

figure; hold on
plot(0:Kmax, xlist(1,:), 'k', 'LineWidth', 2);
string = ["Executed Trajectory"];
for t = 1:Kmax
    plot(t-1:t+sys.N-1, xlist_steps(1,:,t), "color", [1-t/Kmax, t/Kmax, 1-t/Kmax]);
    string = [string, "MPC Time "+num2str(t-1)];
end
box on
legend(string)
end