function feasible = visualize_mpc(sys, c, Kmax, integrator)
% VISUALIZE_MPC  Execute MPC for Kmax steps and overlay per-step predictions.
% Uses: step(x,u,dt,dyn,mode)  with dyn from resolve_dynamics(sys.name).
%
% c must provide:
%   u(x,t)    -> control at state/time
%   xref(x0)  -> predicted states over horizon (n x (N+1))
%
% integrator (optional): 'Euler' | 'ZOH' | 'RK4' | 'Direct' (default: 'ZOH')

    if nargin < 4 || isempty(integrator), integrator = 'ZOH'; end

    feasible = 1;
    dt  = sys.dt;
    nx  = length(sys.x0);
    tks = 0:dt:dt*Kmax;

    % Resolve dynamics once; step(...) expects dyn struct
    dyn = resolve_dynamics(sys.name);

    % Infer input dimension
    try
        u0 = c.u(sys.x0, 0);
        m  = size(u0,1);
    catch
        error('visualize_mpc:init', ...
              'No solution for the initial state [%s].', num2str(sys.x0', '%g '));
    end

    % Preallocate
    xlist       = zeros(nx, Kmax+1);  xlist(:,1) = sys.x0;
    ulist       = zeros(m,  Kmax);
    xlist_steps = NaN(nx, sys.N+1, Kmax+1);  % store predicted rollouts (if available)

    % Main loop
    for k = 1:Kmax
        try
            % control at current state/time
            uk = c.u(xlist(:,k), tks(k));
            ulist(:,k) = uk;

            % controller's predicted trajectory from current state (optional)
            try
                xlist_steps(:,:,k) = reshape(c.xref(xlist(:,k)), nx, sys.N+1);
            catch
                % leave NaNs if xref is unavailable
            end

        catch
            warning('visualize_mpc:stepFail', ...
                    'No solution at time step k=%d. Using u_k = 0.', k);
            feasible    = -1;
            ulist(:,k)  = zeros(m,1);
        end

        % propagate one step with the requested integrator
        xlist(:,k+1) = step(xlist(:,k), ulist(:,k), dt, dyn, integrator);
    end

    % (Optional) final prediction from terminal executed state
    try
        xlist_steps(:,:,Kmax+1) = reshape(c.xref(xlist(:,Kmax+1)), nx, sys.N+1);
    catch
        % ignore if not provided
    end

    % -------- Plot executed trajectory + rolling predictions (state 1) --------
    figure('Name','MPC Execution and Rolling Predictions'); hold on; box on; grid on;
    plot(0:Kmax, xlist(1,:), 'k', 'LineWidth', 2, 'DisplayName','Executed Trajectory');

    leg = {'Executed Trajectory'};
    for t = 1:Kmax
        if all(~isnan(xlist_steps(1,:,t)))
            plot(t-1:t+sys.N-1, xlist_steps(1,:,t), ...
                'Color', [1 - t/Kmax, t/Kmax, 1 - t/Kmax], ...
                'DisplayName', sprintf('MPC Time %d', t-1));
            leg{end+1} = sprintf('MPC Time %d', t-1); 
        end
    end
    xlabel('k'); ylabel('x_1'); legend(leg, 'Location','best');
end
