classdef RLagent < handle
    properties
        control   % control law derived from the Q function
        W         % parameter for the Q function
        ft        % features that depend on state and control
        Q         % Q function
        gradQ     % gradient of the Q function wrt W
        V         % Value function (greedy w.r.t. u)
        tcondition
        sim       % cached sim struct for roll_out
        opts      % cached opts for roll_out
    end
    methods
        % Initialization of the Learning Agent
        function obj = RLagent(sys)
            obj.ft    = @(x,u) [x;u];
            obj.W     = sys.W0;
            obj.Q     = @(x,u,W) obj.ft(x,u)'*W*obj.ft(x,u)/2;
            obj.gradQ = @(x,u)  (obj.ft(x,u)*obj.ft(x,u)')/2;
            nx        = numel(sys.x0);
            obj.V     = @(x,W) x'*(W(1:nx,1:nx) - W(1:nx,nx+1:end)/W(nx+1:end,nx+1:end)*W(nx+1:end,1:nx))*x/2;
            obj.tcondition = sys.tcondition;

            % epsilon-greedy control (closure captures sys.epsilon)
            obj.control = @(x,W) greedy(x,W,sys.epsilon);

            % cache a DT sim config compatible with roll_out
            obj.sim = struct( ...
                'type','DT', ...
                'K', sys.kmax, ...
                'dt', sys.dt, ...
                'integrator','Direct', ...
                'stop', @(xk,tk,log) obj.tcondition(xk,tk) ...
            );
            obj.opts = struct();
        end

        % Update for one episode
        function [x_list,u_list,l_list] = update(obj, alg, alpha, sys)
            switch alg
                case 'MC' % Monte Carlo with episode roll_out
                    ufun = @(x,t) obj.control(x, obj.W);
                    [t_list, x_list, u_list, log] = roll_out(sys, ufun, sys.x0, obj.sim, obj.opts);

                    % recover l_list (stage costs) — from log if present, else compute
                    if isstruct(log) && isfield(log,'l')
                        l_list = log.l;
                    else
                        K = size(x_list,2);
                        l_list = zeros(1,K-1);
                        for k = 1:K-1
                            l_list(k) = sys.h(x_list(:,k), u_list(:,k));
                        end
                    end

                    % MC weight update
                    dW = zeros(size(obj.W));
                    for k = 1:numel(l_list)
                        G = sum(l_list(k:end));                 % return from k
                        xk = x_list(:,k); uk = u_list(:,k);
                        dW = dW + alpha * (G - obj.Q(xk,uk,obj.W)) * obj.gradQ(xk,uk);
                    end
                    obj.W = obj.W + dW;

                case 'SARSA'
                    x = sys.x0; x_list = x; u_list = []; l_list = []; k = 0;
                    u = obj.control(x, obj.W);
                    while ~obj.tcondition(x, k*sys.dt)
                        % one-step transition under DT dynamics
                        xnew = sys.f_dt(x, u);
                        l    = sys.h(x, u);
                        unew = obj.control(xnew, obj.W);
                        obj.W = obj.W + alpha * (l + sys.delta*obj.Q(xnew,unew,obj.W) - obj.Q(x,u,obj.W)) * obj.gradQ(x,u);
                        % log & advance
                        x = xnew; u = unew; k = k+1;
                        x_list = [x_list, x];
                        u_list = [u_list, obj.control(x_list(:,end-1), obj.W)]; 
                        l_list = [l_list, l]; 
                    end

                case 'QLearning'
                    x = sys.x0; x_list = x; u_list = []; l_list = []; k = 0;
                    while ~obj.tcondition(x, k*sys.dt)
                        u = obj.control(x, obj.W);
                        xnew = sys.f_dt(x, u);
                        l    = sys.h(x, u);
                        obj.W = obj.W + alpha * (l + sys.delta*obj.V(xnew,obj.W) - obj.Q(x,u,obj.W)) * obj.gradQ(x,u);
                        % log & advance
                        x = xnew; k = k+1;
                        x_list = [x_list, x];    
                        u_list = [u_list, u];    
                        l_list = [l_list, l];    
                    end
            end
        end
    end
end

function u = greedy(x, W, epsilon)
    sample = rand(1);
    nx = numel(x);
    Guu = W(nx+1:end, nx+1:end);
    Gux = W(nx+1:end, 1:nx);
    Gxu = W(1:nx, nx+1:end);
    % greedy minimizer for quadratic Q(x,u) = 1/2 [x;u]^T W [x;u]
    u_star = - (Guu \ ((Gux + Gxu')/2)) * x;
    if sample < epsilon
        % exploratory linear action around state
        u = u_star + (-rand(size(Guu,1), nx))*x;
    else
        u = u_star;
    end
end