classdef PGagent < handle
    properties
        W % parameter for the value function
        theta % parameter for the policy
        gradP
        V
        gradV
        control
        tcondition
    end
    methods
        function obj = PGagent(sys)
            switch sys.policy.name
                case 'Gaussian'
                    obj.theta = sys.policy.init;
                    obj.W = sys.W0;
                    % stochastic linear-Gaussian in gain space (scalar case)
                    obj.control = @(x,mu,sigma) mu*x + randn() * sigma*x;
                    % grad wrt [mu; sigma] (stacked)
                    obj.gradP = @(x,u,mu,sigma) [ (u-mu*x)/(sigma*x)^2 * x; ...
                                                  -x + (u-mu*x)^2/(sigma*x)^2 * x];
                    obj.V = @(x,W) x'*W*x/2;
                    obj.gradV = @(x) x^2/2;
                    obj.tcondition = sys.tcondition;
                case 'Softmax'
                    % (not used here)
            end
        end

        function [x_list,u_list,l_list] = update(obj,alg,alpha,sys)
            % shared sim/opts for DT rollout
            simDT  = struct('type','DT','K',sys.kmax,'dt',sys.dt,'integrator','Direct', ...
                            'stop', @(xk,tk,log) obj.tcondition(xk,tk));
            optsDT = struct();
            optsDT.cost     = @(x,u,xnext,t) sys.h(x,u);

            switch alg
                case 'REINFORCE'
                    ufun = @(x,t) obj.control(x,obj.theta.mu,obj.theta.sigma);
                    [t_list,x_list,u_list,log] = roll_out(sys, ufun, sys.x0, simDT, optsDT);
                    l_list = log.cost;

                    dtheta = [zeros(size(obj.theta.mu)); zeros(size(obj.theta.sigma))];
                    for k = 1:length(t_list)-1
                        G = sum(l_list(k:end));
                        xk = x_list(:,k); uk = u_list(:,k);
                        dtheta = dtheta - alpha * G * obj.gradP(xk,uk,obj.theta.mu,obj.theta.sigma);
                    end
                    obj.theta.mu    = obj.theta.mu    + dtheta(1:size(obj.theta.mu,1), :);
                    obj.theta.sigma = obj.theta.sigma + dtheta(size(obj.theta.mu,1)+1:end,:);

                case 'REINFORCEbl'
                    ufun = @(x,t) obj.control(x,obj.theta.mu,obj.theta.sigma);
                    [t_list,x_list,u_list,log] = roll_out(sys, ufun, sys.x0, simDT, optsDT);
                    l_list = log.cost;

                    dtheta = [zeros(size(obj.theta.mu)); zeros(size(obj.theta.sigma))];
                    dW = zeros(size(obj.W));
                    for k = 1:length(t_list)-1
                        G = sum(l_list(k:end));
                        xk = x_list(:,k); uk = u_list(:,k);
                        d  = G - obj.V(xk,obj.W);
                        dtheta = dtheta - alpha.p * G * obj.gradP(xk,uk,obj.theta.mu,obj.theta.sigma);
                        dW     = dW     + alpha.v * d * obj.gradV(xk);
                    end
                    obj.theta.mu    = obj.theta.mu    + dtheta(1:size(obj.theta.mu,1), :);
                    obj.theta.sigma = obj.theta.sigma + dtheta(size(obj.theta.mu,1)+1:end,:);
                    obj.W = obj.W + dW;

                case 'ActorCritic'
                    % on-policy TD(0)-style update without full roll_out call
                    x = sys.x0; x_list = x; u_list = []; l_list = [];
                    k = 0;
                    while ~obj.tcondition(x, k/sys.dt)
                        u = obj.control(x,obj.theta.mu,obj.theta.sigma);
                        u_list = [u_list u];
                        xnew = sys.f_dt(x,u,k/sys.dt);
                        l    = sys.h(x,u);
                        l_list = [l_list l];

                        d  = l + sys.delta*obj.V(xnew,obj.W) - obj.V(x,obj.W);
                        D  = obj.gradP(x,u,obj.theta.mu,obj.theta.sigma);
                        obj.theta.mu    = obj.theta.mu    - alpha.p * d * D(1:size(obj.theta.mu,1), :);
                        obj.theta.sigma = obj.theta.sigma - alpha.p * d * D(size(obj.theta.mu,1)+1:end,:);

                        obj.W = obj.W + alpha.v * d * obj.gradV(x);

                        x = xnew; x_list = [x_list x];
                        k = k + 1;
                    end
            end
        end
    end
end