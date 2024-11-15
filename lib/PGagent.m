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
                    obj.control = @(x,mu,sigma) mu*x + randn() * sigma*x; 
                    obj.gradP = @(x,u,mu,sigma) [(u-mu*x)/(sigma*x)^2*x; 
                                            -x + (u-mu*x)^2/(sigma*x)^2*x];
                    obj.V = @(x,W) x'*W*x/2;
                    obj.gradV = @(x) x^2/2;
                    obj.tcondition = sys.tcondition;
                case 'Softmax'
            end

        end

        function [x_list,u_list,l_list] = update(obj,alg,alpha,sys)
            switch alg
                case 'REINFORCE'
                    u = @(x,t) obj.control(x,obj.theta.mu,obj.theta.sigma);
                    [t_list,x_list,u_list,l_list] = roll_out(sys, u, sys.x0, 'DT', obj.tcondition, sys.dt,'Direct');
                    
                    dtheta = [zeros(size(obj.theta.mu));zeros(size(obj.theta.sigma))];
                    for k = 1:length(t_list)-1
                        G = sum(l_list(k:end));
                        x = x_list(:,k); u = u_list(:,k);
                        dtheta = dtheta - alpha * G * obj.gradP(x,u,obj.theta.mu,obj.theta.sigma);
                    end
                    
                    obj.theta.mu = obj.theta.mu + dtheta(1:size(obj.theta.mu,1), :);
                    obj.theta.sigma = obj.theta.sigma + dtheta(size(obj.theta.mu,1)+1:end,:);

                case 'REINFORCEbl'
                    u = @(x,t) obj.control(x,obj.theta.mu,obj.theta.sigma);
                    [t_list,x_list,u_list,l_list] = roll_out(sys, u, sys.x0, 'DT', obj.tcondition, sys.dt,'Direct');

                    dtheta = [zeros(size(obj.theta.mu));zeros(size(obj.theta.sigma))];
                    dW = zeros(size(obj.W));
                    for k = 1:length(t_list)-1
                        G = sum(l_list(k:end));
                        x = x_list(:,k); u = u_list(:,k);
                        d = G - obj.V(x,obj.W);
                        dtheta = dtheta - alpha.p * G * obj.gradP(x,u,obj.theta.mu,obj.theta.sigma);
                        dW = dW + alpha.v * d * obj.gradV(x);
                    end
                    
                    obj.theta.mu = obj.theta.mu + dtheta(1:size(obj.theta.mu,1), :);
                    obj.theta.sigma = obj.theta.sigma + dtheta(size(obj.theta.mu,1)+1:end,:);
                    obj.W = obj.W + dW;
                case 'ActorCritic'
                    x = sys.x0; x_list = [x];
                    u_list = [];
                    k = 0; l_list = [];
                    while ~obj.tcondition(x, k/sys.dt)
                        u = obj.control(x,obj.theta.mu,obj.theta.sigma);u_list = [u_list u];
                        xnew = sys.f(x, u); l = sys.h(x,u);
                        
                        d = l + sys.delta*obj.V(xnew,obj.W) - obj.V(x,obj.W);
                        D = obj.gradP(x,u,obj.theta.mu,obj.theta.sigma);
                        obj.theta.mu = obj.theta.mu - alpha.p * d * D(1:size(obj.theta.mu,1), :);
                        obj.theta.sigma = obj.theta.sigma - alpha.p * d * D(size(obj.theta.mu,1)+1:end,:);
                        
                        obj.W = obj.W + alpha.v * d * obj.gradV(x);

                        x = xnew; x_list = [x_list x];
                        k = k+1; l_list = [l_list l];
                    end
            end
        end
    end
end