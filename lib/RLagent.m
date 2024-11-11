classdef RLagent < handle
    properties
        control % the control law derived from the Q function
        W % the parameter for the Q function
        ft % the features that depend on state and control
        Q % the Q function
        gradQ % the gradient of the Q function
        V % the Value fuction
        tcondition % Termination condition for each episode
    end
    methods
        % Initialization of the Learning Agent
        function obj = RLagent(sys) 
            obj.ft = @(x,u) [x;u];
            obj.W = sys.W0;
            obj.Q = @(x,u,W) obj.ft(x,u)'*W*obj.ft(x,u)/2;
            obj.gradQ = @(x,u) obj.ft(x,u)*obj.ft(x,u)'/2;
            nx = length(sys.x0);
            obj.V = @(x,W) x'*(W(1:nx,1:nx)-W(1:nx,nx+1:end)*inv(W(nx+1:end, nx+1:end))*W(nx+1:end, 1:nx))*x/2;
            obj.tcondition = sys.tcondition;
            obj.control = @(x,W) greedy(x,W,sys.epsilon);
        end
        
        % Update for one episode
        function [x_list,u_list,l_list] = update(obj,alg,alpha,sys)
            switch alg
                case 'MC' % Monte Carlo
                    u = @(x,t) obj.control(x,obj.W);

                    [t_list,x_list,u_list,l_list] = roll_out(sys, u, sys.x0, 'DT', obj.tcondition, sys.dt,'Direct');

                    dW = zeros(size(obj.W));
                    for k = 1:length(t_list)-1
                        G = sum(l_list(k:end));
                        x = x_list(:,k); u = u_list(:,k);
                        dW = dW + alpha * (G - obj.Q(x,u,obj.W)) * obj.gradQ(x,u);
                    end
                    obj.W = obj.W + dW;

                case 'SARSA'
                    x = sys.x0; x_list = [x];
                    u = obj.control(x, obj.W); u_list = [u];
                    k = 0; l_list = [];
                    while ~obj.tcondition(x, k/sys.dt)
                        xnew = sys.f(x, u); l = sys.h(x,u);
                        unew = obj.control(xnew, obj.W); 
                        obj.W = obj.W + alpha * (l + sys.delta * obj.Q(xnew,unew,obj.W) - obj.Q(x,u,obj.W)) * obj.gradQ(x,u);
                        x = xnew; x_list = [x_list x];
                        u = unew; u_list = [u_list u];
                        k = k+1; l_list = [l_list l];
                    end

                case 'QLearning'
                    x = sys.x0; x_list = [x];
                    u_list = [];
                    k = 0; l_list = [];
                    while ~obj.tcondition(x, k/sys.dt)
                        u = obj.control(x, obj.W); u_list = [u_list u];
                        xnew = sys.f(x, u); l = sys.h(x,u);
                        obj.W = obj.W + alpha * (l + sys.delta * obj.V(xnew,obj.W) - obj.Q(x,u,obj.W)) * obj.gradQ(x,u);
                        x = xnew; x_list = [x_list x];
                        k = k+1; l_list = [l_list l];
                    end
            end
        end

    end
end

function u = greedy(x, W, epsilon)
sample = rand(1);
nx = length(x);
nu = size(W,1)-nx;
u = - inv(W(nx+1:end, nx+1:end))*(W(nx+1:end, 1:nx) + W(1:nx, nx+1:end))/2*x;
if sample < epsilon
    u = -rand(nu)*x;
end
end