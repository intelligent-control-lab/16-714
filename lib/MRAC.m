classdef MRAC < handle
    properties
        name
        Ahat
        Bhat
        F
        phi
        ref
        dt
        Astar
        n
        data
    end

    methods
        function obj = MRAC(sys,varargin)
            obj.name = 'MRACx';
            refsys.f = @(x,u) sys.Astar * x;
            [~,obj.ref] = roll_out(refsys, @(x,t) 0, sys.x0, 'DT', sys.N, sys.dt, 'Direct');
            obj.Ahat(:,:,1) = sys.Ahat; 
            obj.Bhat(:,:,1) = sys.Bhat;
            obj.F = sys.F;
            obj.n = size(sys.Ahat,1);
            obj.dt = sys.dt;
            obj.Astar = sys.Astar;
            if ~isempty(varargin) && strcmp(varargin{1}, 'record')
                obj.data.Ahat(:,:,1) = obj.Ahat;
                obj.data.Bhat(:,:,1) = obj.Bhat;
                obj.data.F(:,:,1) = obj.F;
                obj.data.phi = [];
                obj.data.error(:,1) = sys.x0;
            end
        end

        function input = u(obj,x,t)
            if t > 0
                error = x - obj.ref(:,round(t/obj.dt)+1);
                obj.F = inv(inv(obj.F) + obj.phi*obj.phi');              
                obj.Ahat = obj.Ahat + error*obj.phi'*obj.F(:,1:obj.n);
                obj.Bhat = obj.Bhat + error*obj.phi'*obj.F(:,obj.n+1:end);
                if ~isempty(obj.data) obj.refresh(error); end
            end
            K = pinv(obj.Bhat)*(obj.Astar - obj.Ahat);
            input = K*x;
            obj.phi = [x;input]; %This is to be used in the next update
        end

        function refresh(obj,error)
            obj.data.Ahat(:,:,end+1) = obj.Ahat;
            obj.data.Bhat(:,:,end+1) = obj.Bhat;
            obj.data.F(:,:,end+1) = obj.F;
            obj.data.phi(:,end+1) = obj.phi;
            obj.data.error(:,end+1) = error;
        end

        function analysis(obj, A, B, N)
            for k = 1:N
                ABtilde = [obj.data.Ahat(:,:,k) - A obj.data.Bhat(:,:,k) -  B];
                if k > 1
                    obj.data.error_pos(:,k) = (1-obj.data.phi(:,k-1)'*obj.data.F(:,:,k)*obj.data.phi(:,k-1))*obj.data.error(:,k);
                else
                    obj.data.error_pos(:,k) = obj.data.error(:,k);
                end
                obj.data.V(k) = norm(obj.data.error_pos(:,k))^2 + trace(ABtilde * inv(obj.data.F(:,:,k)) * ABtilde');
            end
        end

    end
end