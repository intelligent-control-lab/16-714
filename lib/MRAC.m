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
        function obj = MRAC(sys, varargin)
            obj.name  = 'MRACx';

            % Build reference system and simulate with new roll_out (DT)
            refsys.f_dt = @(x,u) sys.Astar * x;

            % Use the same horizon and step as the plant
            sim_ref = struct('type','DT', 'K',sys.N, 'dt',sys.dt, 'integrator','Direct');
            opts    = struct();

            % Reference input is zero
            [~, obj.ref] = roll_out(refsys, @(x,t) 0, sys.x0, sim_ref, opts);

            % Initialize estimates
            obj.Ahat(:,:,1) = sys.Ahat;
            obj.Bhat(:,:,1) = sys.Bhat;
            obj.F           = sys.F;

            obj.n     = size(sys.Ahat,1);
            obj.dt    = sys.dt;
            obj.Astar = sys.Astar;

            % Optional recording
            if ~isempty(varargin) && strcmp(varargin{1}, 'record')
                obj.data.Ahat(:,:,1)  = obj.Ahat;
                obj.data.Bhat(:,:,1)  = obj.Bhat;
                obj.data.F(:,:,1)     = obj.F;
                obj.data.phi          = [];
                obj.data.error(:,1)   = sys.x0;
            end
        end

        function input = u(obj, x, t)
            % Adapt parameters after the first step
            if t > 0
                % Align reference index with time t and dt
                kref   = max(1, min(size(obj.ref,2), round(t/obj.dt)+1));
                error  = x - obj.ref(:, kref);

                % RLS-style update of covariance and parameters
                % (F is the information inverse; ensure dimensions consistent)
                obj.F    = inv( inv(obj.F) + obj.phi*obj.phi.' );
                obj.Ahat = obj.Ahat + error * (obj.phi.' * obj.F(:, 1:obj.n));
                obj.Bhat = obj.Bhat + error * (obj.phi.' * obj.F(:, obj.n+1:end));

                if ~isempty(obj.data); obj.refresh(error); end
            end

            % Control law K = Bhat^\dagger (A* - Ahat)
            % Use pinv for robustness if Bhat is poorly conditioned
            K     = pinv(obj.Bhat) * (obj.Astar - obj.Ahat);
            input = K * x;

            % Regressor for the next update
            obj.phi = [x; input];
        end

        function refresh(obj, error)
            obj.data.Ahat(:,:,end+1) = obj.Ahat;
            obj.data.Bhat(:,:,end+1) = obj.Bhat;
            obj.data.F(:,:,end+1)    = obj.F;
            obj.data.phi(:,end+1)    = obj.phi;
            obj.data.error(:,end+1)  = error;
        end

        function analysis(obj, A, B, N)
            for k = 1:N
                ABtilde = [obj.data.Ahat(:,:,k) - A,  obj.data.Bhat(:,:,k) - B];

                if k > 1
                    phi_prev = obj.data.phi(:, k-1);
                    Fk       = obj.data.F(:,:,k);
                    obj.data.error_pos(:,k) = (1 - phi_prev.' * Fk * phi_prev) * obj.data.error(:,k);
                else
                    obj.data.error_pos(:,k) = obj.data.error(:,k);
                end

                obj.data.V(k) = norm(obj.data.error_pos(:,k))^2 + trace(ABtilde * inv(obj.data.F(:,:,k)) * ABtilde.');
            end
        end
    end
end
