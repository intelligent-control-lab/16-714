function c = synthesis(name, sys)
switch name
    case 'LQR' % infinite horizon LQR
        c.name = 'LQR';
        if sys.dt <= 0 % Continuous time system
            [K, P] = lqr(sys.A, sys.B, sys.Q, sys.R);
        else
            [K, P] = dlqr(sys.A, sys.B, sys.Q, sys.R);
        end
        c.u = @(x,t) -K*x;
        c.K = K;
        c.P = P;
    case 'LQRn' % finite horizon LQR
        c.name = 'LQRn';
        nx = length(sys.x0); nu = size(sys.B,2);
        P = zeros(nx, nx, sys.N+1);
        P(:,:,sys.N+1) = sys.S;
        K = zeros(nu, nx, sys.N);
        for i = sys.N:-1:1
            K(:,:,i) = inv(sys.R + sys.B' * P(:,:,i+1) * sys.B) * sys.B' * P(:,:,i+1) * sys.A;
            P(:,:,i) = sys.Q + sys.A' * P(:,:,i+1) * sys.A - sys.A' * P(:,:,i+1) * sys.B * K(:,:,i);
        end
        c.u = @(x,t) -K(:,:,t/sys.dt+1)*x;
        c.K = K;
        c.P = P;
    case 'nMPC' % unconstrained linear quadratic MPC
        c.name = 'nMPC';
        nx = length(sys.x0); nu = size(sys.B,2);
        P = zeros(nx, nx, sys.N+1);
        P(:,:,sys.N+1) = sys.S;
        K = zeros(nu, nx, sys.N);
        for i = sys.N:-1:1
            K(:,:,i) = inv(sys.R + sys.B' * P(:,:,i+1) * sys.B) * sys.B' * P(:,:,i+1) * sys.A;
            P(:,:,i) = sys.Q + sys.A' * P(:,:,i+1) * sys.A - sys.A' * P(:,:,i+1) * sys.B * K(:,:,i);
        end
        c.u = @(x,t) -K(:,:,1)*x;
        c.K = K;
        c.P = P;
    case 'uMPC' % control constrained linear quadratic MPC
        c.name = 'uMPC';
        nx = length(sys.x0); nu = size(sys.B,2);
        P = zeros(nx, nx, sys.N+1);
        P(:,:,sys.N+1) = sys.S;
        K = zeros(nu, nx, sys.N);
        for i = sys.N:-1:1
            K(:,:,i) = inv(sys.R + sys.B' * P(:,:,i+1) * sys.B) * sys.B' * P(:,:,i+1) * sys.A;
            P(:,:,i) = sys.Q + sys.A' * P(:,:,i+1) * sys.A - sys.A' * P(:,:,i+1) * sys.B * K(:,:,i);
        end
        c.u = @(x,t) max(min(-K(:,:,1)*x, sys.umax), sys.umin);
        c.K = K;
        c.P = P;
    case 'qMPC' % general quadratic MPC with both control and state constraints
        % Extended dynamics
        n = size(sys.A,1);
        m = size(sys.B,2);
        bar_A = eye(n);bar_B = [];
        for i=1:sys.N
            bar_B = [zeros(n, i); [bar_A * sys.B bar_B]];
            bar_A = [eye(n); bar_A * sys.A];
        end
        % State and control constraints
        max_U = repmat(sys.umax, sys.N, 1); min_U = repmat(sys.umin, sys.N, 1);
        max_X = repmat(sys.xmax, sys.N+1, 1); min_X = repmat(sys.xmin, sys.N+1, 1);
        % Extended costs
        bar_Q = kron(eye(sys.N+1), sys.Q);
        bar_Q(end-1:end, end-1:end) = sys.S;
        bar_R = kron(eye(sys.N), sys.R);
        QQ = bar_B'*bar_Q*bar_B + bar_R;
        AA = [bar_B;-bar_B;eye(sys.N);-eye(sys.N)];
        bb = @(x) [max_X - bar_A*x; -min_X + bar_A*x; max_U; -min_U];
        if isfield(sys,'xterminal')
            AA = [AA;[zeros(size(sys.xterminal.A,1),n*(sys.N)) sys.xterminal.A]*bar_B];
            bb = @(x) [max_X - bar_A*x; -min_X + bar_A*x; max_U; -min_U; sys.xterminal.b - [zeros(size(sys.xterminal.A,1),n*(sys.N)) sys.xterminal.A]*bar_A*x];
        end
        CC = @(x) (bar_A*x)'*bar_Q*bar_B;
        c.uref = @(x) quadprog(QQ, CC(x)', AA, bb(x)); % return all controls
        c.xref = @(x) bar_A*x + bar_B*c.uref(x); %return all predicted states;          note this implementation could be improved by only calling one quadprog
        c.u = @(x,t) get_first_u(c.uref, x, m);
end
end

function u = get_first_u(uref, x, m)
    ulist = uref(x);
    u = ulist(1:m);
end