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
        [bar_A, bar_B] = lift_dynamics(sys.A, sys.B, sys.N);
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
        c.xref = @(x) bar_A*x + bar_B*c.uref(x); %return all predicted states; note this implementation could be improved by only calling one quadprog
        c.u = @(x,t) get_first_u(c.uref, x, m);
    case 'ILCw' % Frequency domain ILC
        % PS = P/(1+PC) = b / (a + b * C);
        c.name = 'ILCw';
        c.L.b = add(sys.a, conv(sys.b, -sys.c));
        c.L.a = sys.b;
        c.Q.b = [1]; c.Q.a = [1]; % Design Q to ensure robustness when there is noise
        c.u = @(error, ffold) filter(c.Q.b, c.Q.a, ffold + iclfilter(c.L.b, c.L.a, error));
    case 'ILCt' % Time domain ILC
        c.name = 'ILCt';
        Acl = sys.A - sys.B*sys.K; % get closed loop A matrix
        [~, bar_B] = lift_dynamics(Acl, sys.B, sys.N);
        clist = mat2cell(repmat(pinv(sys.C),1,sys.N+1)',ones(1,sys.N+1)); % assume single output
        bar_inv_C = blkdiag(clist{:})';
        c.L = pinv(bar_B) * bar_inv_C;
        c.u = @(error, ffold) ffold + error * c.L';
    case 'MRACx' % Model reference adaptive control with full state feedback
        c = MRAC(sys,'record');
end
end

function u = get_first_u(uref, x, m)
    ulist = uref(x);
    u = ulist(1:m);
end

% Note the difference between this filter and the original matlab filter
% function is that this filter allows non-causal transformation
function y = iclfilter(b, a, error) 
    nshift = length(b) - length(a);
    y = filter(b, a, error);
    y(1:end-nshift) = y(nshift+1:end);
end

function c = add(a, b)
n = max(length(a),length(b));
c = zeros(n,1);
for i = n:-1:1
    if n-i < length(a)
        c(i) = c(i) + a(end - n + i);
    end
    if n-i < length(b)
        c(i) = c(i) + b(end - n + i);
    end
end
end

function [bar_A, bar_B] = lift_dynamics(A,B,N)
n = size(A, 1);
bar_A = eye(n);
bar_B = [];
for i=1:N
    bar_B = [zeros(n, i); [bar_A * B bar_B]];
    bar_A = [eye(n); bar_A * A];
end
end
