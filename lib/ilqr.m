function ilqr(model, x0, xg, opts)
% iLQR with:
% - First-iteration nominal states from linear interpolation between x0 and xg
% - roll_out.m for all forward passes
% - resolve_dynamics.getAB for linearization (only works for unicycle3 for now)
% - f.m for computing defect d_k using the same DT step as simulation
%
% Usage:
%   ilqr("unicycle3", [0;0;0], [5;5;pi/2], struct('N',100,'dt',0.1))
% ---------- options ----------
if nargin < 4, opts = struct(); end
N        = get_opt(opts, 'N', 100);        % horizon steps
dt       = get_opt(opts, 'dt', 0.1);
Q        = get_opt(opts, 'Q', diag([1,1,0.1]));
R        = get_opt(opts, 'R', diag([1,0.5]));
S        = get_opt(opts, 'S', diag([50,50,10]));
u_ref    = get_opt(opts, 'u_ref', zeros(2,1));
alphas   = get_opt(opts, 'alphas', [1.0 0.5 0.25 0.1 0.05]);
maxIter  = get_opt(opts, 'maxIter', 100);
tolJ     = get_opt(opts, 'tolJ', 1e-1);
lambda0  = get_opt(opts, 'lambda0', 1e-8);
lambdaMax= get_opt(opts, 'lambdaMax', 1e8);
integr   = get_opt(opts, 'integrator', 'ZOH');   % 'ZOH'|'Euler'|'RK4'|'Direct'
improve_tol = 1e-12;

% ---------- prepare dynamics interface ----------
dyn = resolve_dynamics(model);

% ---------- initial nominal control (zeros) ----------
nu = numel(u_ref);
ubar = zeros(nu, N);

% ---------- initial nominal state trajectory: linear interpolation ----------
% xbar(:,1) = x0, xbar(:,N+1) = xg, and linearly spaced in-between.
nx = numel(x0);
xbar = zeros(nx, N+1);
for k = 0:N
    tau = k / N;
    xbar(:,k+1) = (1 - tau)*x0 + tau*xg;
end

% (Optional) If you prefer a dynamically consistent starting nominal, you can
% uncomment the following 4 lines; but we keep the linear reference per request.
% ctrl_replay = @(x,t) replay_u(t, ubar, dt);
% simDT = struct('type','DT','K',N,'dt',dt,'integrator',integr);
% [~, xbar, ubar_] = roll_out(model, ctrl_replay, x0, simDT);
% ubar = ubar_;  % keep zeros normally

% ---------- initial cost ----------
J = total_cost(Q,R,S, xbar, ubar, xg, u_ref);
fprintf('Iter %2d: J = %.6f\n', 0, J);
prevJ = J;

lambda = lambda0;

% ---------- main iLQR loop ----------
for it = 1:maxIter

    % ----- build A_k, B_k, d_k and Δ-cost linear terms -----
    A = cell(1,N); B = cell(1,N); d = cell(1,N);
    q = cell(1,N); r = cell(1,N);
    for k = 1:N
        xk = xbar(:,k);
        uk = ubar(:,k);

        % pointwise linearization via resolve_dynamics (DT, ZOH)
        [A{k}, B{k}] = dyn.getAB('DT', struct('dim',nx,'dt',dt,'dmode','ZOH','x',xk,'u',uk));

        % defect d_k computed with SAME discrete step as used for linearization
        xk_plus = step(xk, uk, dt, dyn, 'Euler');
        d{k} = xk_plus - xbar(:,k+1);

        % Δ-cost linear terms (tracking to xg and u_ref)
        q{k} = Q*(xbar(:,k) - xg);
        r{k} = R*(uk - u_ref);
    end

    % ----- backward pass (P_k, s_k) & gains (K_k, k_k) -----
    P = cell(1,N+1); s = cell(1,N+1);
    K = cell(1,N);  kfeed = cell(1,N);
    P{N+1} = S;     s{N+1} = zeros(nx,1);
    diverged = false;

    for k = N:-1:1
        Ak = A{k};  Bk = B{k}; dk = d{k};
        Pn = P{k+1}; sn = s{k+1};
        qk = q{k};   rk = r{k};

        Quu = R + Bk.'*Pn*Bk + lambda*eye(size(R));
        Qux = Bk.'*Pn*Ak;
        gu  = rk + Bk.'*(Pn*dk + sn);   % includes defect d_k

        [L,p] = chol(Quu,'lower');
        if p>0
            diverged = true; break;
        end
        invQuu = L'\(L\eye(size(Quu)));

        % gains
        K{k}     = -invQuu * Qux;
        kfeed{k} = -invQuu * gu;

        % Riccati-like recursion (PMP form)
        P{k} = Q + Ak.'*Pn*Ak - Qux.'*invQuu*Qux;
        s{k} = qk + Ak.'*(Pn*dk + sn) - Qux.'*invQuu*gu;
        P{k} = 0.5*(P{k}+P{k}.'); % symmetrize
    end

    if diverged
        lambda = min(lambda*10, lambdaMax);
        if lambda >= lambdaMax
            warning('Backward pass failed (Quu not PD). Stopping.');
            break;
        end
        continue
    end

    % ----- forward pass / line search (each candidate via roll_out) -----
    simDT = struct('type','DT','K',N,'dt',dt,'integrator',integr);
    best = struct('J', inf, 'alpha', NaN, 'x', [], 'u', []);
    for a = alphas
        ctrl_affine = @(x,t) affine_ctrl(t, x, xbar, ubar, kfeed, K, dt, a);
        [~, xnew, unew] = roll_out(model, ctrl_affine, x0, simDT);
        Jcand = total_cost(Q,R,S, xnew, unew, xg, u_ref);
        if Jcand < best.J
            best.J = Jcand; best.alpha = a; best.x = xnew; best.u = unew;
        end
    end

    % ----- accept / reject -----
    if (J - best.J) > max(improve_tol, 1e-12*abs(J))
        xbar = best.x; ubar = best.u; J = best.J;
        lambda = max(lambda/5, 1e-12);
    else
        lambda = min(lambda*10, lambdaMax);
        if lambda >= lambdaMax
            warning('No improvement from line search. Stopping.');
            break;
        end
        continue
    end

    fprintf('Iter %2d: J = %.6f  (alpha=%.2f, lambda=%.1e)\n', it, J, best.alpha, lambda);
    if abs(prevJ - J) < tolJ, break; end
    prevJ = J;
end

% ---------- reports & plots ----------
fprintf('Final state: [%g %g %g]^T\n', xbar(:,end));
fprintf('Goal error:  [%g %g %g]^T\n', xbar(:,end)-xg);

figure; plot(xbar(1,:), xbar(2,:), 'LineWidth', 2); hold on;
plot(x0(1), x0(2), 'ko','MarkerFaceColor','k'); plot(xg(1), xg(2),'r*','MarkerSize',10);
axis equal; grid on; xlabel('x'); ylabel('y'); title(sprintf('%s iLQR (PMP) trajectory', model));

figure;
subplot(2,1,1); plot(0:N-1, ubar(1,:), 'LineWidth', 1.5); grid on; ylabel('u_1 = v');
subplot(2,1,2); plot(0:N-1, ubar(2,:), 'LineWidth', 1.5); grid on; ylabel('u_2 = \\omega'); xlabel('k');

% ============================================================
% helpers
% ============================================================
function val = get_opt(s, field, default)
    if isfield(s, field), val = s.(field); else, val = default; end
end

function u = replay_u(t, U, dt)
    % Replay open-loop U(:,k) as a function of time t (DT).
    % k = 1 at t in [0,dt). Use floor for DT alignment.
    k = 1 + floor(t/dt);
    k = max(1, min(size(U,2), k));
    u = U(:,k);
end

function u = affine_ctrl(t, x, xbar, ubar, kfeed, K, dt, alpha)
    % u_k = ubar_k + alpha*k_k + K_k * (x - xbar_k)
    k = 1 + floor(t/dt);
    k = max(1, min(size(ubar,2), k));
    u = ubar(:,k) + alpha*kfeed{k} + K{k}*(x - xbar(:,k));
end

function J = total_cost(Q,R,S, x, u, xg, uref)
    Nloc = size(u,2);
    J = 0;
    for kk = 1:Nloc
        dx = x(:,kk)   - xg;
        du = u(:,kk)   - uref;
        J = J + dx.'*Q*dx + du.'*R*du;
    end
    dT = x(:,Nloc+1) - xg;
    J = J + dT.'*S*dT;
end

end
