%% ===== make an estimator update closure with persistent Z =====
function upd = estimator_update_handle(est_obj)
% Returns upd(xhat,u,y,t) that keeps its own covariance-like state.
    Z_init = est_obj.Z(:,:,1);   % initial steady-state (or user-provided) value
    function xhat_next = upd_impl(xhat,u,y,~)
        % keep a private Z across calls (persistent)
        persistent Zp
        if isempty(Zp), Zp = Z_init; end
        xhat_next = est_obj.update_x(xhat, u, y, Zp);
        Zp        = est_obj.update_Z(Zp);
    end
    upd = @upd_impl;
end