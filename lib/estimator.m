function est = estimator(name, sys)
est.xhat(:,1) = sys.x0hat;
est.Z(:,:,1) = sys.X0;
switch name
    case 'KF' % KF
        % Update covariance
        M = @(Z) sys.A*Z*sys.A' + sys.Bw*sys.W*sys.Bw';
        est.update_Z = @(Z) M(Z) - M(Z)*sys.C'*inv(sys.V+sys.C*M(Z)*sys.C')*sys.C*M(Z);
        % Update estimate
        xprior = @(xhat,u) sys.fnominal(xhat,u);
        est.KF_gain = @(Z) M(Z)*sys.C'*inv(sys.V+sys.C*M(Z)*sys.C');
        error = @(y,xhat,u) y - sys.C*xprior(xhat,u);
        est.update_x = @(xhat,u,y,Z) xprior(xhat, u) + est.KF_gain(Z) * error(y,xhat,u);
    case 'KFss' % Steady State KF
        % Update covariance
        Ms = dare(sys.A', sys.C', sys.Bw*sys.W*sys.Bw', sys.V);
        est.update_Z = @(Z) Ms - Ms*sys.C'*inv(sys.V+sys.C*Ms*sys.C')*sys.C*Ms;
        % Update estimate
        xprior = @(xhat,u) sys.fnominal(xhat,u);
        est.KF_gain = Ms*sys.C'*inv(sys.V+sys.C*Ms*sys.C');
        error = @(y,xhat,u) y - sys.C*xprior(xhat,u);
        est.update_x = @(xhat,u,y,Z) xprior(xhat, u) + est.KF_gain * error(y,xhat,u);
    case 'EKF' % Note this function takes in an external linearization method
        % Update covariance
        M = @(Z,xhat,u) sys.A(xhat,u)*Z*sys.A(xhat,u)' + sys.Bw*sys.W*sys.Bw';
        est.update_Z = @(Z,xhat,u) M(Z,xhat,u) - M(Z,xhat,u)*sys.C(xhat)'*inv(sys.C(xhat)*M(Z,xhat,u)*sys.C(xhat)'+sys.V)*sys.C(xhat)*M(Z,xhat,u);
        % Update estimate
        xprior = @(xhat,u) sys.fnominal(xhat, u);
        est.KF_gain = @(Z,xhat,u) M(Z,xhat,u)*sys.C(xhat)'*inv(sys.V+sys.C(xhat)*M(Z,xhat,u)*sys.C(xhat)');
        error = @(y,xhat,u) y - sys.hnominal(xhat);
        est.update_x = @(xhat,u,y,Z) xprior(xhat, u) + est.KF_gain(Z,xhat,u) * error(y,xhat,u);
    case 'UKF'
        %Dynamic update
        xlist = @(xhat, Z, u) dynamic_update_sigma_points(xhat, Z, u, sys);
        xprior = @(xhat, Z, u) weighted_mean(xlist(xhat, Z, u));
        varprior = @(xhat, Z, u) weighted_variance(xlist(xhat, Z, u));
        M = @(Z, xhat, u) varprior(xhat, Z, u) + sys.Bw*sys.W*sys.Bw';
        %Measurement update
        ylist = @(xhat, Z, u) measurement_update_sigma_points(xprior(xhat, Z, u), M(Z, xhat, u), sys);
        yprior = @(xhat, Z, u) weighted_mean(ylist(xhat, Z, u));
        yvar = @(xhat, Z, u) weighted_variance(ylist(xhat, Z, u));
        xyvar = @(xhat, Z, u) weighted_cov(xprior(xhat, Z, u), M(Z, xhat, u), ylist(xhat, Z, u));

        est.KF_gain = @(Z,xhat,u) xyvar(xhat, Z, u) * inv(sys.V + yvar(xhat, Z, u));
        error = @(y,xhat,u,Z) y - yprior(xhat, Z, u);
        est.update_x = @(xhat,u,y,Z) xprior(xhat, Z, u) + est.KF_gain(Z,xhat,u) * error(y,xhat,u,Z);
        est.update_Z = @(Z,xhat,u) M(Z, xhat, u) - xyvar(xhat, Z, u) * inv(sys.V + yvar(xhat, Z, u)) * xyvar(xhat, Z, u)';
end
end

function xlist = dynamic_update_sigma_points(xhat, Z, u, sys)
    xlist = get_sigma_points(xhat, Z);
    for j = 1:length(xlist)
        x = xlist{j}.x;
        xlist{j}.x = sys.fnominal(x, u);
    end
end

function ylist = measurement_update_sigma_points(xprior, M, sys)
    xlist = get_sigma_points(xprior, M);
    ylist = cell(length(xlist),1);
    for j = 1:length(xlist)
        x = xlist{j}.x;
        ylist{j}.x = sys.hnominal(x);
        ylist{j}.w = xlist{j}.w;
    end
end

function mean = weighted_mean(xlist)
    sum = zeros(length(xlist{1}.x),1);
    weight = 0;
    for j = 1:length(xlist)
        sum = sum + xlist{j}.x * xlist{j}.w;
        weight = weight + xlist{j}.w;
    end
    mean = sum ./ weight;
end
% Note these two functions could be improved
function var = weighted_variance(xlist)
    sum = zeros(length(xlist{1}.x),1);
    weight = 0;
    for j = 1:length(xlist)
        sum = sum + xlist{j}.x * xlist{j}.w;
        weight = weight + xlist{j}.w;
    end
    mean = sum ./ weight;
    weight = 0;
    sum = zeros(length(xlist{1}.x));
    for j = 1:length(xlist)
        sum = sum + (xlist{j}.x - mean) * (xlist{j}.x - mean)' * xlist{j}.w;
        weight = weight + xlist{j}.w;
    end
    var = sum ./ weight;
end

function covariance = weighted_cov(xprior, M, ylist)    
    xlist = get_sigma_points(xprior, M);    
    sum = zeros(length(xlist{1}.x),1);
    weight = 0;
    for j = 1:length(xlist)
        sum = sum + xlist{j}.x * xlist{j}.w;
        weight = weight + xlist{j}.w;
    end
    mx = sum ./ weight;
    sum = zeros(length(ylist{1}.x),1);
    weight = 0;
    for j = 1:length(ylist)
        sum = sum + ylist{j}.x * ylist{j}.w;
        weight = weight + ylist{j}.w;
    end
    my = sum ./ weight;
    sum = zeros(length(xlist{1}.x), length(ylist{1}.x));
    weight = 0;
    for j = 1:length(xlist)
        dx = xlist{j}.x - mx;
        dy = ylist{j}.x - my;
        sum = sum + dx * dy' * xlist{j}.w;
        weight = weight + xlist{j}.w;
    end
    covariance = sum ./ weight;
end

function xlist = get_sigma_points(mean, var)
    kappa = 2; n = length(mean);
    xlist = cell(1+2*n,1);
    xlist{1}.x = mean;
    xlist{1}.w = kappa / (n+kappa);
    
    std = sqrtm((n+kappa).*var);
    for i = 1:n
        xlist{2*i}.x = mean + std(:,i);
        xlist{2*i}.w = 1/2/(n+kappa);
        xlist{2*i+1}.x = mean - std(:,i);
        xlist{2*i+1}.w = 1/2/(n+kappa);
    end
end