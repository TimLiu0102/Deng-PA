function W = AO_W(S, X, theta, phi, W, scene, params)
% 论文WMMSE子问题：固定(S,X,theta,phi)更新W
params.scene = scene;
H = Channel_model('composite_channels', S, X, theta, phi, [], params);
Nt = size(H,1); K = size(H,2);
if isempty(W)
    W = H; W = W/max(norm(W,'fro'),1e-12)*sqrt(params.P_max);
end

R_prev = Signal_model('sum_rate', S, X, W, theta, phi, scene, params);
for it = 1:params.I_W
    % W -> u
    u = zeros(K,1);
    for k=1:K
        hk = H(:,k);
        denom = sum(abs(hk' * W).^2) + params.sigma2;
        u(k) = (hk' * W(:,k))/denom;
    end

    % u -> v
    v = zeros(K,1);
    for k=1:K
        hk = H(:,k);
        e = abs(u(k))^2*(sum(abs(hk' * W).^2)+params.sigma2) - 2*real(u(k)*hk'*W(:,k)) + 1;
        v(k) = 1/max(real(e),1e-12);
    end

    % v -> W (mu二分满足功率约束)
    A = zeros(Nt,Nt); b = zeros(Nt,K);
    for j=1:K
        hj = H(:,j);
        A = A + v(j)*abs(u(j))^2*(hj*hj');
        b(:,j) = v(j)*conj(u(j))*hj;
    end
    W = solve_with_mu(A,b,params.P_max);

    R_now = Signal_model('sum_rate', S, X, W, theta, phi, scene, params);
    if abs(R_now - R_prev) < params.eps_W
        break;
    end
    R_prev = R_now;
end

end

function W = solve_with_mu(A,b,Pmax)
Nt = size(A,1); I = eye(Nt);
W0 = (A + 0*I) \ b;
if real(trace(W0*W0')) <= Pmax
    W = W0; return;
end
muL = 0; muU = 1;
while real(trace(((A+muU*I)\b)*((A+muU*I)\b)')) > Pmax
    muU = 2*muU;
    if muU > 1e8, break; end
end
for it=1:50
    mu = 0.5*(muL+muU);
    Wm = (A+mu*I)\b;
    p = real(trace(Wm*Wm'));
    if p > Pmax, muL = mu; else, muU = mu; end
end
W = (A+muU*I)\b;
end
