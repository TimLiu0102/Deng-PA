function X = AO_X(S, X, W, theta, phi, scene, params)
% 论文位置子问题：按波导逐条更新 + L-BFGS-B风格方向 + 投影
N = params.N;
for n = 1:N
    xn = X(n,:).';
    g = num_grad(n,xn,S,X,W,theta,phi,scene,params);

    Hinv = eye(numel(xn)); % L-BFGS-B style初始逆Hessian近似
    d = -Hinv*g;
    if norm(d) < 1e-12
        continue;
    end

    f0 = local_obj(n,xn,S,X,W,theta,phi,scene,params);
    alpha = params.line_search_alpha0;
    xbar = xn;
    for it = 1:params.line_search_max_iter
        xt = xn + alpha*d;
        xt = Constraint_Checker('project', xt, params);
        ft = local_obj(n,xt,S,X,W,theta,phi,scene,params);
        if ft >= f0
            xbar = xt; break;
        else
            alpha = alpha*params.line_search_beta;
        end
    end

    xhat = Constraint_Checker('project', xbar, params);
    fhat = local_obj(n,xhat,S,X,W,theta,phi,scene,params);
    if fhat >= f0 + params.eps_X
        s = xhat - xn;
        g2 = num_grad(n,xhat,S,X,W,theta,phi,scene,params);
        y = g2 - g;
        Hinv = bfgs_inv_update(Hinv,s,y); %#ok<NASGU>
        X(n,:) = xhat.';
    end
end
end

function f = local_obj(n,xn,S,X,W,theta,phi,scene,params)
Xt = X; Xt(n,:) = xn.';
f = Signal_model('sum_rate', S, Xt, W, theta, phi, scene, params);
end

function g = num_grad(n,xn,S,X,W,theta,phi,scene,params)
M = numel(xn); g = zeros(M,1);
for i=1:M
    e = zeros(M,1); e(i)=1;
    xp = Constraint_Checker('project', xn + params.step_fd*e, params);
    xm = Constraint_Checker('project', xn - params.step_fd*e, params);
    fp = local_obj(n,xp,S,X,W,theta,phi,scene,params);
    fm = local_obj(n,xm,S,X,W,theta,phi,scene,params);
    g(i) = (fp-fm)/(2*params.step_fd);
end
end

function Hn = bfgs_inv_update(H,s,y)
rho = 1/(y'*s + 1e-12);
I = eye(size(H));
Hn = (I-rho*s*y')*H*(I-rho*y*s') + rho*(s*s');
end
