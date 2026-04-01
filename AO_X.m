function X_new = AO_X(params, scene, model, state)
% AO_X：固定 (S,W,theta,phi) 后，按波导逐条更新位置 X

if nargin < 3
    model = struct(); %#ok<NASGU>
end

X_new = state.X;
N = size(X_new,1);

% 逐条波导更新，不并行
for n = 1:N
    X_new = update_single_waveguide(n, X_new, params, scene, state);
end

end

%% ======================== 内部子函数 ========================
function X_cur = update_single_waveguide(n, X_cur, params, scene, state)
% 第n条波导的位置子问题更新
x_n = X_cur(n,:).';

% 局部目标值 f_n(x_n)
f_old = local_objective_fn(n, x_n, X_cur, params, scene, state);

% 数值微分近似梯度
g = numerical_gradient_position(n, x_n, X_cur, params, scene, state);

% L-BFGS-B style方向（简洁实现：可用负梯度作为退化情形）
d = lbfgsb_direction(g);
if norm(d) < 1e-12
    return;
end

% 一次局部线搜索得到 raw candidate
xbar_n = line_search_position(n, x_n, d, X_cur, params, scene, state);

% 投影到可行域 chi_n（必须调用 Constraint_Checker）
xhat_n = Constraint_Checker('project_position', params, xbar_n);

% 非下降接受准则（真实sum rate）
f_new = local_objective_fn(n, xhat_n, X_cur, params, scene, state);
if f_new >= f_old + params.eps_X
    X_cur = replace_waveguide_position(X_cur, n, xhat_n);
end
end

function f = local_objective_fn(n, x_n, X_ref, params, scene, state)
% 局部目标函数：f_n(x_n)=R_sum(S, X_-n, x_n, W, theta, phi)
X_tmp = replace_waveguide_position(X_ref, n, x_n);
state_tmp = state;
state_tmp.X = X_tmp;
f = Signal_model('sum_rate', params, scene, state_tmp, struct());
end

function g = numerical_gradient_position(n, x_n, X_ref, params, scene, state)
% 数值微分近似梯度（中心差分）
M = numel(x_n);
g = zeros(M,1);

for i = 1:M
    e = zeros(M,1);
    e(i) = 1;

    x_p = x_n + params.step_fd * e;
    x_m = x_n - params.step_fd * e;

    x_p = Constraint_Checker('project_position', params, x_p);
    x_m = Constraint_Checker('project_position', params, x_m);

    f_p = local_objective_fn(n, x_p, X_ref, params, scene, state);
    f_m = local_objective_fn(n, x_m, X_ref, params, scene, state);

    g(i) = (f_p - f_m) / (2*params.step_fd);
end
end

function d = lbfgsb_direction(g)
% L-BFGS-B style inverse Hessian近似方向（最简可读形式）
% 历史信息不足时退化为负梯度方向
Hinv = eye(numel(g));
d = -Hinv * g;
end

function xbar = line_search_position(n, x_n, d, X_ref, params, scene, state)
% 局部线搜索：alpha从alpha0开始，按beta回缩
alpha = params.line_search_alpha0;
f0 = local_objective_fn(n, x_n, X_ref, params, scene, state);

xbar = x_n;
for t = 1:params.line_search_max_iter
    x_try = x_n + alpha * d;
    x_try = Constraint_Checker('project_position', params, x_try);

    f_try = local_objective_fn(n, x_try, X_ref, params, scene, state);
    if f_try >= f0
        xbar = x_try;
        return;
    end
    alpha = alpha * params.line_search_beta;
end
end

function X_out = replace_waveguide_position(X_in, n, x_n)
% 用第n条波导候选位置替换整体位置矩阵中的对应行
X_out = X_in;
X_out(n,:) = x_n(:).';
end
