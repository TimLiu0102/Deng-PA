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
% 第n条波导的位置子问题更新（数值梯度 + L-BFGS方向 + 线搜索 + 投影）
x_cur = X_cur(n,:).';
R_old = local_objective_fn(n, x_cur, X_cur, params, scene, state);

% L-BFGS有限记忆历史：s_k = x_{k+1}-x_k, y_k = g_{k+1}-g_k
s_hist = [];
y_hist = [];

for it = 1:params.I_X
    % 这里对应论文中通过数值微分近似位置变量梯度
    g_vec = numerical_gradient_position(n, x_cur, X_cur, params, scene, state);

    % two-loop recursion 计算 quasi-Newton 方向 d = -H_k g_k
    d_vec = compute_lbfgs_direction(g_vec, s_hist, y_hist);

    % 回溯线搜索：X_trial = Projection(X + alpha*d)，以真实sum rate是否提升为准
    [x_trial, R_trial] = line_search_position(n, x_cur, d_vec, X_cur, params, scene, state, R_old);

    % 对应论文中位置块停止条件：提升小于 eps_X 则停止
    if (R_trial - R_old) < params.eps_X
        break;
    end

    g_new = numerical_gradient_position(n, x_trial, X_cur, params, scene, state);

    s_k = x_trial - x_cur;
    y_k = g_new - g_vec;

    [s_hist, y_hist] = push_lbfgs_history(s_hist, y_hist, s_k, y_k, params.lbfgs_mem);

    x_cur = x_trial;
    R_old = R_trial;
end

X_cur = replace_waveguide_position(X_cur, n, x_cur);
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

function d = compute_lbfgs_direction(g, s_hist, y_hist)
% L-BFGS two-loop recursion：用有限记忆历史近似逆Hessian作用
% s_k = x_{k+1}-x_k, y_k = g_{k+1}-g_k
% 算法含义：不显式构造Hessian，而是由历史曲率信息近似 d = -H_k g_k
k = size(s_hist,2);
q = g;

alpha = zeros(k,1);
rho = zeros(k,1);
for i = k:-1:1
    rho(i) = 1 / (y_hist(:,i)' * s_hist(:,i));
    alpha(i) = rho(i) * (s_hist(:,i)' * q);
    q = q - alpha(i) * y_hist(:,i);
end

if k >= 1
    s_last = s_hist(:,k);
    y_last = y_hist(:,k);
    gamma = (s_last' * y_last) / (y_last' * y_last);
else
    gamma = 1;
end
r = gamma * q;

for i = 1:k
    beta = rho(i) * (y_hist(:,i)' * r);
    r = r + s_hist(:,i) * (alpha(i) - beta);
end

d = -r;
end

function [xbar, fbar] = line_search_position(n, x_n, d, X_ref, params, scene, state, f0)
% 局部线搜索：alpha从alpha0开始，按beta回缩；每个试探点先投影再验收
alpha = params.line_search_alpha0;

xbar = x_n;
fbar = f0;
for t = 1:params.line_search_max_iter
    x_try = x_n + alpha * d;
    norm(x_try - x_n)x_try = Constraint_Checker('project_position', params, x_try);

    f_try = local_objective_fn(n, x_try, X_ref, params, scene, state);
    if f_try > f0
        xbar = x_try;
        fbar = f_try;
        return;
    end
    alpha = alpha * params.line_search_beta;
end
end

function [s_hist, y_hist] = push_lbfgs_history(s_hist, y_hist, s_k, y_k, mem_len)
% 维护最近 mem_len 组 (s_k, y_k)
s_hist = [s_hist, s_k];
y_hist = [y_hist, y_k];

if size(s_hist,2) > mem_len
    s_hist = s_hist(:, end-mem_len+1:end);
    y_hist = y_hist(:, end-mem_len+1:end);
end
end

function X_out = replace_waveguide_position(X_in, n, x_n)
% 用第n条波导候选位置替换整体位置矩阵中的对应行
X_out = X_in;
X_out(n,:) = x_n(:).';
end
