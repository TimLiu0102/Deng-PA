function [X_new, DEBUG_X_info] = AO_X(params, scene, model, state)
% AO_X：固定 (S,W,theta,phi) 后，按波导逐条更新位置 X

if nargin < 3
    model = struct(); %#ok<NASGU>
end

X_new = state.X;
N = size(X_new,1);

%% ======================== DEBUG_X START ========================
DEBUG_X_info = struct();
DEBUG_X_info.DEBUG_X_waveguides = cell(N,1);
%% ======================== DEBUG_X END ==========================

% 逐条波导更新，不并行
for n = 1:N
    [X_new, DEBUG_X_waveguide] = update_single_waveguide(n, X_new, params, scene, state);

    %% ======================== DEBUG_X START ========================
    DEBUG_X_info.DEBUG_X_waveguides{n,1} = DEBUG_X_waveguide;
    %% ======================== DEBUG_X END ==========================
end

end

%% ======================== 内部子函数 ========================
function [X_cur, DEBUG_X_waveguide] = update_single_waveguide(n, X_cur, params, scene, state)
% 第n条波导的位置子问题更新（数值梯度 + L-BFGS方向 + 线搜索 + 投影）
x_cur = X_cur(n,:).';
[R_old, feasible_old] = local_objective_fn(n, x_cur, X_cur, params, scene, state);
if ~feasible_old
    R_old = local_reff_value(n, x_cur, X_cur, params, scene, state);
end

%% ======================== DEBUG_X START ========================
DEBUG_X_waveguide = struct();
DEBUG_X_waveguide.DEBUG_X_iters = cell(params.I_X,1);
DEBUG_X_waveguide.DEBUG_X_g_vec = cell(params.I_X,1);
DEBUG_X_waveguide.DEBUG_X_d_vec = cell(params.I_X,1);
DEBUG_X_waveguide.DEBUG_X_d_norm_raw = nan(params.I_X,1);
DEBUG_X_waveguide.DEBUG_X_R_old = nan(params.I_X,1);
DEBUG_X_waveguide.DEBUG_X_R_trial = nan(params.I_X,1);
%% ======================== DEBUG_X END ==========================

% L-BFGS有限记忆历史：s_k = x_{k+1}-x_k, y_k = g_{k+1}-g_k
s_hist = [];
y_hist = [];

for it = 1:params.I_X
    % 这里对应论文中通过数值微分近似位置变量梯度
    g_vec = numerical_gradient_position(n, x_cur, X_cur, params, scene, state);

    % two-loop recursion 计算 quasi-Newton ascent 方向 d = H_k g_k
    d_vec = compute_lbfgs_direction(g_vec, s_hist, y_hist);

    % ======================== DIRECTION SCALE CONTROL START ========================
    % d_vec 的原始尺度来自数值梯度/L-BFGS，可能非常大。
    % 这里仅保留方向，将其归一化，使 line_search_alpha0 近似表示本轮最大移动距离。
    d_norm_raw = norm(d_vec);
    if d_norm_raw > 0
        d_vec = d_vec / d_norm_raw;
    end
    % ======================== DIRECTION SCALE CONTROL END ==========================

    % 回溯线搜索：X_trial = Projection(X + alpha*d)，以真实R_eff是否提升为准
    [x_trial, R_trial, DEBUG_X_linesearch] = line_search_position(n, x_cur, d_vec, X_cur, params, scene, state, R_old);

    %% ======================== DEBUG_X START ========================
    DEBUG_X_waveguide.DEBUG_X_iters{it,1} = DEBUG_X_linesearch;
    DEBUG_X_waveguide.DEBUG_X_g_vec{it,1} = g_vec;
    DEBUG_X_waveguide.DEBUG_X_d_vec{it,1} = d_vec;
    DEBUG_X_waveguide.DEBUG_X_d_norm_raw(it,1) = d_norm_raw;
    DEBUG_X_waveguide.DEBUG_X_R_old(it,1) = R_old;
    DEBUG_X_waveguide.DEBUG_X_R_trial(it,1) = R_trial;
    %% ======================== DEBUG_X END ==========================

    % 对应论文中位置块停止条件：有效速率提升小于 eps_X 则停止
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

function [f, feasible] = local_objective_fn(n, x_n, X_ref, params, scene, state)
% 局部目标函数：f_n(x_n)=R_eff(S, X_-n, x_n, W, theta, phi)
X_tmp = replace_waveguide_position(X_ref, n, x_n);
state_tmp = state;
state_tmp.X = X_tmp;
[R_eff, detail] = Effective_rate_model(params, scene, state_tmp, []);
feasible = detail.time_feasible;
if feasible
    f = R_eff;
else
    f = -inf;
end
end

function f = local_reff_value(n, x_n, X_ref, params, scene, state)
% 仅用于当前点基准值：不做可行性过滤，直接返回R_eff
X_tmp = replace_waveguide_position(X_ref, n, x_n);
state_tmp = state;
state_tmp.X = X_tmp;
[f, ~] = Effective_rate_model(params, scene, state_tmp, []);
end

function g = numerical_gradient_position(n, x_n, X_ref, params, scene, state)
% 数值微分近似梯度（中心差分）
M = numel(x_n);
g = zeros(M,1);
[f0, feasible0] = local_objective_fn(n, x_n, X_ref, params, scene, state);
if ~feasible0
    f0 = local_reff_value(n, x_n, X_ref, params, scene, state);
end

for i = 1:M
    e = zeros(M,1);
    e(i) = 1;

    x_p = x_n + params.step_fd * e;
    x_m = x_n - params.step_fd * e;

    x_p = Constraint_Checker('project_position', params, x_p);
    x_m = Constraint_Checker('project_position', params, x_m);

    [f_p, feasible_p] = local_objective_fn(n, x_p, X_ref, params, scene, state);
    [f_m, feasible_m] = local_objective_fn(n, x_m, X_ref, params, scene, state);

    if feasible_p && feasible_m
        g(i) = (f_p - f_m) / (2*params.step_fd);
    elseif feasible_p
        g(i) = (f_p - f0) / params.step_fd;
    elseif feasible_m
        g(i) = (f0 - f_m) / params.step_fd;
    else
        g(i) = 0;
    end
end
end

function d = compute_lbfgs_direction(g, s_hist, y_hist)
% L-BFGS two-loop recursion：用有限记忆历史近似逆Hessian作用
% s_k = x_{k+1}-x_k, y_k = g_{k+1}-g_k
% 算法含义：不显式构造Hessian，而是由历史曲率信息近似 d = H_k g_k（ascent direction）
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

d = r;
end

function [xbar, fbar, DEBUG_X_linesearch] = line_search_position(n, x_n, d, X_ref, params, scene, state, f0)
% 局部线搜索：alpha从alpha0开始，按beta回缩；每个试探点先投影再验收
alpha = params.line_search_alpha0;
if ~isfinite(f0)
    f0 = local_reff_value(n, x_n, X_ref, params, scene, state);
end

xbar = x_n;
fbar = f0;

%% ======================== DEBUG_X START ========================
L = params.line_search_max_iter;
M = numel(x_n);
DEBUG_X_linesearch = struct();
DEBUG_X_linesearch.DEBUG_X_alpha = nan(L,1);
DEBUG_X_linesearch.DEBUG_X_x_raw = nan(M,L);
DEBUG_X_linesearch.DEBUG_X_x_proj = nan(M,L);
DEBUG_X_linesearch.DEBUG_X_f_try = nan(L,1);
DEBUG_X_linesearch.DEBUG_X_delta_f = nan(L,1);
DEBUG_X_linesearch.DEBUG_X_accepted_step = nan(L,1);
DEBUG_X_linesearch.DEBUG_X_x_start = x_n;
DEBUG_X_linesearch.DEBUG_X_f0 = f0;
%% ======================== DEBUG_X END ==========================

for t = 1:params.line_search_max_iter
    x_try_raw = x_n + alpha * d;
    x_try = Constraint_Checker('project_position', params, x_try_raw);

    [f_try, feasible_try] = local_objective_fn(n, x_try, X_ref, params, scene, state);

    %% ======================== DEBUG_X START ========================
    DEBUG_X_linesearch.DEBUG_X_alpha(t,1) = alpha;
    DEBUG_X_linesearch.DEBUG_X_x_raw(:,t) = x_try_raw;
    DEBUG_X_linesearch.DEBUG_X_x_proj(:,t) = x_try;
    DEBUG_X_linesearch.DEBUG_X_f_try(t,1) = f_try;
    DEBUG_X_linesearch.DEBUG_X_delta_f(t,1) = f_try - f0;
    DEBUG_X_linesearch.DEBUG_X_accepted_step(t,1) = double(feasible_try && (f_try >= f0 + params.eps_X));
    %% ======================== DEBUG_X END ==========================

    if feasible_try && (f_try >= f0 + params.eps_X)
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
