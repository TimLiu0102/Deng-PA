function W_new = AO_W(params, scene, model, state)
% AO_W：固定 (S,X,theta,phi) 后，按WMMSE更新预编码矩阵W

if nargin < 3
    model = struct(); %#ok<NASGU>
end
if isfield(state, 'run_id')
    run_id = state.run_id;
else
    run_id = 'no_run_id';
end

if isfield(state, 't')
    t_id = state.t;
else
    t_id = -1;
end

if isfield(state, 'ao_w_call_id')
    call_id = state.ao_w_call_id;
else
    call_id = -1;
end

%% 1) 构造当前服务用户复合信道矩阵 H (NM x Kc)
H = build_service_channel_matrix(params, scene, state);
[Nt, Kc] = size(H);

%% 2) 生成/继承WMMSE初值（列顺序与 state.S 一致）
W = initialize_precoder_if_needed(params, state, H);

%% 3) WMMSE 内循环：W -> u -> v -> W
state_eval = state;
state_eval.W = W;
R_prev = Signal_model('sum_rate', params, scene, state_eval, struct());
W_best = W;
R_best = R_prev;

state_in = state;
R_in = Signal_model('sum_rate', params, scene, state_in, struct());
pow_in = real(trace(state.W * state.W'));
pow_init = real(trace(W * W'));
fprintf('AO_W start: run=%s, t=%d, call=%d, R_in=%.10f, R_prev=%.10f, dR_init=%.10e, pow_in=%.10f, pow_init=%.10f, ||W_init-W_in||_F=%.10e\n', ...
    run_id, t_id, call_id, R_in, R_prev, R_prev - R_in, pow_in, pow_init, norm(W - state.W, 'fro'));

for it = 1:params.I_W
    % 3.1 MMSE接收器更新 u_k
    u = update_receiver_u(H, W, params.sigma2);  % Kc x 1

    % 3.2 权重更新 v_k
    v = update_weight_v(H, W, u, params.sigma2); % Kc x 1

    % 3.3 固定u,v更新W，并用mu二分满足功率约束
    W_prev_iter = W;
    W = update_precoder_given_uv(H, u, v, params.P_max);

    % 3.4 用真实sum rate做内循环停止判断
    state_eval.W = W;
    R_now = Signal_model('sum_rate', params, scene, state_eval, struct());
    fprintf('AO_W it=%d: run=%s, t=%d, call=%d, R_now=%.10f, R_best=%.10f, dR_now_prev=%.10e, ||W-W_prev||_F=%.10e\n', ...
        it, run_id, t_id, call_id, R_now, R_best, R_now - R_prev, norm(W - W_prev_iter, 'fro'));
    if R_now > R_best
        W_best = W;
        R_best = R_now;
        fprintf('AO_W it=%d: run=%s, t=%d, call=%d, accept as new best\n', ...
            it, run_id, t_id, call_id);
    end
    if abs(R_now - R_prev) < params.eps_W
        break;
    end
    R_prev = R_now;
end

state_out = state;
state_out.W = W_best;
R_out = Signal_model('sum_rate', params, scene, state_out, struct());
fprintf('AO_W end: run=%s, t=%d, call=%d, R_out=%.10f, R_best=%.10f, R_out-R_in=%.10e, ||W_best-W_in||_F=%.10e\n', ...
    run_id, t_id, call_id, R_out, R_best, R_out - R_in, norm(W_best - state.W, 'fro'));

W_new = W_best;
end

%% ======================== 内部子函数 ========================
function H = build_service_channel_matrix(params, scene, state)
% 当前服务用户复合信道，仅调用Channel_model，不重写信道公式
extra_ch = struct();
extra_ch.use_all = false;
ch_out = Channel_model('all_users', params, scene, state, extra_ch);
H = ch_out.H;
end

function W0 = initialize_precoder_if_needed(params, state, H)
% 若state.W不可用，则用简单MRT型方向作为初值并缩放到P_max
[~, Kc] = size(H);
Nt = size(H,1);

need_init = true;
if isfield(state,'W') && ~isempty(state.W)
    W_in = state.W;
    if size(W_in,1)==Nt && size(W_in,2)==Kc
        need_init = false;
        W0 = W_in;
    end
end

if need_init
    W0 = zeros(Nt, Kc);
    for k = 1:Kc
        hk = H(:,k);
        nrm = norm(hk);
        if nrm > 0
            W0(:,k) = hk / nrm;
        end
    end
end

% 统一缩放满足总功率约束
pow = real(trace(W0*W0'));
if pow > params.P_max && pow > 0
    W0 = W0 * sqrt(params.P_max / pow);
end
end

function u = update_receiver_u(H, W, sigma2)
% 对应论文：u_k = (h_k^H w_k) / (sum_j |h_k^H w_j|^2 + sigma2)
Kc = size(H,2);
u = zeros(Kc,1);
for k = 1:Kc
    hk = H(:,k);
    denom = sum(abs(hk' * W).^2) + sigma2;
    u(k) = (hk' * W(:,k)) / denom;
end
end

function v = update_weight_v(H, W, u, sigma2)
% 对应论文：e_k = |u_k|^2*(sum_j|h_k^H w_j|^2+sigma2)-2Re(u_k h_k^H w_k)+1, v_k=1/e_k
Kc = size(H,2);
v = zeros(Kc,1);
for k = 1:Kc
    hk = H(:,k);
    total_pow = sum(abs(hk' * W).^2) + sigma2;
    ek = abs(u(k))^2 * total_pow - 2*real(u(k) * hk' * W(:,k)) + 1;
    v(k) = 1 / max(real(ek), 1e-12);
end
end

function W = update_precoder_given_uv(H, u, v, Pmax)
% 对应论文：w_k = (sum_j v_j|u_j|^2 h_j h_j^H + mu I)^(-1) * v_k * conj(u_k) * h_k
Nt = size(H,1);
Kc = size(H,2);

A = zeros(Nt,Nt);
B = zeros(Nt,Kc);
for j = 1:Kc
    hj = H(:,j);
    A = A + v(j)*abs(u(j))^2*(hj*hj');
    B(:,j) = v(j)*conj(u(j))*hj;
end

W = solve_mu_bisection(A, B, Pmax);
end

function W = solve_mu_bisection(A, B, Pmax)
% 一维二分搜索mu以满足 tr(WW^H)<=Pmax
Nt = size(A,1);
I = eye(Nt);
eps_reg = 1e-10;

% 数值对称化，避免浮点误差导致病态
A = (A + A')/2;

W0 = (A + eps_reg*I) \ B;
p0 = real(trace(W0*W0'));
if p0 <= Pmax
    W = W0;
    return;
end

mu_l = 0;
mu_u = 1;
while true
    Wu = (A + (mu_u + eps_reg)*I) \ B;
    pu = real(trace(Wu*Wu'));
    if pu <= Pmax
        break;
    end
    mu_u = 2 * mu_u;
    if mu_u > 1e12
        % 极端病态时回退到大正则解，避免死循环和奇异警告
        W = (A + (mu_u + eps_reg)*I) \ B;
        return;
    end
end

for it = 1:50
    mu_m = 0.5*(mu_l + mu_u);
    Wm = (A + (mu_m + eps_reg)*I) \ B;
    pm = real(trace(Wm*Wm'));
    if pm > Pmax
        mu_l = mu_m;
    else
        mu_u = mu_m;
    end
end

W = (A + (mu_u + eps_reg)*I) \ B;
end
