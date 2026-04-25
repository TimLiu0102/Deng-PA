function [state_best, history_hg] = HG_multiuser(params, scene, model, state0)
% HG_multiuser：Hybrid Hungarian-Greedy 多用户一次性分配 + 配置生成 + WMMSE

if nargin < 3
    model = struct(); %#ok<NASGU>
end

% 1) 构造 Gopt / Yopt / Emax
[Gopt, Yopt, Emax] = build_opt_gain_matrix(params, scene);

% 2) primary matching
[primary_pairs, S, assoc] = solve_primary_matching(params, Gopt);

% 3) greedy allocation for remaining PAs
[assoc, greedy_pairs] = greedy_allocate_remaining_pa(params, Gopt, S, assoc, primary_pairs);

% 4) build X and project
[X, theta, phi] = build_state_from_assignment(params, scene, S, assoc, Yopt);

% 5) build theta/phi（已在 build_state_from_assignment 内完成）

% 6) simple phase fine-tuning
state_tmp = struct();
state_tmp.S = S;
state_tmp.X = X;
state_tmp.theta = theta;
state_tmp.phi = phi;
state_tmp.W = AO_W(params, scene, model, state_tmp);

[X_tuned, tune_info] = phase_finetune_simple(params, scene, model, state_tmp, assoc, S, Gopt);
[X, theta, phi] = build_state_from_assignment(params, scene, S, assoc, Yopt, X_tuned);

% 7) final WMMSE
state_best = struct();
state_best.S = S;
state_best.X = X;
state_best.theta = theta;
state_best.phi = phi;

state_best.C = 1:scene.K;
state_best.Emax = Emax;
state_best.Gpot = Gopt;
state_best.matching = [];
state_best.y_star = Yopt;
state_best.y_ref = [];
state_best.swap_flag = false;

state_best.W = AO_W(params, scene, model, state_best);

% 8) build history
history_hg = build_history_hg(params, scene, state0, state_best, assoc, primary_pairs, greedy_pairs, tune_info);

end

%% ======================== 内部子函数 ========================
function [Gopt, Yopt, Emax] = build_opt_gain_matrix(params, scene)
% 构造潜在最优增益矩阵与名义最优位置
N = params.N;
M = params.M;
K = scene.K;

Gopt = zeros(K, N*M);
Yopt = zeros(K, N, M);

den = (2*log(params.alphaL))^2 - params.alphaW^2;

for k = 1:K
    qk = scene.user_pos(:,k);
    xk = qk(1);
    yk = qk(2);
    zk = qk(3);

    for n = 1:N
        A_kn = (xk - scene.xW(n))^2 + (zk - params.d)^2;
        gamma_star = sqrt(A_kn * params.alphaW^2 / den);

        for m = 1:M
            y_star = yk - gamma_star;
            Yopt(k,n,m) = y_star;

            p_star = [scene.xW(n); y_star; params.d];
            d_star = norm(qk - p_star);

            h_abs = sqrt(1/M) ...
                * exp(-(params.alphaW/2) * y_star) ...
                * (params.alphaL^d_star) ...
                * (params.lambda * params.n_refr * params.v * sqrt(2*params.a*params.b)) ...
                  / (2 * d_star);

            col = (n-1)*M + m;
            Gopt(k,col) = h_abs;
        end
    end
end

Emax = max(Gopt, [], 2);
end

function [primary_pairs, S, assoc] = solve_primary_matching(params, Gopt)
% 一对一主匹配，选出 K_serv 个服务用户和主关联 PA
K_serv = params.K_serv;
NPA = params.N * params.M;
U = Gopt;

maxPairs = min([K_serv, size(U,1), size(U,2)]);

if exist('matchpairs','file') == 2
    [pairs, ~] = matchpairs(-U, -1e12, 'max');
else
    pairs = greedy_match(U, maxPairs);
end

if isempty(pairs)
    pairs = greedy_match(U, maxPairs);
end

if size(pairs,1) > maxPairs
    gain_pairs = zeros(size(pairs,1),1);
    for i = 1:size(pairs,1)
        gain_pairs(i) = Gopt(pairs(i,1), pairs(i,2));
    end
    [~, idx_sort] = sort(gain_pairs, 'descend');
    pairs = pairs(idx_sort(1:maxPairs), :);
end

if size(pairs,1) > K_serv
    gain_pairs = zeros(size(pairs,1),1);
    for i = 1:size(pairs,1)
        gain_pairs(i) = Gopt(pairs(i,1), pairs(i,2));
    end
    [~, idx_sort] = sort(gain_pairs, 'descend');
    pairs = pairs(idx_sort(1:K_serv), :);
end

primary_pairs = pairs;
S = primary_pairs(:,1).';
S = S(:).';

assoc = zeros(numel(S), NPA);
for r = 1:size(primary_pairs,1)
    k = primary_pairs(r,1);
    col = primary_pairs(r,2);
    row = find(S == k, 1);
    assoc(row, col) = 1;
end
end

function pairs = greedy_match(U, maxPairs)
% 简单贪心后备匹配
[nr, nc] = size(U);
Uwork = U;
pairs = zeros(0,2);

for t = 1:maxPairs
    [best, idx] = max(Uwork(:));
    if ~isfinite(best)
        break;
    end
    [r, c] = ind2sub([nr, nc], idx);
    pairs(end+1,:) = [r, c]; %#ok<AGROW>
    Uwork(r,:) = -inf;
    Uwork(:,c) = -inf;
end
end

function [assoc, greedy_pairs] = greedy_allocate_remaining_pa(params, Gopt, S, assoc, primary_pairs)
% 对未被主匹配占用的 PA 做增量速率贪心分配
NPA = params.N * params.M;
remaining_cols = setdiff(1:NPA, primary_pairs(:,2).');

greedy_pairs = zeros(0,2);

for ii = 1:numel(remaining_cols)
    col = remaining_cols(ii);
    best_row = 1;
    best_delta = -inf;

    for row = 1:numel(S)
        cols_now = (assoc(row,:) == 1);
        cur_power = sum(Gopt(S(row), cols_now).^2);
        new_power = cur_power + Gopt(S(row), col)^2;
        delta_rate = log2(1 + new_power/params.sigma2) - log2(1 + cur_power/params.sigma2);

        if delta_rate > best_delta
            best_delta = delta_rate;
            best_row = row;
        end
    end

    assoc(best_row, col) = 1;
    greedy_pairs(end+1,:) = [S(best_row), col]; %#ok<AGROW>
end
end

function [X, theta, phi] = build_state_from_assignment(params, scene, S, assoc, Yopt, X_fixed)
% 由 assignment 构造位置与角度
N = params.N;
M = params.M;
NPA = N*M;

if nargin < 6 || isempty(X_fixed)
    X = zeros(N, M);

    for col = 1:NPA
        n = ceil(col / M);
        m = col - (n-1)*M;

        row = find(assoc(:,col) == 1, 1);
        k = S(row);
        X(n,m) = Yopt(k,n,m);
    end

    for n = 1:N
        x_proj = Constraint_Checker('project_position', params, X(n,:).');
        X(n,:) = x_proj(:).';
    end
else
    X = X_fixed;
end

theta = pi * ones(N, M);
phi = zeros(N, M);

for col = 1:NPA
    n = ceil(col / M);
    m = col - (n-1)*M;

    row = find(assoc(:,col) == 1, 1);
    k = S(row);

    p_nm = [scene.xW(n); X(n,m); params.d];
    qk = scene.user_pos(:,k);
    v = qk - p_nm;

    [th, ph] = compute_angle_from_vector(v);
    [th, ph] = project_angle_pair(th, ph);
    theta(n,m) = th;
    phi(n,m) = ph;
end
end

function [theta, phi] = compute_angle_from_vector(v)
% 根据方向向量计算天线朝向
r = norm(v);
if r < 1e-12
    theta = pi;
    phi = 0;
    return;
end

theta = acos(v(3)/r);
phi = atan2(v(2), v(1));
end

function [theta, phi] = project_angle_pair(theta, phi)
% 投影到可行角度范围
if theta < pi/2
    theta = pi/2;
elseif theta > pi
    theta = pi;
end

while phi <= -pi
    phi = phi + 2*pi;
end
while phi > pi
    phi = phi - 2*pi;
end
end

function [X_new, tune_info] = phase_finetune_simple(params, scene, ~, state_tmp, assoc, S, Gopt)
% 简化版相位微调：对非最强关联PA做一维小范围位置扫描
X_new = state_tmp.X;
N = params.N;
M = params.M;
offsets = linspace(-params.lambda, params.lambda, 9);

state_eval = state_tmp;
state_eval.X = X_new;
[~, state_eval.theta, state_eval.phi] = build_state_from_assignment(params, scene, S, assoc, zeros(scene.K,N,M), X_new);
base_rate = Signal_model('sum_rate', params, scene, state_eval, []);

n_try = 0;
n_accept = 0;
accept_log = zeros(0,5);

for row = 1:numel(S)
    cols_user = find(assoc(row,:) == 1);
    if numel(cols_user) <= 1
        continue;
    end

    gains = Gopt(S(row), cols_user);
    [~, idx_strong] = max(gains);
    strong_col = cols_user(idx_strong);

    cols_tune = setdiff(cols_user, strong_col);

    for jj = 1:numel(cols_tune)
        col = cols_tune(jj);
        n = ceil(col / M);
        m = col - (n-1)*M;
        y_old = X_new(n,m);

        best_local_rate = base_rate;
        best_local_x = X_new;

        for oo = 1:numel(offsets)
            n_try = n_try + 1;
            X_cand = X_new;
            X_cand(n,m) = y_old + offsets(oo);
            x_proj = Constraint_Checker('project_position', params, X_cand(n,:).');
            X_cand(n,:) = x_proj(:).';

            state_cand = state_tmp;
            state_cand.X = X_cand;
            [~, theta_cand, phi_cand] = build_state_from_assignment(params, scene, S, assoc, zeros(scene.K,N,M), X_cand);
            state_cand.theta = theta_cand;
            state_cand.phi = phi_cand;
            state_cand.W = state_tmp.W;

            R_cand = Signal_model('sum_rate', params, scene, state_cand, []);

            if R_cand > best_local_rate
                best_local_rate = R_cand;
                best_local_x = X_cand;
            end
        end

        if best_local_rate > base_rate
            X_new = best_local_x;
            base_rate = best_local_rate;
            n_accept = n_accept + 1;
            accept_log(end+1,:) = [row, col, n, m, base_rate]; %#ok<AGROW>
        end
    end
end

tune_info = struct();
tune_info.offsets = offsets;
tune_info.num_try = n_try;
tune_info.num_accept = n_accept;
tune_info.accept_log = accept_log;
tune_info.base_rate_after_tune = base_rate;

end

function history_hg = build_history_hg(params, scene, state0, state_best, assoc, primary_pairs, greedy_pairs, tune_info)
% 构造与现有流程兼容的历史输出
history_hg = struct();

history_hg.X0 = state0.X;
history_hg.theta0 = state0.theta;
history_hg.phi0 = state0.phi;
history_hg.S0 = state0.S;
history_hg.rates0 = Signal_model('individual_rates', params, scene, state0, []);

R0 = Signal_model('sum_rate', params, scene, state0, []);
R_final = Signal_model('sum_rate', params, scene, state_best, []);
history_hg.R_sum = [R0; R_final];

history_hg.R_after_W = R_final;
history_hg.R_after_angle = [];
history_hg.R_after_X = [];
history_hg.R_after_S = [];
history_hg.R_before_final_W = [];
history_hg.R_after_final_W = R_final;

history_hg.S_cells = {state_best.S};
history_hg.X_cells = {state_best.X};
history_hg.theta_cells = {state_best.theta};
history_hg.phi_cells = {state_best.phi};
history_hg.DEBUG_X_cells = {};
history_hg.X_update_mode = 'hungarian_greedy';
history_hg.swap_flag = false;

history_hg.assoc = assoc;
history_hg.primary_pairs = primary_pairs;
history_hg.greedy_pairs = greedy_pairs;
history_hg.tune_info = tune_info;
end
