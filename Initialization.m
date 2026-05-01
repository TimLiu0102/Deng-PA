function state = Initialization(params, scene, model)
% Initialization：论文初始化：候选池 -> 初始服务用户 -> 位移折扣效用 ->
% PA关联 -> 位置投影 -> 角度 -> MRT预编码

if nargin < 3
    model = struct(); %#ok<NASGU>
end

N = params.N; M = params.M; K = scene.K;

%% Step 1) Potential Gain Matrix and Candidate User Pool
[Gpot, y_star, Emax] = build_potential_gain_matrix(params, scene);
C = build_candidate_pool(Emax, params.K_serv);

%% Step 2) Initial User Set S^(0)
S = build_initial_service_set_by_emax(C, Emax, params.K_serv);

%% Step 3) Reference Positions
y_ref = build_reference_positions(params);

%% Step 4) Displacement-discounted Utility and PA Association
[mu0, X_bar] = build_pa_association_and_nominal_position(S, Gpot, y_star, y_ref, params);

%% Step 5) Initial Position Projection
X = build_initial_positions_from_nominal(X_bar, params);

%% Step 6) Initial Orientation Angles
[theta, phi] = build_initial_angles_from_assoc(mu0, X, scene, params);

%% Step 7) Initial Precoders W^(0) by MRT
tmp_state = struct();
tmp_state.S = S;
tmp_state.X = X;
tmp_state.theta = theta;
tmp_state.phi = phi;
W = build_initial_precoder(params, scene, tmp_state);

% 输出初始化状态与中间量
state = struct();
state.S = S;
state.X = X;
state.theta = theta;
state.phi = phi;
state.W = W;

state.C = C;
state.Emax = Emax;
state.Gpot = Gpot;
state.matching = [];
state.y_star = y_star;
state.y_ref = y_ref;
state.mu0 = mu0;
state.assoc_user = mu0;
assoc_count = zeros(1, numel(S));
for i = 1:numel(S)
    assoc_count(i) = sum(mu0(:) == S(i));
end
state.assoc_count = assoc_count;
state.init_mode = 'paper_reff';
end

%% ======================== 内部子函数 ========================
function [Gpot, y_star, Emax] = build_potential_gain_matrix(params, scene)
% 对应论文：Potential Gain Matrix + Emax
N = params.N; M = params.M; K = scene.K;

Gpot = zeros(K, N*M);      % K x (N*M), 元素为 |h_k,n,m^(0)|
y_star = zeros(K, N, M);   % y_star(k,n,m)

for k = 1:K
    qk = scene.user_pos(:,k);
    xk = qk(1); yk = qk(2); zk = qk(3);

    for n = 1:N
        A_kn = (xk - scene.xW(n))^2 + (zk - params.d)^2;
        den = (2*log(params.alphaL))^2 - params.alphaW^2;
        gamma_star = sqrt(A_kn * params.alphaW^2 / den);

        for m = 1:M
            y_star(k,n,m) = yk - gamma_star;
            p_star = [scene.xW(n); y_star(k,n,m); params.d];
            d_star = norm(qk - p_star);

            h_abs = sqrt(1/M) ...
                * exp(-(params.alphaW/2) * y_star(k,n,m)) ...
                * (params.alphaL^d_star) ...
                * (params.lambda * params.n_refr * params.v * sqrt(2*params.a*params.b)) ...
                  / (2 * d_star);

            col = (n-1)*M + m;
            Gpot(k,col) = h_abs;
        end
    end
end

Emax = max(Gpot, [], 2); % Kx1
end

function C = build_candidate_pool(Emax, K_serv)
% 对应论文：按Emax降序取前2*K_serv作为候选池C
K = numel(Emax);
[~, idx] = sort(Emax, 'descend');
C = idx(1:min(2*K_serv, K));
C = C(:).';
end

function y_ref = build_reference_positions(params)
% 对应论文参考位置 y_ref_{n,m}
N = params.N; M = params.M;
y_ref = zeros(N,M);
if M == 1
    y_ref(:,1) = 0;
else
    for n = 1:N
        for m = 1:M
            y_ref(n,m) = (m-1)*params.Delta + ((m-1)/(M-1))*(params.Dy - (M-1)*params.Delta);
        end
    end
end
end

function S = build_initial_service_set_by_emax(C, Emax, K_serv)
% S^(0)：候选池C中按Emax降序取前K_serv
[~, ord] = sort(Emax(C), 'descend');
S = C(ord(1:min(K_serv, numel(C))));
S = S(:).';
end

function [mu0, X_bar] = build_pa_association_and_nominal_position(S, Gpot, y_star, y_ref, params)
% 两步分配：
% Step A 先覆盖所有服务用户（每个用户至少一个PA主关联）
% Step B 剩余PA再按原效用贪心分配
N = params.N; M = params.M;
NPA = N * M;
Kserv = numel(S);
mu0 = zeros(N,M);
X_bar = y_ref;

% U_serv(i,col): i为服务用户索引，col为PA展平索引
U_serv = -inf(Kserv, NPA);
for is = 1:Kserv
    k = S(is);
    for n = 1:N
        for m = 1:M
            col = (n-1)*M + m;
            D = abs(y_star(k,n,m) - y_ref(n,m)) / params.Dy;
            U_serv(is,col) = max(1 - params.rho * D, 0) * Gpot(k, col);
        end
    end
end

% Step A: 用户覆盖分配（每个用户至少一次，且PA不重复）
pairs = zeros(0,2); % 每行 [i_user, col_pa]
if exist('matchpairs', 'file') == 2
    tmp_pairs = matchpairs(U_serv, -1e12, 'max');
    if ~isempty(tmp_pairs)
        take_num = min(size(tmp_pairs,1), Kserv);
        pairs = tmp_pairs(1:take_num, :);
    end
else
    U_work = U_serv;
    for t = 1:min(Kserv, NPA)
        [mx, idx] = max(U_work(:));
        if ~isfinite(mx)
            break;
        end
        [i_user, col_pa] = ind2sub(size(U_work), idx);
        pairs(end+1,:) = [i_user, col_pa]; %#ok<AGROW>
        U_work(i_user,:) = -inf;
        U_work(:,col_pa) = -inf;
    end
end

used_pa = false(1, NPA);
covered_user = false(1, Kserv);

% 写入已找到的覆盖配对
for r = 1:size(pairs,1)
    i_user = pairs(r,1);
    col_pa = pairs(r,2);
    if i_user < 1 || i_user > Kserv || col_pa < 1 || col_pa > NPA
        continue;
    end
    if covered_user(i_user) || used_pa(col_pa)
        continue;
    end
    k = S(i_user);
    n = floor((col_pa-1)/M) + 1;
    m = mod(col_pa-1, M) + 1;
    mu0(n,m) = k;
    X_bar(n,m) = y_star(k,n,m);
    covered_user(i_user) = true;
    used_pa(col_pa) = true;
end

% 极端情况下补齐未覆盖用户
for i_user = 1:Kserv
    if covered_user(i_user)
        continue;
    end
    free_cols = find(~used_pa);
    if isempty(free_cols)
        break;
    end
    [~, loc] = max(U_serv(i_user, free_cols));
    col_pa = free_cols(loc);
    k = S(i_user);
    n = floor((col_pa-1)/M) + 1;
    m = mod(col_pa-1, M) + 1;
    mu0(n,m) = k;
    X_bar(n,m) = y_star(k,n,m);
    covered_user(i_user) = true;
    used_pa(col_pa) = true;
end

% Step B: 剩余PA按原效用贪心
for col_pa = find(~used_pa)
    [~, i_best] = max(U_serv(:, col_pa));
    k_best = S(i_best);
    n = floor((col_pa-1)/M) + 1;
    m = mod(col_pa-1, M) + 1;
    mu0(n,m) = k_best;
    X_bar(n,m) = y_star(k_best,n,m);
end
end

function X = build_initial_positions_from_nominal(X_bar, params)
% 对每条波导调用Constraint_Checker投影
N = params.N; M = params.M;
X = X_bar;
for n = 1:N
    xbar_n = X(n,:).';
    x0_n = Constraint_Checker('project_position', params, xbar_n);
    X(n,:) = x0_n(:).';
end
end

function [theta, phi] = build_initial_angles_from_assoc(mu0, X, scene, params)
% 每个PA按关联用户mu0(n,m)设置初始角度
N = params.N; M = params.M;
theta = pi * ones(N,M);
phi = zeros(N,M);

for n = 1:N
    for m = 1:M
        k = mu0(n,m);
        qk = scene.user_pos(:,k);
        pnm0 = [scene.xW(n); X(n,m); params.d];
        v = qk - pnm0;
        [th, ph] = compute_angle_from_vector(v);
        [theta(n,m), phi(n,m)] = project_angle_pair(th, ph);
    end
end
end

function [theta, phi] = compute_angle_from_vector(v)
% 几何方向角：与后续Channel_model坐标变换约定一致
r = norm(v);
if r < 1e-12
    theta = pi;
    phi = 0;
    return;
end

theta = acos(v(3)/r);
theta = min(max(theta, pi/2), pi);
phi = atan2(v(2), v(1));
if phi <= -pi
    phi = phi + 2*pi;
end
if phi > pi
    phi = phi - 2*pi;
end
end

function [theta, phi] = project_angle_pair(theta, phi)
theta = min(max(theta, pi/2), pi);
if phi <= -pi
    phi = phi + 2*pi;
end
if phi > pi
    phi = phi - 2*pi;
end
end

function W0 = build_initial_precoder(params, scene, state)
% 对应初始化显式 W^(0)：按当前服务用户复合信道做 MRT 并统一功率缩放
extra_ch = struct();
extra_ch.use_all = false;
ch_out = Channel_model('all_users', params, scene, state, extra_ch);
H = ch_out.H;

Nt = size(H,1);
Kserv = size(H,2);
W0 = zeros(Nt, Kserv);
for k = 1:Kserv
    hk = H(:,k);
    nrm = norm(hk);
    if nrm > 0
        W0(:,k) = hk / nrm;
    end
end

p = real(trace(W0*W0'));
if p > params.P_max && p > 0
    W0 = W0 * sqrt(params.P_max/p);
end
end
