function state = Initialization(params, scene, model)
% Initialization：严格对应论文初始化四步，不包含W初始化

if nargin < 3
    model = struct(); %#ok<NASGU>
end

N = params.N; M = params.M; K = scene.K;

%% Step 1) Potential Gain Matrix and Candidate User Pool
[Gpot, y_star, Emax] = build_potential_gain_matrix(params, scene);
C = build_candidate_pool(Emax, params.K_serv);

%% Step 2) Global Matching and Initial User Set S^(0)
y_ref = build_reference_positions(params);
matching = solve_global_matching(C, Gpot, y_star, y_ref, params);
S = build_initial_service_set(C, matching, params.K_serv);

%% Step 3) Initial Position Matrix X^(0)
X = build_initial_positions(C, matching, y_star, y_ref, params);

%% Step 4) Initial Orientation Angles
[theta, phi] = build_initial_angles(C, matching, X, scene, params);

% 输出初始化状态与中间量
state = struct();
state.S = S;
state.X = X;
state.theta = theta;
state.phi = phi;

state.C = C;
state.Emax = Emax;
state.Gpot = Gpot;
state.matching = matching;
state.y_star = y_star;
state.y_ref = y_ref;
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
for n = 1:N
    for m = 1:M
        y_ref(n,m) = (m-1)*params.Delta + ((m-1)/(M-1))*(params.Dy - (M-1)*params.Delta);
    end
end
end

function matching = solve_global_matching(C, Gpot, y_star, y_ref, params)
% 对应论文：maximum weight matching
% 输出 matching.M_bin(ic,n,m)=1 表示候选用户C(ic)与PA(n,m) primary association

N = params.N; M = params.M;
Nc = numel(C);
NPA = N*M;

% 构造收益矩阵 U (Nc x NPA)
U = -inf(Nc, NPA);
for ic = 1:Nc
    k = C(ic);
    for n = 1:N
        for m = 1:M
            col = (n-1)*M + m;
            d_mov = abs(y_ref(n,m) - y_star(k,n,m));
            U(ic,col) = Gpot(k,col) - params.lambda_mov * d_mov;
        end
    end
end

% 目标：在“每个用户/每个PA至多匹配一次”下，恰好选择 K_serv 条匹配
maxPairs = min([params.K_serv, Nc, NPA]);
pairs = maximum_weight_matching_hungarian(U, maxPairs);

M_bin = zeros(Nc, N, M);
pair_table = zeros(size(pairs,1), 3); % [ic, n, m]
for r = 1:size(pairs,1)
    ic = pairs(r,1);
    col = pairs(r,2);
    n = ceil(col/M);
    m = col - (n-1)*M;
    M_bin(ic,n,m) = 1;
    pair_table(r,:) = [ic, n, m];
end

matching = struct();
matching.M_bin = M_bin;
matching.pairs_ic_col = pairs;
matching.pair_table = pair_table;
matching.U = U;
end

function pairs = maximum_weight_matching_hungarian(U, maxPairs)
% 对应论文初始化主关联：maximum weight matching
% 最大化 sum M_{k,n,m} U_{k,n,m} 等价于最小化 assignment cost
% 通过引入 Nc-maxPairs 个虚拟PA列（收益=0），把“恰好maxPairs条真实匹配”写成标准指派问题
[Nc, NPA] = size(U);
Ndummy = Nc - maxPairs;
U_ext = [U, zeros(Nc, Ndummy)];

% 收益最大化 -> 代价最小化：cost = c0 - U_ext
c0 = max(U_ext(:));
cost = c0 - U_ext;

% 匈牙利算法求解：每个候选用户恰好分配1列，每列至多被1个用户使用
assign_col = hungarian_rectangular(cost);

real_mask = assign_col <= NPA;
ic_idx = find(real_mask);
col_idx = assign_col(real_mask);
pairs = [ic_idx(:), col_idx(:)];
end

function assign_col = hungarian_rectangular(cost)
% 匈牙利算法（最小化版本）：输入 n x m, 要求 n<=m
% 输出 assign_col(i)=j，表示第i行分配到第j列
[n, m] = size(cost);
if n > m
    error('Initialization: Hungarian requires n<=m');
end

u = zeros(n+1,1);
v = zeros(m+1,1);
p = zeros(m+1,1);
way = zeros(m+1,1);

for i = 1:n
    p(1) = i;
    j0 = 1;
    minv = inf(m+1,1);
    used = false(m+1,1);

    while true
        used(j0) = true;
        i0 = p(j0);
        delta = inf;
        j1 = 1;

        for j = 2:(m+1)
            if ~used(j)
                cur = cost(i0, j-1) - u(i0+1) - v(j);
                if cur < minv(j)
                    minv(j) = cur;
                    way(j) = j0;
                end
                if minv(j) < delta
                    delta = minv(j);
                    j1 = j;
                end
            end
        end

        for j = 1:(m+1)
            if used(j)
                u(p(j)+1) = u(p(j)+1) + delta;
                v(j) = v(j) - delta;
            else
                minv(j) = minv(j) - delta;
            end
        end

        j0 = j1;
        if p(j0) == 0
            break;
        end
    end

    while true
        j1 = way(j0);
        p(j0) = p(j1);
        j0 = j1;
        if j0 == 1
            break;
        end
    end
end

assign_col = zeros(n,1);
for j = 2:(m+1)
    if p(j) ~= 0
        assign_col(p(j)) = j-1;
    end
end
end

function S = build_initial_service_set(C, matching, K_serv)
% 对应论文：S^(0) = {k in C: 存在匹配}
ic_used = unique(matching.pairs_ic_col(:,1));
S = C(ic_used);
S = S(1:min(K_serv, numel(S)));
S = S(:).';
end

function X = build_initial_positions(C, matching, y_star, y_ref, params)
% 对应论文：先名义位置，再投影到chi_n
N = params.N; M = params.M;
X = y_ref;

% 将匹配到的PA位置设为对应 y_star
for r = 1:size(matching.pair_table,1)
    ic = matching.pair_table(r,1);
    n = matching.pair_table(r,2);
    m = matching.pair_table(r,3);
    k = C(ic);
    X(n,m) = y_star(k,n,m);
end

% 每条波导调用Constraint_Checker投影
for n = 1:N
    xbar_n = X(n,:).';
    x0_n = Constraint_Checker('project_position', params, xbar_n);
    X(n,:) = x0_n(:).';
end
end

function [theta, phi] = build_initial_angles(C, matching, X, scene, params)
% 对应论文：匹配PA指向其主关联用户；未匹配默认(theta,phi)=(pi,0)
N = params.N; M = params.M;
theta = pi * ones(N,M);
phi = zeros(N,M);

for r = 1:size(matching.pair_table,1)
    ic = matching.pair_table(r,1);
    n = matching.pair_table(r,2);
    m = matching.pair_table(r,3);
    k = C(ic);

    qk = scene.user_pos(:,k);
    pnm0 = [scene.xW(n); X(n,m); params.d];
    v = qk - pnm0;
    [th, ph] = compute_angle_from_vector(v);
    theta(n,m) = th;
    phi(n,m) = ph;
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
