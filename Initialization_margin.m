function state = Initialization_margin(params, scene, model)
% Initialization_margin：论文版初始化 + 未匹配 PA 的边际速率增益分配

if nargin < 3
    model = struct(); %#ok<NASGU>
end

N = params.N; M = params.M;

%% Step 1) Potential Gain Matrix and Candidate User Pool
[Gpot, y_star, Emax] = build_potential_gain_matrix(params, scene);
C = build_candidate_pool(Emax, params);

%% Step 2) Global Matching and Initial User Set S^(0)
y_ref = build_reference_positions(params);
matching = solve_global_matching(C, Gpot, y_star, y_ref, params);
S = build_initial_service_set(C, matching, params.K_serv);

%% Step 3) Primary association -> initial assoc_user / X / angles
assoc_user = zeros(N, M);
for r = 1:size(matching.pair_table,1)
    ic = matching.pair_table(r,1);
    n = matching.pair_table(r,2);
    m = matching.pair_table(r,3);
    k = C(ic);
    assoc_user(n,m) = k;
end

X = build_positions_from_assoc(assoc_user, y_star, y_ref, params);
[theta, phi] = build_angles_from_assoc(assoc_user, X, scene, params);

%% Step 4) Build current state by MRT and evaluate current rate
state_cur = struct();
state_cur.S = S;
state_cur.X = X;
state_cur.theta = theta;
state_cur.phi = phi;
state_cur.W = build_initial_precoder(params, scene, state_cur);

R_cur = Signal_model('sum_rate', params, scene, state_cur, []);

%% Step 5) Marginal rate gain assignment for remaining PAs
primary_cols = zeros(1, size(matching.pair_table,1));
for r = 1:size(matching.pair_table,1)
    n = matching.pair_table(r,2);
    m = matching.pair_table(r,3);
    primary_cols(r) = (n-1)*M + m;
end
remaining_cols = setdiff(1:N*M, primary_cols);

eps_margin = 0;
if isfield(params, 'eps_init_margin')
    eps_margin = params.eps_init_margin;
end

margin_log = zeros(0,6); % [col, n, m, k_best, best_delta, R_after]
for ii = 1:numel(remaining_cols)
    col = remaining_cols(ii);
    n = ceil(col/M);
    m = col - (n-1)*M;

    best_delta = -inf;
    best_assoc = [];
    best_state = [];
    best_R = R_cur;
    k_best = 0;

    for kk = 1:numel(S)
        k = S(kk);

        assoc_try = assoc_user;
        assoc_try(n,m) = k;

        X_try = state_cur.X;
        X_try(n,m) = y_star(k,n,m);
        x_proj = Constraint_Checker('project_position', params, X_try(n,:).');
        X_try(n,:) = x_proj(:).';

        [theta_try, phi_try] = build_angles_from_assoc(assoc_try, X_try, scene, params);

        state_try = struct();
        state_try.S = S;
        state_try.X = X_try;
        state_try.theta = theta_try;
        state_try.phi = phi_try;
        state_try.W = build_initial_precoder(params, scene, state_try);

        R_try = Signal_model('sum_rate', params, scene, state_try, []);
        delta_R = R_try - R_cur;

        if delta_R > best_delta
            best_delta = delta_R;
            best_assoc = assoc_try;
            best_state = state_try;
            best_R = R_try;
            k_best = k;
        end
    end

    if best_delta >= eps_margin
        assoc_user = best_assoc;
        state_cur = best_state;
        R_cur = best_R;
        margin_log(end+1,:) = [col, n, m, k_best, best_delta, R_cur]; %#ok<AGROW>
    end
end

%% Step 6) Output
state = state_cur;
state.C = C;
state.Emax = Emax;
state.Gpot = Gpot;
state.matching = matching;
state.y_star = y_star;
state.y_ref = y_ref;
state.assoc_user = assoc_user;
state.init_mode = 'margin_rate';
state.margin_log = margin_log;
end

%% ======================== 内部子函数 ========================
function [Gpot, y_star, Emax] = build_potential_gain_matrix(params, scene)
N = params.N; M = params.M; K = scene.K;

Gpot = zeros(K, N*M);
y_star = zeros(K, N, M);

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

Emax = max(Gpot, [], 2);
end

function C = build_candidate_pool(Emax, params)
% 按 Emax 降序选候选池 C
% 原论文默认 C_factor=2，即 |C|=2*K_serv
K = numel(Emax);

C_factor = 2;
if isfield(params, 'C_factor')
    C_factor = params.C_factor;
end

Nc = min(C_factor * params.K_serv, K);
[~, idx] = sort(Emax, 'descend');
C = idx(1:Nc);
C = C(:).';
end

function y_ref = build_reference_positions(params)
N = params.N; M = params.M;
y_ref = zeros(N,M);
for n = 1:N
    for m = 1:M
        y_ref(n,m) = (m-1)*params.Delta + ((m-1)/(M-1))*(params.Dy - (M-1)*params.Delta);
    end
end
end

function matching = solve_global_matching(C, Gpot, y_star, y_ref, params)
N = params.N; M = params.M;
Nc = numel(C);
NPA = N*M;

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

maxPairs = min([params.K_serv, Nc, NPA]);

if exist('matchpairs','file') == 2
    [pairs,~] = matchpairs(U, -1e12, 'max');
    if size(pairs,1) > maxPairs
        gain_pairs = zeros(size(pairs,1),1);
        for i = 1:size(pairs,1)
            gain_pairs(i) = U(pairs(i,1), pairs(i,2));
        end
        [~, idx_sort] = sort(gain_pairs, 'descend');
        pairs = pairs(idx_sort(1:maxPairs), :);
    end
else
    pairs = greedy_match(U, maxPairs);
end

if isempty(pairs)
    pairs = greedy_match(U, maxPairs);
end

M_bin = zeros(Nc, N, M);
pair_table = zeros(size(pairs,1), 3);
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

function pairs = greedy_match(U, maxPairs)
[nr,nc] = size(U);
Uwork = U;
pairs = zeros(0,2);
for t = 1:maxPairs
    [best,id] = max(Uwork(:));
    if ~isfinite(best), break; end
    [r,c] = ind2sub([nr,nc], id);
    pairs(end+1,:) = [r,c]; %#ok<AGROW>
    Uwork(r,:) = -inf;
    Uwork(:,c) = -inf;
end
end

function S = build_initial_service_set(C, matching, K_serv)
ic_used = unique(matching.pairs_ic_col(:,1), 'stable');
S = C(ic_used);

if numel(S) < K_serv
    rest = setdiff(C, S, 'stable');
    S = [S, rest(1:min(K_serv-numel(S), numel(rest)))];
end
S = S(1:min(K_serv, numel(S)));
S = S(:).';
end

function X = build_positions_from_assoc(assoc_user, y_star, y_ref, params)
N = params.N;
X = y_ref;

for n = 1:params.N
    for m = 1:params.M
        k = assoc_user(n,m);
        if k > 0
            X(n,m) = y_star(k,n,m);
        end
    end
end

for n = 1:N
    xbar_n = X(n,:).';
    x0_n = Constraint_Checker('project_position', params, xbar_n);
    X(n,:) = x0_n(:).';
end
end

function [theta, phi] = build_angles_from_assoc(assoc_user, X, scene, params)
N = params.N;
M = params.M;
theta = pi * ones(N,M);
phi = zeros(N,M);

for n = 1:N
    for m = 1:M
        k = assoc_user(n,m);
        if k > 0
            qk = scene.user_pos(:,k);
            pnm = [scene.xW(n); X(n,m); params.d];
            v = qk - pnm;
            [th, ph] = compute_angle_from_vector(v);
            theta(n,m) = th;
            phi(n,m) = ph;
        end
    end
end
end

function [theta, phi] = compute_angle_from_vector(v)
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

function W0 = build_initial_precoder(params, scene, state)
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
