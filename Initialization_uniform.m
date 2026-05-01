function state = Initialization_uniform(params, scene, model)
% Initialization_uniform：均匀参考构型初始化（不做proposed/random初始化）

if nargin < 3
    model = struct(); %#ok<NASGU>
end

% 1) 参考位置 X_ref（与论文参考构型一致）
X_ref = build_reference_positions(params);

% 2) 参考角度（无指向性初始化）
theta = pi * ones(params.N, params.M);
phi = zeros(params.N, params.M);

% 3) 随机服务用户集合（可复现，并恢复随机状态）
rng_state = rng;
if isfield(params,'seed') && ~isempty(params.seed)
    rng(params.seed + 7001);
else
    rng(7001);
end
S = randperm(scene.K, params.K_serv);
rng(rng_state);
S = S(:).';

% 4) 简单可行预编码（不调用AO_W，不做MRT/WMMSE）
Nt = params.N * params.M;
W = zeros(Nt, params.K_serv);
for k = 1:params.K_serv
    idx = mod(k-1, Nt) + 1;
    W(idx,k) = 1;
end
p = real(trace(W*W'));
if p > 0
    W = W * sqrt(params.P_max / p);
end

% 5) 构造兼容 AO_S 所需中间量（仅记录，不用于选S）
[Gpot, y_star, Emax] = build_potential_gain_matrix(params, scene);

state = struct();
state.S = S;
state.X = X_ref;
state.theta = theta;
state.phi = phi;
state.W = W;

state.C = 1:scene.K;
state.Emax = Emax;
state.Gpot = Gpot;
state.y_star = y_star;
state.y_ref = X_ref;
state.matching = [];
state.mu0 = [];
state.assoc_user = [];
state.assoc_count = [];
state.swap_flag = false;
state.init_mode = 'uniform_reference';
end

%% ======================== 内部子函数 ========================
function X_ref = build_reference_positions(params)
N = params.N; M = params.M;
X_ref = zeros(N,M);
if M == 1
    X_ref(:,1) = 0;
else
    for n = 1:N
        for m = 1:M
            X_ref(n,m) = (m-1)*params.Delta + ((m-1)/(M-1))*(params.Dy - (M-1)*params.Delta);
        end
    end
end
end

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
