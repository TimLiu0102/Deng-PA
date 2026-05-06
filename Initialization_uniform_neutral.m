function state = Initialization_uniform_neutral(params, scene, model)
% Initialization_uniform_neutral：中性均匀初始化（不构造移动天线辅助量）

if nargin < 3
    model = struct(); %#ok<NASGU>
end

% 1) 参考位置 X_ref
X_ref = build_reference_positions(params);

% 2) 中性角度
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

% 4) 简单可行预编码
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

state = struct();
state.S = S;
state.X = X_ref;
state.theta = theta;
state.phi = phi;
state.W = W;

state.C = 1:scene.K;
state.Emax = [];
state.Gpot = [];
state.y_star = [];
state.matching = [];
state.y_ref = X_ref;
state.mu0 = [];
state.assoc_user = [];
state.assoc_count = [];
state.swap_flag = false;
state.init_mode = 'uniform_neutral';
end

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
