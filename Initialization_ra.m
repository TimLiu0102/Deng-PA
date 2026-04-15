function state = Initialization_ra(params, scene, model)
% Initialization_ra：随机/默认初始化，不使用论文候选池匹配初始化

if nargin < 3
    model = struct(); %#ok<NASGU>
end

N = params.N;
M = params.M;
K_serv = params.K_serv;

% 1) 初始服务用户集合：从全体用户中随机选 K_serv 个
S_all = 1:scene.K;
S_rand = randperm(numel(S_all), K_serv);
state.S = S_all(S_rand);
state.S = state.S(:).';

% 2) 初始位置矩阵：每条波导独立随机生成可行位置
span = params.Dy - (M-1)*params.Delta;
state.X = zeros(N, M);
for n = 1:N
    u = sort(rand(1, M));
    base = (0:M-1) * params.Delta;
    x_row = base + span * u;
    state.X(n,:) = x_row;
end

% 3) 初始角度：统一竖直朝下
state.theta = pi * ones(N, M);
state.phi = zeros(N, M);

% 4) 初始预编码 W^(0)：MRT + 统一功率缩放
state.W = build_initial_precoder(params, scene, state);

% 5) 兼容字段
state.C = [];
state.Emax = [];
state.Gpot = [];
state.matching = [];
state.y_star = [];
state.y_ref = [];

end

function W0 = build_initial_precoder(params, scene, state)
% 显式构造 W^(0)：按当前服务用户复合信道做 MRT
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
