function state = Initialization_raX(params, scene, model)
% Initialization_raX：仅随机化X初始化，其他初始化逻辑沿用原始Initialization

if nargin < 3
    model = struct(); %#ok<NASGU>
end

% 1) 先执行原始初始化
state = Initialization(params, scene, model);

% 2) 读取系统规模
N = params.N;
M = params.M;

% 3) 重置初始位置矩阵：每条波导独立随机生成可行位置
span = params.Dy - (M-1)*params.Delta;
for n = 1:N
    u = sort(rand(1, M));
    base = (0:M-1) * params.Delta;
    x_row = base + span * u;
    state.X(n,:) = x_row;
end

% 4) 基于新的(S, X, theta, phi)重算初始预编码W
state.W = build_initial_precoder(params, scene, state);

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
