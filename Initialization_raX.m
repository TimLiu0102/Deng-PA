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

% 4) 基于原始匹配关系与新X重算初始角度
[state.theta, state.phi] = rebuild_initial_angles_from_matching(state, scene, params);

% 5) 基于新的(S, X, theta, phi)重算初始预编码W
state.W = build_initial_precoder(params, scene, state);

end

function [theta, phi] = rebuild_initial_angles_from_matching(state, scene, params)
% 按原始Initialization语义：匹配PA指向主关联用户，未匹配默认(theta,phi)=(pi,0)
N = params.N;
M = params.M;

theta = pi * ones(N,M);
phi = zeros(N,M);
for r = 1:size(state.matching.pair_table,1)
    ic = state.matching.pair_table(r,1);
    n = state.matching.pair_table(r,2);
    m = state.matching.pair_table(r,3);
    k = state.C(ic);

    qk = scene.user_pos(:,k);
    p_nm = [scene.xW(n); state.X(n,m); params.d];
    v = qk - p_nm;
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
