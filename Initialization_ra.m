function state = Initialization_ra(params, scene, model)
% Initialization_ra：随机/默认初始化，不使用论文候选池匹配初始化

if nargin < 3
    model = struct(); %#ok<NASGU>
end

% 1) 初始服务用户集合：从全体用户中随机选 K_serv 个
state = struct();
state.S = build_random_service_set(scene, params);

% 2) 初始位置矩阵：每条波导随机采样并投影到可行域
state.X = build_random_positions(params);

% 3) 初始角度：每个PA对准当前服务用户集合中的最近用户
[state.theta, state.phi] = build_user_aligned_angles(state, scene, params);

% 4) 初始预编码 W^(0)：随机复矩阵 + 统一功率缩放
state.W = build_initial_precoder(params, scene, state);

% 5) 兼容字段
state.C = [];
state.Emax = [];
state.Gpot = [];
state.matching = [];
state.y_star = [];
state.y_ref = [];

end

function S = build_random_service_set(scene, params)
% 从全体用户中随机选 K_serv 个
K_serv = params.K_serv;
S_all = 1:scene.K;
idx = randperm(numel(S_all), K_serv);
S = S_all(idx);
S = S(:).';
end

function X = build_random_positions(params)
% 每条波导独立随机生成位置并投影到可行域
N = params.N;
M = params.M;
X = zeros(N, M);

for n = 1:N
    xbar = params.Dy * rand(M,1);
    xbar = sort(xbar, 'ascend');
    x0 = Constraint_Checker('project_position', params, xbar);
    X(n,:) = x0(:).';
end
end

function [theta, phi] = build_user_aligned_angles(state, scene, params)
% 每个PA指向当前服务用户集合中与其几何距离最近的用户
N = params.N;
M = params.M;
theta = pi * ones(N, M);
phi = zeros(N, M);

for n = 1:N
    for m = 1:M
        p_nm = [scene.xW(n); state.X(n,m); params.d];

        best_dist = inf;
        k_star = state.S(1);
        for ii = 1:numel(state.S)
            k = state.S(ii);
            d = norm(scene.user_pos(:,k) - p_nm);
            if d < best_dist
                best_dist = d;
                k_star = k;
            end
        end

        qk = scene.user_pos(:,k_star);
        v = qk - p_nm;
        [theta(n,m), phi(n,m)] = compute_angle_from_vector(v);
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

function W0 = build_initial_precoder(params, scene, state)
% 显式构造 W^(0)：随机复矩阵并统一功率缩放
Nt = params.N * params.M;
Kserv = numel(state.S);
W0 = randn(Nt, Kserv) + 1j*randn(Nt, Kserv);

p = real(trace(W0*W0'));
W0 = W0 * sqrt(params.P_max / p);
end
