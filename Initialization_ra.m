function state = Initialization_ra(params, scene, model)
% Initialization_ra：弱初始化 / 去匹配初始化

if nargin < 3
    model = struct(); %#ok<NASGU>
end

% 1) 计算潜在增益相关量（保留给后续AO_S使用）
[Gpot, y_star, Emax] = build_potential_gain_matrix_ra(params, scene);

% 2) 候选池固定为全部用户
C = 1:scene.K;

% 3) 随机初始服务用户集合
S = build_random_service_set_ra(params, scene);

% 4) 参考可行位置作为初始位置
y_ref = build_reference_positions_ra(params);
X = y_ref;

% 5) 初始角度：指向当前随机服务用户集合中的最近用户
[theta, phi] = build_initial_angles_from_random_service_set_ra(S, X, scene, params);

% 输出初始化状态
state = struct();
state.S = S;
state.X = X;
state.theta = theta;
state.phi = phi;

state.C = C;
state.Emax = Emax;
state.Gpot = Gpot;
state.matching = [];
state.y_star = y_star;
state.y_ref = y_ref;

end

%% ======================== 内部子函数 ========================
function [Gpot, y_star, Emax] = build_potential_gain_matrix_ra(params, scene)
% Potential Gain Matrix + Emax（复用原初始化对应物理公式）
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

function y_ref = build_reference_positions_ra(params)
% 参考可行位置（复用原初始化参考位置逻辑）
N = params.N; M = params.M;
y_ref = zeros(N,M);
for n = 1:N
    for m = 1:M
        y_ref(n,m) = (m-1)*params.Delta + ((m-1)/(M-1))*(params.Dy - (M-1)*params.Delta);
    end
end
end

function S = build_random_service_set_ra(params, scene)
% 从全部用户中随机选取初始服务用户
S_all = 1:scene.K;
idx = randperm(numel(S_all), params.K_serv);
S = S_all(idx);
S = S(:).';
end

function [theta, phi] = build_initial_angles_from_random_service_set_ra(S, X, scene, params)
% 每个PA指向当前随机服务用户集合中的最近用户
N = params.N; M = params.M;

theta = pi * ones(N,M);
phi = zeros(N,M);

for n = 1:N
    for m = 1:M
        p_nm = [scene.xW(n); X(n,m); params.d];

        d_min = inf;
        k_near = S(1);
        for i = 1:numel(S)
            k = S(i);
            qk = scene.user_pos(:,k);
            d_now = norm(qk - p_nm);
            if d_now < d_min
                d_min = d_now;
                k_near = k;
            end
        end

        v = scene.user_pos(:,k_near) - p_nm;
        [th, ph] = compute_angle_from_vector_ra(v);
        theta(n,m) = th;
        phi(n,m) = ph;
    end
end
end

function [theta, phi] = compute_angle_from_vector_ra(v)
% 几何方向角
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