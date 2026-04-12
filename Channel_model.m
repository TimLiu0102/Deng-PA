function out = Channel_model(mode, params, scene, state, extra)
% Channel Model：仅包含几何与信道建模

if nargin < 3 || isempty(scene), scene = struct(); end
if nargin < 4 || isempty(state), state = struct(); end
if nargin < 5 || isempty(extra), extra = struct(); end

switch lower(mode)
    case 'build_scene'
        out = build_scene(params);

    case 'guided_wave'
        % 输入: extra.x_n (M x 1)
        out = guided_wave(extra.x_n, params);

    case 'radiation_pattern'
        % 输入: extra.local_xyz = [x;y;z]
        xyz = extra.local_xyz(:);
        out = radiation_pattern(xyz(1), xyz(2), xyz(3), params);

    case 'global_to_local'
        % 输入: extra.q_k, extra.p_nm, extra.theta_nm, extra.phi_nm
        out = global_to_local(extra.q_k, extra.p_nm, extra.theta_nm, extra.phi_nm);

    case 'free_space'
        % 输入: extra.k, extra.n, 依赖 state.X/theta/phi
        out = free_space_channel(extra.k, extra.n, params, scene, state);

    case 'composite_user'
        % 输入: extra.k, 依赖 state.X/theta/phi
        out = composite_user_channel(extra.k, params, scene, state);

    case 'all_users'
        % extra.use_all=true 时计算1:K，否则优先用 state.S
        out = all_users_channel(params, scene, state, extra);

    otherwise
        error('Channel_model: unsupported mode');
end

end

%% ======================== 内部子函数 ========================
function scene = build_scene(params)
% 对应论文几何定义：用户位置、波导横向坐标、馈电点
N = params.N; M = params.M; K = params.K;

xW = ((2*(1:N)-1)/(2*N)) * params.Dx;      % 1 x N
feed_pos = [xW; zeros(1,N); params.d*ones(1,N)]; % 3 x N

if isfield(params,'user_x_rng'), xr = params.user_x_rng; else, xr = [0, params.Dx]; end
if isfield(params,'user_y_rng'), yr = params.user_y_rng; else, yr = [0, params.Dy]; end
xk = xr(1) + (xr(2)-xr(1))*rand(1,K);
yk = yr(1) + (yr(2)-yr(1))*rand(1,K);
zk = zeros(1,K);

scene = struct();
scene.user_pos = [xk; yk; zk];    % 3 x K
scene.xW = xW;                    % 1 x N
scene.feed_pos = feed_pos;        % 3 x N
scene.waveguide_axis = 'y';
scene.K = K;
scene.N = N;
scene.M = M;
end

function g_n = guided_wave(x_n, params)
% 对应论文导波信道 g_n(x_n)
lambda_g = params.lambda / params.n_eff;
M = params.M;
y = x_n(:);
g_n = sqrt(1/M) * exp(-(params.alphaW/2)*y) .* exp(-1j*2*pi/lambda_g*y);
end

function ups = radiation_pattern(x, y, z, params)
% 对应论文局部方向图 Upsilon(x,y,z)
% y_eff = max(y, 1e-9);  % 仅最基本数值保护，不改变公式结构
eps0 = 1e-9;  % 测试用平滑 y_eff 保护
y_eff = sqrt(y.^2 + eps0^2);
k0 = 2*pi/params.lambda;

w1 = params.v * params.a * params.lambda;
w2 = params.v * params.b * params.lambda;
B = sqrt(2/(pi*w1*w2)); % 由 ∫∫|f(x,y,z)|^2 dx dz = 1 得到
W1 = params.lambda * y_eff / (pi * params.n_refr * w1);
W2 = params.lambda * y_eff / (pi * params.n_refr * w2);
R1 = y_eff;
R2 = y_eff;
Theta1 = atan(params.lambda * y_eff / (pi * params.n_refr * w1));
Theta2 = atan(params.lambda * y_eff / (pi * params.n_refr * w2));

% ups = sqrt((w1*w2)/(W1*W2)) * B * exp( ...
%     -(x^2/W1^2 + z^2/W2^2) ...
%     -1j * k0 * params.n_refr * (x^2/(2*R1) + z^2/(2*R2) + y_eff) ...
%     +1j * (Theta1 + Theta2)/2 );
% 测试版本：只保留 y 项、去掉 x,z 横向笔形项
ups = sqrt((w1*w2)/(W1*W2)) * B * exp( ...
    -1j * k0 * params.n_refr * y_eff ...
    +1j * (Theta1 + Theta2)/2 );
end

function local_xyz = global_to_local(q_k, p_nm, theta_nm, phi_nm)
% 对应论文全局到局部坐标变换
v = q_k(:) - p_nm(:);
Rz = [cos(pi/2 - phi_nm), -sin(pi/2 - phi_nm), 0; ...
      sin(pi/2 - phi_nm),  cos(pi/2 - phi_nm), 0; ...
      0,                   0,                  1];
Rx = [1, 0, 0; ...
      0, cos(theta_nm - pi/2), -sin(theta_nm - pi/2); ...
      0, sin(theta_nm - pi/2),  cos(theta_nm - pi/2)];
local_xyz = Rx * Rz * v;
end

function h_tilde_kn = free_space_channel(k, n, params, scene, state)
% 对应论文自由空间信道 h_tilde_{k,n}(x_n)
M = params.M;
q_k = scene.user_pos(:,k);
x_n = state.X(n,:).';
theta_n = state.theta(n,:).';
phi_n = state.phi(n,:).';

h_tilde_kn = zeros(M,1);
for m = 1:M
    p_nm = [scene.xW(n); x_n(m); params.d];
    d_k_nm = norm(q_k - p_nm);

    local_xyz = global_to_local(q_k, p_nm, theta_n(m), phi_n(m));
    ups = radiation_pattern(local_xyz(1), local_xyz(2), local_xyz(3), params);

    % 严格对应论文自由空间信道公式：sqrt(eta) * alphaL^d * Upsilon
    h_tilde_kn(m) = sqrt(params.eta) * (params.alphaL^d_k_nm) * ups;
end
end

function h_k = composite_user_channel(k, params, scene, state)
% 对应论文总复合信道 h_k(X)
N = params.N; M = params.M;
h_k = zeros(N*M,1);

for n = 1:N
    x_n = state.X(n,:).';
    g_n = guided_wave(x_n, params);
    h_tilde_kn = free_space_channel(k, n, params, scene, state);

    idx = (n-1)*M + (1:M);
    h_k(idx) = g_n .* h_tilde_kn;
end
end

function out = all_users_channel(params, scene, state, extra)
% 计算全部用户或服务用户集合的复合信道矩阵
if isfield(extra, 'use_all') && extra.use_all
    user_idx = 1:scene.K;
elseif isfield(state, 'S') && ~isempty(state.S)
    user_idx = state.S(:).';
else
    user_idx = 1:scene.K;
end

NM = params.N * params.M;
Ksel = numel(user_idx);
H = zeros(NM, Ksel);
for ii = 1:Ksel
    H(:,ii) = composite_user_channel(user_idx(ii), params, scene, state);
end

out = struct();
out.H = H;
out.user_idx = user_idx;
end
