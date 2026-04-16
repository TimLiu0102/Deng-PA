function state = Initialization_ra(params, scene, model)
% Initialization_ra：全随机初始化，保持与现有 AO 流程兼容

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

% 2) 初始位置矩阵：每条波导随机生成并投影到可行域
state.X = zeros(N, M);
for n = 1:N
    x_row = sort(rand(1, M) * params.Dy);
    x_proj = Constraint_Checker('project_position', params, x_row(:));
    state.X(n, :) = x_proj(:).';
end

% 3) 初始角度：在可行域内随机生成
state.theta = pi/2 + (pi/2) * rand(N, M);
state.phi = 2*pi*rand(N, M) - pi;
state.phi(state.phi <= -pi) = pi;

% 4) 初始预编码 W^(0)：复随机 + 统一功率归一化
W0 = randn(N*M, K_serv) + 1j * randn(N*M, K_serv);
p = real(trace(W0 * W0'));
if p > 0
    W0 = W0 * sqrt(params.P_max / p);
end
state.W = W0;

% 5) 与 AO_S 兼容字段
state.C = 1:scene.K;
state.Emax = ones(scene.K, 1);

% 6) 其余兼容字段
state.Gpot = [];
state.matching = [];
state.y_star = [];
state.y_ref = [];

end
