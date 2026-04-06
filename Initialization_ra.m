function state = Initialization_ra(params, scene, model)
% Initialization_ra：随机/默认初始化（不使用论文候选池与匹配初始化）

if nargin < 3
    model = struct(); %#ok<NASGU>
end

N = params.N;
M = params.M;
K_serv = params.K_serv;

% 1) 初始用户集合：从全部用户中随机选 K_serv 个
S_all = 1:scene.K;
S_rand = randperm(numel(S_all), K_serv);
state.S = S_all(S_rand);
state.S = state.S(:).';

% 2) 初始位置矩阵：每条波导上等间隔分布
x_row = linspace(0, params.Dy, M);
state.X = repmat(x_row, N, 1);

% 3) 初始角度：统一竖直向下
state.theta = pi * ones(N, M);
state.phi = zeros(N, M);

% 4) 兼容字段
state.C = [];
state.Emax = [];
state.Gpot = [];
state.matching = [];
state.y_star = [];
state.y_ref = [];

end
