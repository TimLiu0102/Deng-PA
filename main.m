function result = main()
% main：论文复现主流程（参数设置 -> 初始化 -> 外层AO -> 输出）

%% 第1部分：清空环境
clc; clear; close all;

%% 第2部分：参数设置
params = struct();

% 1) 系统规模参数
params.N = 4;
params.M = 4;
params.K = 32;
params.NRF = 4;
params.K_max = 4;
params.K_serv = min(params.NRF, params.K_max);

% 2) 几何参数
params.Dx = 10;
params.Dy = 10;
params.d = 5;
params.Delta = 0.5;

% 3) 信道参数
params.lambda = 0.01;
params.n_eff = 1.6;
params.alphaW = 0.01;
params.alphaL = 0.96;
params.a = 10;
params.b = 6;
params.v = 1.1;
params.n_refr = 1.5;
% 对应论文自由空间传播常数公式：eta = lambda^2 / (4*pi)
params.eta = params.lambda^2 / (4*pi);
params.P_max = 1.0;
params.sigma2 = 5e-9;

% 4) 初始化参数
params.lambda_mov = 0.05;

% 5) WMMSE 参数
params.I_W = 40;
params.eps_W = 1e-4;

% 6) 角度更新参数
params.I_theta = 6;
params.Delta_theta0 = 0.08;
params.Delta_phi0 = 0.08;
params.beta_theta = 0.6;
params.beta_phi = 0.6;
params.eps_theta = 1e-5;

% 7) 位置更新参数
params.step_fd = 1e-3;
params.line_search_alpha0 = 0.5;
params.line_search_beta = 0.5;
params.line_search_max_iter = 8;
params.eps_X = 1e-5;
params.I_X = 6;
params.lbfgs_mem = 5;

% 8) 用户集更新参数
params.T_S = 2;
params.L_in = 2;
params.L_out = 4;
params.eps_S = 1e-5;
params.max_swaps = 1;

% 9) 外层停止参数
params.T_max = 50;
params.eps_outer = 1e-4;

% 10) 随机种子
params.seed = 7;
rng(params.seed);

%% 第3部分：场景生成与问题定义
scene = Channel_model('build_scene', params, [], [], []);
model = Problem_formulation(params, scene);

%% 第4部分：初始化
state = Initialization(params, scene, model);

% 初始化不负责W，这里仅保证字段存在
if ~isfield(state, 'W')
    state.W = [];
end
if ~isfield(state, 'swap_flag')
    state.swap_flag = false;
end

%% 第5部分：初始性能与历史量
R_old = Signal_model('sum_rate', params, scene, state, []);

history = struct();

% 初始快照
history.X0 = state.X;
history.theta0 = state.theta;
history.phi0 = state.phi;
history.S0 = state.S;

% 初始 sum rate
history.R_sum = R_old;

% 每轮四块后的中间 sum rate
history.R_after_W = [];
history.R_after_angle = [];
history.R_after_X = [];
history.R_after_S = [];

% 每轮变量快照
history.S_cells = {};
history.X_cells = {};
history.theta_cells = {};
history.phi_cells = {};

% 交换标记历史（保留原有语义：首个元素对应初始化）
history.swap_flag = false;

%% 第6部分：外层交替优化主循环（固定顺序）
for t = 1:params.T_max
    % 当前外层迭代编号，供 AO_S 周期触发判断
    state.t = t;

    % 1) 更新 W
    state.W = AO_W(params, scene, model, state);
    R_after_W = Signal_model('sum_rate', params, scene, state, []);

    % 2) 更新角度
    [state.theta, state.phi] = AO_angle(params, scene, model, state);
    R_after_angle = Signal_model('sum_rate', params, scene, state, []);

    % 3) 更新位置
    state.X = AO_X(params, scene, model, state);
    R_after_X = Signal_model('sum_rate', params, scene, state, []);

    % 4) 更新用户集合
    [state.S, state.swap_flag] = AO_S(params, scene, model, state);
    R_after_S = Signal_model('sum_rate', params, scene, state, []);

    % 5) 保存每轮四块更新后的中间 sum rate
    history.R_after_W(end+1,1) = R_after_W;
    history.R_after_angle(end+1,1) = R_after_angle;
    history.R_after_X(end+1,1) = R_after_X;
    history.R_after_S(end+1,1) = R_after_S;

    % 6) 该轮最终真实 sum rate
    R_new = R_after_S;
    history.R_sum(end+1,1) = R_new;

    % 7) 保存每轮变量快照
    history.S_cells{t,1} = state.S;
    history.X_cells{t,1} = state.X;
    history.theta_cells{t,1} = state.theta;
    history.phi_cells{t,1} = state.phi;
    history.swap_flag(end+1,1) = state.swap_flag;

    % 8) 外层停止判断
    % 停止条件：|R_sum^(t+1)-R_sum^(t)|<eps_outer 或达到T_max
    % 各块更新均按真实sum rate非下降接受，故目标值序列单调非减
    if abs(R_new - R_old) < params.eps_outer
        break;
    end

    % 9) 更新上一轮目标值
    R_old = R_new;
end

%% 第8部分：整理输出并调用结果显示模块
result = struct();
result.state = state;
result.history = history;
result.scene = scene;
result.model = model;

Print_and_Plot(params, scene, model, result);

end
