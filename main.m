function result = main()
% 论文复现主流程：参数 -> 建模 -> 初始化 -> 外层AO -> 输出

%% 第1部分：清空环境
clc; clear; close all;

%% 第2部分：参数设置
params = struct();

% 1) 系统规模参数
params.N = 4;
params.M = 4;
params.K = 32;
params.NRF = 4;
params.K_serv = 4;
params.K_max = 4;

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
params.eta = 1.0;
params.P_max = 1.0;
params.sigma2 = 1e-3;

% 4) 初始化参数
params.lambda_mov = 0.05;

% 5) WMMSE参数
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

% 8) 用户集更新参数
params.T_S = 2;
params.L_in = 2;
params.L_out = 4;
params.eps_S = 1e-5;
params.max_swaps = 1;

% 9) 外层停止参数
params.T_max = 30;
params.eps_outer = 1e-4;

% 可复现实验随机种子
params.seed = 7;
rng(params.seed);

%% 第3部分：场景与问题定义
scene = Channel_model('build_scene', params);
model = Problem_formulation(params, scene);

%% 第4部分：初始化
state = Initialization(params, scene, model);

% 明确W不在Initialization中初始化，首次由AO_W更新
if ~isfield(state, 'W')
    state.W = [];
end
if ~isfield(state, 'swap_flag')
    state.swap_flag = false;
end

%% 第5部分：初始性能与历史量
R_old = Signal_model('sum_rate', params, scene, model, state);

history = struct();
history.R_sum = R_old;
history.S = {state.S};
history.X = {state.X};
history.theta = {state.theta};
history.phi = {state.phi};
history.swap_flag = state.swap_flag;

%% 第6部分：外层交替优化循环（固定顺序）
for t = 1:params.T_max
    state.iter = t;

    % (1) 更新W
    state.W = AO_W(params, scene, model, state);

    % (2) 更新角度
    [state.theta, state.phi] = AO_angle(params, scene, model, state);

    % (3) 更新位置
    state.X = AO_X(params, scene, model, state);

    % (4) 更新用户集
    [state.S, state.swap_flag] = AO_S(params, scene, model, state);

    % (5) 计算真实sum rate
    R_new = Signal_model('sum_rate', params, scene, model, state);

    % (6) 保存历史量
    history.R_sum(end+1,1) = R_new;
    history.S{end+1,1} = state.S;
    history.X{end+1,1} = state.X;
    history.theta{end+1,1} = state.theta;
    history.phi{end+1,1} = state.phi;
    history.swap_flag(end+1,1) = state.swap_flag;

    % (7) 论文外层停止准则：sum rate收敛
    % 由于各块更新均按真实sum rate非下降接受，目标值序列单调非减
    if abs(R_new - R_old) < params.eps_outer
        break;
    end

    % (8) 更新上一轮目标值
    R_old = R_new;
end

%% 第7部分：收尾输出
result = struct();
result.state = state;
result.history = history;
result.model = model;
result.scene = scene;

Print_and_Plot(params, scene, model, result);

end
