function result = main()
% main：论文复现主流程（参数设置 -> 初始化 -> 外层AO -> 输出）

%% 第1部分：清空环境
clc; clear; close all;

%% 第2部分：参数设置
params = struct();

% 1) 系统规模参数
params.N = 8;
params.M = 4;
params.K = 32;
params.NRF = params.N;
params.K_max = params.N;
params.K_serv = min(params.NRF, params.K_max);

% 2) 几何参数
params.Dx = 10;
params.Dy = 10;
params.d = 3.5;
params.Delta = 0.5;


% 用户位置
params.user_x_rng = [1, 20];
params.user_y_rng = [0, 20];

% 3) 信道参数
params.lambda = 0.01;
params.n_eff = 1.6;
params.alphaW = 0.01;
params.alphaL = 0.96;
params.a = 0.5;
params.b = 0.3;
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
params.T_max = 30;
params.eps_outer = 1e-4;

% 9.5) SA 联合优化参数
params.SA_max_iter = params.T_max;

% 10) 随机种子
params.seed = 7;
rng(params.seed);

% ======================== 算法方案开关 ========================
scheme_mode = 'ao_final_w';   % 'ao_final_w' | 'w_only' | 'sa_joint'

%% 第3部分：场景生成与问题定义
scene = Channel_model('build_scene', params, [], [], []);
model = Problem_formulation(params, scene);

%% 第4部分：初始化
state = Initialization(params, scene, model);
% state = Initialization_ra(params, scene, model);

if ~isfield(state, 'swap_flag')
    state.swap_flag = false;
end

%% 第5部分：初始性能与历史量
rates0 = Signal_model('individual_rates', params, scene, state, []);
R_old = sum(rates0);

history = struct();

% 初始快照
history.X0 = state.X;
history.theta0 = state.theta;
history.phi0 = state.phi;
history.S0 = state.S;
history.rates0 = rates0;

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

%% ======================== DEBUG_X START ========================
history.DEBUG_X_cells = {};
%% ======================== DEBUG_X END ==========================

% 交换标记历史（保留原有语义：首个元素对应初始化）
history.swap_flag = false;

%% 第6部分：根据 scheme_mode 执行算法
if strcmp(scheme_mode, 'ao_final_w')
    % 外层交替优化主循环（固定顺序）
    for t = 1:params.T_max
        % 当前外层迭代编号，供 AO_S 周期触发判断
        state.t = t;
        R_after_W = R_old;

        % 2) 更新角度
        [state.theta, state.phi] = AO_angle(params, scene, model, state);
        % [state.theta, state.phi] = AO_angle_ex(params, scene, model, state);
        R_after_angle = Signal_model('sum_rate', params, scene, state, []);

        % 3) 更新位置
        [state.X, DEBUG_X_t] = AO_X(params, scene, model, state);
        history.X_update_mode = 'gradient';
        % [state.X, DEBUG_X_t] = AO_X_ex(params, scene, model, state);
        % history.X_update_mode = 'exhaustive';
        R_after_X = Signal_model('sum_rate', params, scene, state, []);

        % 4) 更新用户集合
        [state.S, state.swap_flag] = AO_S(params, scene, model, state);
        % [state.S, state.swap_flag] = AO_S_ex(params, scene, model, state);
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

        %% ======================== DEBUG_X START ========================
        % history.DEBUG_X_cells{t,1} = DEBUG_X_t;
        %% ======================== DEBUG_X END ==========================

        % 8) 外层停止判断
        % 停止条件：|R_sum^(t+1)-R_sum^(t)|<eps_outer 或达到T_max
        % 各块更新均按真实sum rate非下降接受，故目标值序列单调非减
        if abs(R_new - R_old) < params.eps_outer
            break;
        end

        % 9) 更新上一轮目标值
        R_old = R_new;
    end

    % 角度、位置和用户集合交替优化结束后，最终更新一次 W
    R_before_final_W = Signal_model('sum_rate', params, scene, state, []);

    state.W = AO_W(params, scene, model, state);

    R_after_final_W = Signal_model('sum_rate', params, scene, state, []);

    history.R_before_final_W = R_before_final_W;
    history.R_after_final_W = R_after_final_W;
    history.R_sum(end,1) = R_after_final_W;

elseif strcmp(scheme_mode, 'w_only')
    history.X_update_mode = 'none';
    history.R_before_final_W = [];
    history.R_after_final_W = [];

    for t = 1:params.T_max
        state.t = t;
        state.swap_flag = false;

        % 1) 只更新 W
        state.W = AO_W(params, scene, model, state);
        R_after_W = Signal_model('sum_rate', params, scene, state, []);

        % 2) 角度、位置和用户集合全部冻结
        R_after_angle = R_after_W;
        R_after_X = R_after_W;
        R_after_S = R_after_W;

        history.R_after_W(end+1,1) = R_after_W;
        history.R_after_angle(end+1,1) = R_after_angle;
        history.R_after_X(end+1,1) = R_after_X;
        history.R_after_S(end+1,1) = R_after_S;

        R_new = R_after_S;
        history.R_sum(end+1,1) = R_new;

        history.S_cells{t,1} = state.S;
        history.X_cells{t,1} = state.X;
        history.theta_cells{t,1} = state.theta;
        history.phi_cells{t,1} = state.phi;
        history.swap_flag(end+1,1) = state.swap_flag;

        history.DEBUG_X_cells{t,1} = [];

        if abs(R_new - R_old) < params.eps_outer
            break;
        end

        R_old = R_new;
    end

elseif strcmp(scheme_mode, 'sa_joint')
    [state_best, history_sa] = SA_joint(params, scene, model, state);

    state = state_best;
    history = history_sa;

    if ~isfield(history, 'DEBUG_X_cells')
        history.DEBUG_X_cells = {};
    end
    if ~isfield(history, 'X_update_mode')
        history.X_update_mode = 'none';
    end
    if ~isfield(history, 'R_after_W')
        history.R_after_W = [];
    end
    if ~isfield(history, 'R_after_angle')
        history.R_after_angle = [];
    end
    if ~isfield(history, 'R_after_X')
        history.R_after_X = [];
    end
    if ~isfield(history, 'R_after_S')
        history.R_after_S = [];
    end
    if ~isfield(history, 'R_before_final_W')
        history.R_before_final_W = [];
    end
    if ~isfield(history, 'R_after_final_W')
        history.R_after_final_W = [];
    end
    history.R_sum(end,1) = Signal_model('sum_rate', params, scene, state, []);

else
    error('main: unsupported scheme_mode');
end

%% 第7部分：整理输出并调用结果显示模块
result = struct();
result.state = state;
result.history = history;
result.params = params;
result.scene = scene;
result.model = model;

Print_and_Plot(params, scene, model, result);

end
