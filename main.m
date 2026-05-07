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
% 服务/部署空间大小
params.area_Dx = 20;
params.area_Dy = 20;

% 波导沿 x 方向的部署宽度
params.waveguide_Dx = 20;

% PA 沿波导 y 方向的可移动长度 / 波导长度
params.waveguide_Dy = 20;

% 为兼容旧函数，保留 Dx / Dy
params.Dx = params.waveguide_Dx;
params.Dy = params.waveguide_Dy;
params.d = 3;
params.Delta = 0.5;


% 用户位置
params.user_x_rng = [0, params.area_Dx];
params.user_y_rng = [0, params.area_Dy];

% 3) 信道参数
params.lambda = 0.01;
params.n_eff = 1.6;
params.alphaW = 0.01;
params.alphaL = 0.96;
params.a = 0.3;
params.b = 0.18;
params.v = 1.1;
params.n_refr = 1.5;
% 对应论文自由空间传播常数公式：eta = lambda^2 / (4*pi)
params.eta = params.lambda^2 / (4*pi);
params.P_max = 1.0;
params.sigma2 = 5e-9;
params.SNR_dB = 20;
% 绘图参考 SNR：当前 sigma2 对应 SNR_dB

% 4) 初始化参数
params.lambda_mov = 0.05;
params.init_X_grid_num = 21;
params.init_X_tau_max = 1.0;
% params.init_X_move_max = 1.0;

% 4.5) 有效速率模型参数
params.T_f = 5;             % frame duration, s
params.v_PA = 10;            % PA moving speed, m/s
params.omega_theta = 300;   % elevation rotation speed, rad/s
params.omega_phi = 300;     % azimuth rotation speed, rad/s
params.rho = 1.5;

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
params.I_X_refine = 3;                  % 每条波导的粗到细位置内迭代轮数
params.X_refine_grid_num = 41;          % 每次一维搜索的网格点数
params.X_refine_radius_list = [inf, 0.5, 0.1, 0.02, 0.005];  % 第一轮全区间，后续局部细化
params.eps_X_refine = 0;                % 调试阶段允许任意非下降提升
params.X_refine_use_ystar = true;       % 若 state.y_star 存在，则把 y_star 加入候选点

% 8) 用户集更新参数
params.T_S = 1;
params.L_in = 2;
params.L_out = 4;
params.eps_S = 1e-5;
params.max_swaps = 1;

% 9) 外层停止参数
params.T_max = 10;
params.eps_outer = 1e-4;

% 9.5) SA 联合优化参数（纯启发式联合搜索）
params.SA_max_iter = 5000;
params.SA_T0 = 1.0;
params.SA_alpha = 0.985;
params.SA_step_X = 0.5;
params.SA_step_theta = 0.08;
params.SA_step_phi = 0.08;
params.SA_step_W = 0.02;
params.SA_beta_W = 0.3;

% 9.6) PSO 联合启发式搜索参数
params.PSO_num_particles = 40;
params.PSO_max_iter = 500;
params.PSO_w = 0.7;
params.PSO_c1 = 1.5;
params.PSO_c2 = 1.5;
params.PSO_vmax_X = 1.0;
params.PSO_vmax_theta = 0.15;
params.PSO_vmax_phi = 0.15;
params.PSO_prob_S = 0.4;
params.PSO_beta_W = 0.3;
params.PSO_step_W = 0.02;
params.PSO_init_X_jitter = 2.0;

% 10) 随机种子
params.seed = 7;
rng(params.seed);

% ======================== 算法方案开关 ========================
scheme_mode = 'pso_joint';   % 'AO' | 'sa_joint' | 'pso_joint' | 'hg_multiuser' | 'fixed_antenna_ws' | 'fixedX' | 'w_only'

%% 第3部分：场景生成与问题定义
scene = Channel_model('build_scene', params, [], [], []);
model = Problem_formulation(params, scene);

%% 第4部分：初始化
init_mode = 'uniform_neutral';   % 'paper' | 'uniform_neutral' | 'uniform_fixed' | 'fixedX' | 'reffX' | 'margin' | 'random' | 'uniform'

if strcmp(init_mode, 'paper')
    state = Initialization(params, scene, model);
elseif strcmp(init_mode, 'margin')
    state = Initialization_margin(params, scene, model);
elseif strcmp(init_mode, 'random')
    state = Initialization_ra(params, scene, model);
elseif strcmp(init_mode, 'uniform')
    state = Initialization_uniform(params, scene, model);
elseif strcmp(init_mode, 'uniform_neutral')
    state = Initialization_uniform_neutral(params, scene, model);
elseif strcmp(init_mode, 'uniform_fixed')
    state = Initialization_uniform_fixed(params, scene, model);
elseif strcmp(init_mode, 'fixedX')
    state = Initialization_fixedX(params, scene, model);
elseif strcmp(init_mode, 'reffX')
    state = Initialization_reffX(params, scene, model);
else
    error('main: unsupported init_mode');
end

if ~isfield(state, 'swap_flag')
    state.swap_flag = false;
end

%% 第5部分：初始性能与历史量
rates0 = Signal_model('individual_rates', params, scene, state, []);
R_sum0 = sum(rates0);
[R_eff0, detail0] = Effective_rate_model(params, scene, state, []);
R_old = R_eff0;

history = struct();

% 初始快照
history.X0 = state.X;
history.theta0 = state.theta;
history.phi0 = state.phi;
history.S0 = state.S;
history.rates0 = rates0;

% 初始 sum rate 与有效速率
history.R_sum = R_sum0;
history.R_eff = R_eff0;
history.T_X = detail0.T_X;
history.T_theta = detail0.T_theta;
history.T_phi = detail0.T_phi;
history.T_rec = detail0.T_rec;
history.time_factor = detail0.time_factor;

% 每轮四块后的中间 sum rate
history.R_after_W = [];
history.R_after_angle = [];
history.R_after_X = [];
history.R_after_S = [];
history.R_eff_after_W = [];
history.R_eff_after_angle = [];
history.R_eff_after_X = [];
history.R_eff_after_S = [];

% 每轮变量快照
history.S_cells = {};
history.X_cells = {};
history.theta_cells = {};
history.phi_cells = {};

%% ======================== DEBUG_X START ========================
% history.DEBUG_X_cells = {};
%% ======================== DEBUG_X END ==========================

% 交换标记历史（保留原有语义：首个元素对应初始化）
history.swap_flag = false;

%% 第6部分：根据 scheme_mode 执行算法
if strcmp(scheme_mode, 'AO')
    % 外层交替优化主循环：W -> angle -> X -> S
    for t = 1:params.T_max
        % 当前外层迭代编号，供 AO_S 周期触发判断
        state.t = t;

         state.W = AO_W(params, scene, model, state);
         R_after_W = Signal_model('sum_rate', params, scene, state, []);
         [R_eff_after_W, ~] = Effective_rate_model(params, scene, state, []);

        % 2) 更新角度
        [state.theta, state.phi] = AO_angle(params, scene, model, state);
        % [state.theta, state.phi] = AO_angle_ex(params, scene, model, state);
        R_after_angle = Signal_model('sum_rate', params, scene, state, []);
        [R_eff_after_angle, ~] = Effective_rate_model(params, scene, state, []);

        % 3) 更新位置
        [state.X, DEBUG_X_t] = AO_X(params, scene, model, state);
        history.X_update_mode = 'gradient';

        % [state.X, DEBUG_X_t] = AO_X_ex(params, scene, model, state);
        % history.X_update_mode = 'exhaustive';

        % [state.X, DEBUG_X_t] = AO_X_grid_refine(params, scene, model, state);
        % history.X_update_mode = 'grid_refine';
        R_after_X = Signal_model('sum_rate', params, scene, state, []);
        [R_eff_after_X, ~] = Effective_rate_model(params, scene, state, []);

        % 4) 更新用户集合
        % 按论文思路，S 候选交换评价时固定当前 W，不对每个候选重新 WMMSE；
        % 若发生用户交换，新 W 会在下一轮 W 子问题中重新适配。
        [state.S, state.swap_flag] = AO_S(params, scene, model, state);
        % [state.S, state.swap_flag] = AO_S_ex(params, scene, model, state);
        R_after_S = Signal_model('sum_rate', params, scene, state, []);
        [R_eff_after_S, detail_S] = Effective_rate_model(params, scene, state, []);

        % 5) 保存每轮四块更新后的中间 sum rate
        history.R_after_W(end+1,1) = R_after_W;
        history.R_after_angle(end+1,1) = R_after_angle;
        history.R_after_X(end+1,1) = R_after_X;
        history.R_after_S(end+1,1) = R_after_S;
        history.R_eff_after_W(end+1,1) = R_eff_after_W;
        history.R_eff_after_angle(end+1,1) = R_eff_after_angle;
        history.R_eff_after_X(end+1,1) = R_eff_after_X;
        history.R_eff_after_S(end+1,1) = R_eff_after_S;

        % 6) 该轮最终性能：R_eff 与原始 sum rate
        R_new = R_eff_after_S;
        history.R_eff(end+1,1) = R_new;
        history.R_sum(end+1,1) = R_after_S;
        history.T_X(end+1,1) = detail_S.T_X;
        history.T_theta(end+1,1) = detail_S.T_theta;
        history.T_phi(end+1,1) = detail_S.T_phi;
        history.T_rec(end+1,1) = detail_S.T_rec;
        history.time_factor(end+1,1) = detail_S.time_factor;

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
        % 停止条件：|R_eff^(t+1)-R_eff^(t)|<eps_outer 或达到T_max
        % R_eff 序列作为外层判断指标
        if abs(R_new - R_old) < params.eps_outer
            break;
        end

        % 9) 更新上一轮目标值
        R_old = R_new;
    end


elseif strcmp(scheme_mode, 'fixedX')
    [state, history] = AO_fixedX(params, scene, model, state);

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
        [R_eff_after_W, detail_W] = Effective_rate_model(params, scene, state, []);

        % 2) 角度、位置和用户集合全部冻结
        R_after_angle = R_after_W;
        R_after_X = R_after_W;
        R_after_S = R_after_W;
        R_eff_after_angle = R_eff_after_W;
        R_eff_after_X = R_eff_after_W;
        R_eff_after_S = R_eff_after_W;

        history.R_after_W(end+1,1) = R_after_W;
        history.R_after_angle(end+1,1) = R_after_angle;
        history.R_after_X(end+1,1) = R_after_X;
        history.R_after_S(end+1,1) = R_after_S;
        history.R_eff_after_W(end+1,1) = R_eff_after_W;
        history.R_eff_after_angle(end+1,1) = R_eff_after_angle;
        history.R_eff_after_X(end+1,1) = R_eff_after_X;
        history.R_eff_after_S(end+1,1) = R_eff_after_S;

        R_new = R_eff_after_S;
        history.R_eff(end+1,1) = R_new;
        history.R_sum(end+1,1) = R_after_S;
        history.T_X(end+1,1) = detail_W.T_X;
        history.T_theta(end+1,1) = detail_W.T_theta;
        history.T_phi(end+1,1) = detail_W.T_phi;
        history.T_rec(end+1,1) = detail_W.T_rec;
        history.time_factor(end+1,1) = detail_W.time_factor;

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


elseif strcmp(scheme_mode, 'fixed_antenna_ws')
    history.X_update_mode = 'fixed_antenna_ws';

    X_fixed = state.X;
    theta_fixed = pi * ones(params.N, params.M);
    phi_fixed = zeros(params.N, params.M);
    state.X = X_fixed;
    state.theta = theta_fixed;
    state.phi = phi_fixed;

    for t = 1:params.T_max
        state.t = t;
        state.swap_flag = false;
        state.X = X_fixed;
        state.theta = theta_fixed;
        state.phi = phi_fixed;

        state.W = AO_W(params, scene, model, state);
        R_after_W = Signal_model('sum_rate', params, scene, state, []);
        [R_eff_after_W, ~] = Effective_rate_model(params, scene, state, []);

        R_after_angle = R_after_W;
        R_after_X = R_after_W;
        R_eff_after_angle = R_eff_after_W;
        R_eff_after_X = R_eff_after_W;

        % [state.S, state.swap_flag] = AO_S_fixed(params, scene, model, state);
        [state.S, state.W, state.swap_flag] = AO_S_fixed_reW(params, scene, model, state);

        state.X = X_fixed;
        state.theta = theta_fixed;
        state.phi = phi_fixed;

        R_after_S = Signal_model('sum_rate', params, scene, state, []);
        [R_eff_after_S, detail_S] = Effective_rate_model(params, scene, state, []);

        history.R_after_W(end+1,1) = R_after_W;
        history.R_after_angle(end+1,1) = R_after_angle;
        history.R_after_X(end+1,1) = R_after_X;
        history.R_after_S(end+1,1) = R_after_S;
        history.R_eff_after_W(end+1,1) = R_eff_after_W;
        history.R_eff_after_angle(end+1,1) = R_eff_after_angle;
        history.R_eff_after_X(end+1,1) = R_eff_after_X;
        history.R_eff_after_S(end+1,1) = R_eff_after_S;

        R_new = R_eff_after_S;
        history.R_eff(end+1,1) = R_new;
        history.R_sum(end+1,1) = R_after_S;
        history.T_X(end+1,1) = detail_S.T_X;
        history.T_theta(end+1,1) = detail_S.T_theta;
        history.T_phi(end+1,1) = detail_S.T_phi;
        history.T_rec(end+1,1) = detail_S.T_rec;
        history.time_factor(end+1,1) = detail_S.time_factor;

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
    if ~isfield(history, 'R_eff_after_W')
        history.R_eff_after_W = [];
    end
    if ~isfield(history, 'R_eff_after_angle')
        history.R_eff_after_angle = [];
    end
    if ~isfield(history, 'R_eff_after_X')
        history.R_eff_after_X = [];
    end
    if ~isfield(history, 'R_eff_after_S')
        history.R_eff_after_S = [];
    end
    if ~isfield(history, 'R_before_final_W')
        history.R_before_final_W = [];
    end
    if ~isfield(history, 'R_after_final_W')
        history.R_after_final_W = [];
    end
    R_sum_end = Signal_model('sum_rate', params, scene, state, []);
    [R_eff_end, detail_end] = Effective_rate_model(params, scene, state, []);

    if ~isfield(history,'R_sum') || isempty(history.R_sum)
        history.R_sum = R_sum_end;
    else
        history.R_sum(end,1) = R_sum_end;
    end
    if ~isfield(history,'R_eff') || isempty(history.R_eff)
        history.R_eff = R_eff_end;
    else
        history.R_eff(end,1) = R_eff_end;
    end
    if ~isfield(history,'T_X') || isempty(history.T_X)
        history.T_X = detail_end.T_X;
    else
        history.T_X(end,1) = detail_end.T_X;
    end
    if ~isfield(history,'T_theta') || isempty(history.T_theta)
        history.T_theta = detail_end.T_theta;
    else
        history.T_theta(end,1) = detail_end.T_theta;
    end
    if ~isfield(history,'T_phi') || isempty(history.T_phi)
        history.T_phi = detail_end.T_phi;
    else
        history.T_phi(end,1) = detail_end.T_phi;
    end
    if ~isfield(history,'T_rec') || isempty(history.T_rec)
        history.T_rec = detail_end.T_rec;
    else
        history.T_rec(end,1) = detail_end.T_rec;
    end
    if ~isfield(history,'time_factor') || isempty(history.time_factor)
        history.time_factor = detail_end.time_factor;
    else
        history.time_factor(end,1) = detail_end.time_factor;
    end


elseif strcmp(scheme_mode, 'pso_joint')
    [state_best, history_pso] = PSO_joint(params, scene, model, state);

    state = state_best;
    history = history_pso;

    if ~isfield(history, 'DEBUG_X_cells')
        history.DEBUG_X_cells = {};
    end
    if ~isfield(history, 'X_update_mode')
        history.X_update_mode = 'pso_joint';
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
    if ~isfield(history, 'R_eff_after_W')
        history.R_eff_after_W = [];
    end
    if ~isfield(history, 'R_eff_after_angle')
        history.R_eff_after_angle = [];
    end
    if ~isfield(history, 'R_eff_after_X')
        history.R_eff_after_X = [];
    end
    if ~isfield(history, 'R_eff_after_S')
        history.R_eff_after_S = [];
    end
    if ~isfield(history, 'R_before_final_W')
        history.R_before_final_W = [];
    end
    if ~isfield(history, 'R_after_final_W')
        history.R_after_final_W = [];
    end

    R_sum_end = Signal_model('sum_rate', params, scene, state, []);
    [R_eff_end, detail_end] = Effective_rate_model(params, scene, state, []);

    if ~isfield(history,'R_sum') || isempty(history.R_sum)
        history.R_sum = R_sum_end;
    else
        history.R_sum(end,1) = R_sum_end;
    end
    if ~isfield(history,'R_eff') || isempty(history.R_eff)
        history.R_eff = R_eff_end;
    else
        history.R_eff(end,1) = R_eff_end;
    end
    if ~isfield(history,'T_X') || isempty(history.T_X)
        history.T_X = detail_end.T_X;
    else
        history.T_X(end,1) = detail_end.T_X;
    end
    if ~isfield(history,'T_theta') || isempty(history.T_theta)
        history.T_theta = detail_end.T_theta;
    else
        history.T_theta(end,1) = detail_end.T_theta;
    end
    if ~isfield(history,'T_phi') || isempty(history.T_phi)
        history.T_phi = detail_end.T_phi;
    else
        history.T_phi(end,1) = detail_end.T_phi;
    end
    if ~isfield(history,'T_rec') || isempty(history.T_rec)
        history.T_rec = detail_end.T_rec;
    else
        history.T_rec(end,1) = detail_end.T_rec;
    end
    if ~isfield(history,'time_factor') || isempty(history.time_factor)
        history.time_factor = detail_end.time_factor;
    else
        history.time_factor(end,1) = detail_end.time_factor;
    end

elseif strcmp(scheme_mode, 'hg_multiuser')
    [state_best, history_hg] = HG_multiuser(params, scene, model, state);

    state = state_best;
    history = history_hg;

    if ~isfield(history, 'DEBUG_X_cells')
        history.DEBUG_X_cells = {};
    end
    if ~isfield(history, 'X_update_mode')
        history.X_update_mode = 'hungarian_greedy';
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
    if ~isfield(history, 'R_eff_after_W')
        history.R_eff_after_W = [];
    end
    if ~isfield(history, 'R_eff_after_angle')
        history.R_eff_after_angle = [];
    end
    if ~isfield(history, 'R_eff_after_X')
        history.R_eff_after_X = [];
    end
    if ~isfield(history, 'R_eff_after_S')
        history.R_eff_after_S = [];
    end
    if ~isfield(history, 'R_before_final_W')
        history.R_before_final_W = [];
    end
    if ~isfield(history, 'R_after_final_W')
        history.R_after_final_W = [];
    end
    R_sum_end = Signal_model('sum_rate', params, scene, state, []);
    [R_eff_end, detail_end] = Effective_rate_model(params, scene, state, []);

    if ~isfield(history,'R_sum') || isempty(history.R_sum)
        history.R_sum = R_sum_end;
    else
        history.R_sum(end,1) = R_sum_end;
    end
    if ~isfield(history,'R_eff') || isempty(history.R_eff)
        history.R_eff = R_eff_end;
    else
        history.R_eff(end,1) = R_eff_end;
    end
    if ~isfield(history,'T_X') || isempty(history.T_X)
        history.T_X = detail_end.T_X;
    else
        history.T_X(end,1) = detail_end.T_X;
    end
    if ~isfield(history,'T_theta') || isempty(history.T_theta)
        history.T_theta = detail_end.T_theta;
    else
        history.T_theta(end,1) = detail_end.T_theta;
    end
    if ~isfield(history,'T_phi') || isempty(history.T_phi)
        history.T_phi = detail_end.T_phi;
    else
        history.T_phi(end,1) = detail_end.T_phi;
    end
    if ~isfield(history,'T_rec') || isempty(history.T_rec)
        history.T_rec = detail_end.T_rec;
    else
        history.T_rec(end,1) = detail_end.T_rec;
    end
    if ~isfield(history,'time_factor') || isempty(history.time_factor)
        history.time_factor = detail_end.time_factor;
    else
        history.time_factor(end,1) = detail_end.time_factor;
    end

else
    error('main: unsupported scheme_mode');
end

%% 第7部分：整理输出
result = struct();
result.state = state;
result.history = history;
result.params = params;
result.scene = scene;
result.model = model;

%% 第8部分：结果显示方式切换
% 方式1：新的论文式多方案对比图（默认启用）
% compare_result = Plot_Compare(params, scene);
% result.compare_result = compare_result;

% 方式2：原来的单次仿真结果图
% Print_and_Plot(params, scene, model, result);

end
