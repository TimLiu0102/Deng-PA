function compare_result = Plot_Compare(base_params)
% Plot_Compare：五种方案的论文式对比绘图
% 支持 debug/full 模式切换
% 支持不同图分别开关
% 不改已有算法文件，只在本文件内部重复调用现有初始化、AO模块和SA模块

%% ======================== 画图模式切换 ========================
plot_mode = 'debug';   % 'debug' 或 'full'

%% ======================== 图形开关 ========================
do_snr         = true;   % 图1：频谱效率 vs SNR
do_K           = true;   % 图2：频谱效率 vs 用户数 K
do_N           = true;   % 图3：频谱效率 vs 波导数 N
do_M           = true;   % 图4：频谱效率 vs 每条波导 PA 数 M
do_convergence = true;   % 图5：收敛曲线
do_cdf         = true;   % 图6：Per-user rate CDF
do_geometry    = true;   % 图7：几何图

fprintf('\n================ 多方案对比绘图 ================\n');

% 五种方案列表
schemes = build_scheme_list();

% debug 用于调试，full 用于正式出图
if strcmp(plot_mode, 'debug')
    MC = 1;
    snr_dB_vec = [0 10 20];
    K_vec = [16 32];
    N_vec = [4 8];
    M_vec = [2 4];
else
    MC = 3;
    snr_dB_vec = [-10 -5 0 5 10 15 20 25 30];
    K_vec = [8 16 24 32 48 64];
    N_vec = [2 4 6 8 10 12];
    M_vec = [2 4 6 8];
end

% 保证 K 不小于 K_serv
K_vec = max(K_vec, base_params.K_serv);

compare_result = struct();
compare_result.plot_mode = plot_mode;
compare_result.MC = MC;
compare_result.schemes = schemes;

% 图1：SNR 扫描
if do_snr
    [mean_R, std_R, R_all] = run_sweep(base_params, schemes, snr_dB_vec, 'snr', MC);
    figure('Name', 'Fig1_SNR');
    draw_mean_error_curve(snr_dB_vec, mean_R, std_R, schemes, 'SNR (dB)', 'Average spectral efficiency versus SNR');

    compare_result.snr = struct();
    compare_result.snr.x = snr_dB_vec;
    compare_result.snr.mean_R = mean_R;
    compare_result.snr.std_R = std_R;
    compare_result.snr.R_all = R_all;
end

% 图2：K 扫描
if do_K
    [mean_R, std_R, R_all] = run_sweep(base_params, schemes, K_vec, 'K', MC);
    figure('Name', 'Fig2_K');
    draw_mean_error_curve(K_vec, mean_R, std_R, schemes, 'Number of users K', 'Average spectral efficiency versus number of users');

    compare_result.K = struct();
    compare_result.K.x = K_vec;
    compare_result.K.mean_R = mean_R;
    compare_result.K.std_R = std_R;
    compare_result.K.R_all = R_all;
end

% 图3：N 扫描
if do_N
    [mean_R, std_R, R_all] = run_sweep(base_params, schemes, N_vec, 'N', MC);
    figure('Name', 'Fig3_N');
    draw_mean_error_curve(N_vec, mean_R, std_R, schemes, 'Number of waveguides N', 'Average spectral efficiency versus number of waveguides');

    compare_result.N = struct();
    compare_result.N.x = N_vec;
    compare_result.N.mean_R = mean_R;
    compare_result.N.std_R = std_R;
    compare_result.N.R_all = R_all;
end

% 图4：M 扫描
if do_M
    [mean_R, std_R, R_all] = run_sweep(base_params, schemes, M_vec, 'M', MC);
    figure('Name', 'Fig4_M');
    draw_mean_error_curve(M_vec, mean_R, std_R, schemes, 'Number of PAs per waveguide M', 'Average spectral efficiency versus number of PAs per waveguide');

    compare_result.M = struct();
    compare_result.M.x = M_vec;
    compare_result.M.mean_R = mean_R;
    compare_result.M.std_R = std_R;
    compare_result.M.R_all = R_all;
end

% 图5：收敛曲线
if do_convergence
    conv_results = run_convergence_cases(base_params, schemes);
    draw_convergence(conv_results, schemes);

    compare_result.convergence = conv_results;
end

% 图6：CDF
if do_cdf
    rate_cells = collect_rate_cdf_data(base_params, schemes, MC);
    figure('Name', 'Fig6_CDF');
    draw_rate_cdf(rate_cells, schemes);

    compare_result.cdf = struct();
    compare_result.cdf.rate_cells = rate_cells;
end

% 图7：几何图（只画“所提初始化 + 所提优化算法”）
if do_geometry
    idx_geo = 2;
    seed_geo = base_params.seed + 10000*7 + 1000*1 + 100*idx_geo + 1;
    geo_result = run_one_case(base_params, schemes(idx_geo).init_mode, schemes(idx_geo).alg_mode, seed_geo);

    figure('Name', 'Fig7_Geometry');
    draw_geometry_case(geo_result);

    compare_result.geometry = geo_result;
end

end

%% ======================== 本地子函数 ========================
function schemes = build_scheme_list()
% 方案列表：五种对比方案
schemes = struct('name', {}, 'init_mode', {}, 'alg_mode', {});

schemes(1).name = 'Random init + Proposed AO';
schemes(1).init_mode = 'random';
schemes(1).alg_mode = 'proposed';

schemes(2).name = 'Proposed init + Proposed AO';
schemes(2).init_mode = 'proposed';
schemes(2).alg_mode = 'proposed';

schemes(3).name = 'Random init + W only';
schemes(3).init_mode = 'random';
schemes(3).alg_mode = 'w_only';

schemes(4).name = 'Proposed init + W only';
schemes(4).init_mode = 'proposed';
schemes(4).alg_mode = 'w_only';

schemes(5).name = 'Random init + SA joint';
schemes(5).init_mode = 'random';
schemes(5).alg_mode = 'sa_joint';
end

function [mean_R, std_R, R_all] = run_sweep(base_params, schemes, x_vec, sweep_type, MC)
% 对一个横坐标向量进行批量扫描
ns = numel(schemes);
nx = numel(x_vec);
R_all = zeros(nx, ns, MC);

sweep_id = get_sweep_id(sweep_type);

for idx_x = 1:nx
    x_val = x_vec(idx_x);
    params_case = make_params_for_sweep(base_params, sweep_type, x_val);

    for idx_scheme = 1:ns
        for idx_mc = 1:MC
            seed_now = base_params.seed + 10000*sweep_id + 1000*idx_x + 100*idx_scheme + idx_mc;
            out_case = run_one_case(params_case, schemes(idx_scheme).init_mode, schemes(idx_scheme).alg_mode, seed_now);
            R_all(idx_x, idx_scheme, idx_mc) = out_case.final_R;
        end
    end
end

mean_R = squeeze(mean(R_all, 3));
std_R = squeeze(std(R_all, 0, 3));

if isvector(mean_R)
    mean_R = mean_R(:);
end
if isvector(std_R)
    std_R = std_R(:);
end
end

function params_case = make_params_for_sweep(base_params, sweep_type, x_value)
% 根据扫描类型构造 params_case
params_case = base_params;

if strcmp(sweep_type, 'snr')
    snr_ref_dB = 20;
    sigma2_ref = base_params.sigma2;
    params_case.sigma2 = sigma2_ref * 10.^((snr_ref_dB - x_value)/10);
elseif strcmp(sweep_type, 'K')
    params_case.K = x_value;
elseif strcmp(sweep_type, 'N')
    params_case.N = x_value;
    params_case.NRF = x_value;
    params_case.K_max = x_value;
    params_case.K_serv = min(params_case.NRF, params_case.K_max);
elseif strcmp(sweep_type, 'M')
    params_case.M = x_value;
else
    error('make_params_for_sweep: unsupported sweep_type');
end
end

function out_case = run_one_case(params_case, init_mode, alg_mode, seed_now)
% 单个实验：初始化 + 算法运行 + 指标整理
rng(seed_now);
scene = Channel_model('build_scene', params_case, [], [], []);
model = Problem_formulation(params_case, scene);

% 初始化方式
if strcmp(init_mode, 'random')
    state = Initialization_ra(params_case, scene, model);
else
    state = Initialization(params_case, scene, model);
end

% 历史量初始化
history = init_history(params_case, scene, state);

% 算法模式
if strcmp(alg_mode, 'proposed')
    [state, history] = run_proposed_ao(params_case, scene, model, state, history);
elseif strcmp(alg_mode, 'w_only')
    [state, history] = run_w_only(params_case, scene, model, state, history);
elseif strcmp(alg_mode, 'sa_joint')
    [state, history] = run_sa_joint_case(params_case, scene, model, state);
else
    error('run_one_case: unsupported alg_mode');
end

% 统一输出结构
final_R = Signal_model('sum_rate', params_case, scene, state, []);
rates_final = Signal_model('individual_rates', params_case, scene, state, []);

out_case = struct();
out_case.params = params_case;
out_case.scene = scene;
out_case.model = model;
out_case.state = state;
out_case.history = history;
out_case.final_R = final_R;
out_case.rates_final = rates_final;
end

function history = init_history(params_case, scene, state)
% 复用 main.m 的 history 初始化风格
if ~isfield(state, 'swap_flag')
    state.swap_flag = false;
end

rates0 = Signal_model('individual_rates', params_case, scene, state, []);
R_old = sum(rates0);

history = struct();
history.X0 = state.X;
history.theta0 = state.theta;
history.phi0 = state.phi;
history.S0 = state.S;
history.rates0 = rates0;
history.R_sum = R_old;

history.R_after_W = [];
history.R_after_angle = [];
history.R_after_X = [];
history.R_after_S = [];

history.S_cells = {};
history.X_cells = {};
history.theta_cells = {};
history.phi_cells = {};
history.DEBUG_X_cells = {};
history.swap_flag = false;
end

function [state, history] = run_proposed_ao(params_case, scene, model, state, history)
% 对应 main.m 中 ao_final_w 逻辑
if ~isfield(state, 'swap_flag')
    state.swap_flag = false;
end

R_old = history.R_sum(end,1);
history.X_update_mode = 'gradient';

for t = 1:params_case.T_max
    state.t = t;
    R_after_W = R_old;

    [state.theta, state.phi] = AO_angle(params_case, scene, model, state);
    R_after_angle = Signal_model('sum_rate', params_case, scene, state, []);

    [state.X, DEBUG_X_t] = AO_X(params_case, scene, model, state);
    R_after_X = Signal_model('sum_rate', params_case, scene, state, []);

    [state.S, state.swap_flag] = AO_S(params_case, scene, model, state);
    R_after_S = Signal_model('sum_rate', params_case, scene, state, []);

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
    history.DEBUG_X_cells{t,1} = DEBUG_X_t;

    if abs(R_new - R_old) < params_case.eps_outer
        break;
    end

    R_old = R_new;
end

R_before_final_W = Signal_model('sum_rate', params_case, scene, state, []);
state.W = AO_W(params_case, scene, model, state);
R_after_final_W = Signal_model('sum_rate', params_case, scene, state, []);

history.R_before_final_W = R_before_final_W;
history.R_after_final_W = R_after_final_W;
history.R_sum(end,1) = R_after_final_W;
end

function [state, history] = run_w_only(params_case, scene, model, state, history)
% 对应 main.m 中 w_only 逻辑
if ~isfield(state, 'swap_flag')
    state.swap_flag = false;
end

R_old = history.R_sum(end,1);
history.X_update_mode = 'none';
history.R_before_final_W = [];
history.R_after_final_W = [];

for t = 1:params_case.T_max
    state.t = t;
    state.swap_flag = false;

    state.W = AO_W(params_case, scene, model, state);
    R_after_W = Signal_model('sum_rate', params_case, scene, state, []);

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

    if abs(R_new - R_old) < params_case.eps_outer
        break;
    end

    R_old = R_new;
end
end

function [state, history] = run_sa_joint_case(params_case, scene, model, state)
% 直接调用 SA_joint，并整理成统一历史字段
[state_best, history_sa] = SA_joint(params_case, scene, model, state);
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
history.R_sum(end,1) = Signal_model('sum_rate', params_case, scene, state, []);
end

function conv_results = run_convergence_cases(base_params, schemes)
% 固定参数，每种方案跑一次并收敛曲线对比
ns = numel(schemes);
conv_results = cell(ns,1);

for idx_scheme = 1:ns
    seed_now = base_params.seed + 10000*5 + 1000*1 + 100*idx_scheme + 1;
    out_case = run_one_case(base_params, schemes(idx_scheme).init_mode, schemes(idx_scheme).alg_mode, seed_now);
    conv_results{idx_scheme} = out_case.history.R_sum(:);
end
end

function rate_cells = collect_rate_cdf_data(base_params, schemes, MC)
% 固定参数，多次 Monte Carlo 收集所有服务用户速率
ns = numel(schemes);
rate_cells = cell(ns,1);

for idx_scheme = 1:ns
    rates_all = [];
    for idx_mc = 1:MC
        seed_now = base_params.seed + 10000*6 + 1000*1 + 100*idx_scheme + idx_mc;
        out_case = run_one_case(base_params, schemes(idx_scheme).init_mode, schemes(idx_scheme).alg_mode, seed_now);
        rates_all = [rates_all; out_case.rates_final(:)]; %#ok<AGROW>
    end
    rate_cells{idx_scheme} = rates_all;
end
end

function draw_mean_error_curve(x_vec, mean_R, std_R, schemes, x_label_text, title_text)
% 均值 + 标准差误差棒图
for s = 1:numel(schemes)
    errorbar(x_vec, mean_R(:,s), std_R(:,s), '-o', 'LineWidth', 1.2);
    hold on;
end
xlabel(x_label_text);
ylabel('Average spectral efficiency (bit/s/Hz)');
title(title_text);
legend({schemes.name}, 'Location', 'bestoutside');
grid on;
end

function draw_convergence(conv_results, schemes)
% 五种方案收敛曲线，使用 broken x-axis 思路压缩长迭代区间
% 0~30 正常显示，30 之后压缩显示，避免 SA_joint 把 AO 曲线挤在左侧
figure('Name', 'Fig5_Convergence_BrokenAxis');

break_iter = 30;
compress_ratio = 0.08;

for s = 1:numel(schemes)
    r = conv_results{s};
    it_real = 0:(numel(r)-1);

    % 将真实迭代次数映射到显示用横坐标
    it_plot = it_real;
    idx = it_real > break_iter;
    it_plot(idx) = break_iter + (it_real(idx) - break_iter) * compress_ratio;

    plot(it_plot, r, '-o', 'LineWidth', 1.2);
    hold on;
end

% 设置横坐标刻度：显示真实迭代次数，但位置是压缩后的
real_ticks = [0 10 20 30 100 200 300 400 500];
plot_ticks = real_ticks;
idx_tick = real_ticks > break_iter;
plot_ticks(idx_tick) = break_iter + (real_ticks(idx_tick) - break_iter) * compress_ratio;

set(gca, 'XTick', plot_ticks);
set(gca, 'XTickLabel', string(real_ticks));

% 在断轴位置画虚线
yl = ylim;
plot([break_iter break_iter], yl, 'k--', 'LineWidth', 1.0);
ylim(yl);

xlabel('Iteration index');
ylabel('Spectral efficiency (bit/s/Hz)');
title('Convergence behavior of different schemes with compressed x-axis');
legend({schemes.name}, 'Location', 'bestoutside');
grid on;

% 在图中标注横坐标被压缩
text(break_iter + 2, yl(1) + 0.08*(yl(2)-yl(1)), ...
    'x-axis compressed after 30 iterations', ...
    'FontSize', 9);
end

function draw_rate_cdf(rate_cells, schemes)
% 五种方案的单用户速率 CDF
for s = 1:numel(schemes)
    r = sort(rate_cells{s}(:));
    F = (1:numel(r)) / numel(r);
    plot(r, F, 'LineWidth', 1.2);
    hold on;
end
xlabel('Per-user rate (bit/s/Hz)');
ylabel('CDF');
title('CDF of per-user rate');
legend({schemes.name}, 'Location', 'bestoutside');
grid on;
end

function draw_geometry_case(geo_result)
% 代表性几何图：所提初始化 + 所提优化算法
params_case = geo_result.params;
scene = geo_result.scene;
state = geo_result.state;
history = geo_result.history;

user_pos = scene.user_pos;
S = state.S;
M = params_case.M;

% 全部用户
scatter(user_pos(1,:), user_pos(2,:), 25, 'filled');
hold on;

% 服务用户
scatter(user_pos(1,S), user_pos(2,S), 70);

% 波导、初始PA、最终PA
for n = 1:params_case.N
    line([scene.xW(n), scene.xW(n)], [0, params_case.Dy]);
    plot(scene.xW(n)*ones(1,M), history.X0(n,:), 'o');
    plot(scene.xW(n)*ones(1,M), state.X(n,:), 'x');
end

% 最终PA朝向箭头：画三维方向向量在 x-y 平面的投影
% theta 是与 z 轴的夹角，phi 是 x-y 平面内的方位角
% 因此方向向量为 [sin(theta)cos(phi), sin(theta)sin(phi), cos(theta)]
% 二维几何图只显示 x-y 平面，所以取前两个分量
x_pa = repmat(scene.xW(:), 1, M);
x_pa = x_pa(:);
y_pa = state.X(:);

theta_vec = state.theta(:);
phi_vec = state.phi(:);

u = sin(theta_vec) .* cos(phi_vec);
v = sin(theta_vec) .* sin(phi_vec);

% 只保留方向，统一箭头长度，避免箭头长短影响观察
uv_norm = sqrt(u.^2 + v.^2);
u = u ./ (uv_norm + eps);
v = v ./ (uv_norm + eps);

quiver(x_pa, y_pa, u, v, 0.6, 'LineWidth', 0.8);

xlabel('x (m)');
ylabel('y (m)');
title('Geometry and final PA/user configuration');
legend({'All users', 'Served users', 'Waveguide', 'Initial PA', 'Final PA', 'PA orientation'}, 'Location', 'eastoutside');
grid on;
end

function sweep_id = get_sweep_id(sweep_type)
% 给不同实验分配 sweep_id，便于构造不同随机种子
if strcmp(sweep_type, 'snr')
    sweep_id = 1;
elseif strcmp(sweep_type, 'K')
    sweep_id = 2;
elseif strcmp(sweep_type, 'N')
    sweep_id = 3;
elseif strcmp(sweep_type, 'M')
    sweep_id = 4;
else
    sweep_id = 9;
end
end
