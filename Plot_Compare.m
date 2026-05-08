function compare_result = Plot_Compare(base_params, base_scene)
% Plot_Compare：多方案对比绘图

if nargin < 2 || isempty(base_scene)
    rng(base_params.seed);
    base_scene = Channel_model('build_scene', base_params, [], [], []);
end

plot_mode = 'debug';   % 'debug' 或 'full'
% debug 模式只减少 MC，不减少横轴取值；如果调试 PSO 较慢，可手动关闭 do_N/do_Dy。

do_snr         = false;
do_K           = false;
do_N           = false;
do_M           = false;
do_Dy          = false;
do_convergence = true;
do_cdf         = false;
do_final_bar_ab = false;
do_H2_ab = false;
do_default_geometry = false;
do_default_check = false;

fprintf('\n================ 多方案对比绘图 ================\n');

schemes = build_scheme_list();

if isfield(base_params, 'SNR_dB')
    snr_ref_dB = base_params.SNR_dB;
else
    snr_ref_dB = 20;   % 当前绘图口径：base_params.sigma2 对应参考 SNR=20 dB
end

snr_dB_base_vec = [-10 -5 0 5 10 15 20 25 30];
K_base_vec      = [8 16 24 32 48 64];
N_base_vec      = [2 4 6 8 10 12];
M_base_vec      = [2 4 6 8];
Dy_base_vec     = [4 8 12 16 20];
add_default_point = true;

if isfield(base_params, 'waveguide_Dy')
    Dy_default = base_params.waveguide_Dy;
else
    Dy_default = base_params.Dy;
end
if add_default_point
    snr_dB_vec = unique(sort([snr_dB_base_vec snr_ref_dB]));
    K_vec      = unique(sort([K_base_vec base_params.K]));
    N_vec      = unique(sort([N_base_vec base_params.N]));
    M_vec      = unique(sort([M_base_vec base_params.M]));
    Dy_vec     = unique(sort([Dy_base_vec Dy_default]));
else
    snr_dB_vec = snr_dB_base_vec;
    K_vec      = K_base_vec;
    N_vec      = N_base_vec;
    M_vec      = M_base_vec;
    Dy_vec     = Dy_base_vec;
end

if strcmp(plot_mode, 'debug')
    MC = 3;
else
    MC = 30;
end

K_vec = unique(max(K_vec, base_params.K_serv));

% 注意：
% default_check 使用 base_scene，用于和 main 单次默认结果做 sanity check。
% sweep 图使用 MC user_pos_pools。它们不要求等于 main 单次结果，
% 但同一个 MC 下，不同 sweep 图中的默认参数点应使用相同用户分布和相同 seed，
% 因此这些默认点应彼此一致。
max_K_for_pool = max([K_vec(:); base_params.K]);
user_pos_pools = cell(MC,1);
for idx_mc = 1:MC
    params_pool = base_params;
    params_pool.K = max_K_for_pool;
    scene_seed = base_params.seed + 10000 + idx_mc;
    rng(scene_seed);
    scene_pool = Channel_model('build_scene', params_pool, [], [], []);
    user_pos_pools{idx_mc} = scene_pool.user_pos;
end

compare_result = struct();
compare_result.plot_mode = plot_mode;
compare_result.MC = MC;
compare_result.schemes = schemes;

if do_default_check
    fprintf('\n================ Default sanity check ================\n');
    for idx_scheme = 1:numel(schemes)
        init_seed_case = base_params.seed;
        alg_seed_case  = base_params.seed + 100*idx_scheme;

        out_default = run_one_case(base_params, ...
            schemes(idx_scheme).init_mode, ...
            schemes(idx_scheme).alg_mode, ...
            init_seed_case, alg_seed_case, base_scene);

        fprintf('%-18s: Rsum = %.6f, Reff = %.6f, Trec = %.6f, time_factor = %.6f\n', ...
            schemes(idx_scheme).name, ...
            out_default.final_R_sum, ...
            out_default.final_R_eff, ...
            out_default.final_detail.T_rec, ...
            out_default.final_detail.time_factor);

        compare_result.default_check(idx_scheme).name = schemes(idx_scheme).name;
        compare_result.default_check(idx_scheme).R_sum = out_default.final_R_sum;
        compare_result.default_check(idx_scheme).R_eff = out_default.final_R_eff;
        compare_result.default_check(idx_scheme).T_rec = out_default.final_detail.T_rec;
        compare_result.default_check(idx_scheme).time_factor = out_default.final_detail.time_factor;
    end
end

if do_snr
    [mean_R, std_R, R_all, mean_R_sum, std_R_sum, R_all_sum] = run_sweep(base_params, schemes, snr_dB_vec, 'snr', MC, user_pos_pools);
    figure('Name', 'Fig1_SNR', 'Position', [100 100 760 520]);
    draw_mean_error_curve(snr_dB_vec, mean_R, std_R, schemes, 'SNR (dB)', 'Effective Spectral Efficiency vs. SNR');
    compare_result.snr = pack_sweep_result(snr_dB_vec, mean_R, std_R, R_all, mean_R_sum, std_R_sum, R_all_sum);
end

if do_K
    [mean_R, std_R, R_all, mean_R_sum, std_R_sum, R_all_sum] = run_sweep(base_params, schemes, K_vec, 'K', MC, user_pos_pools);
    figure('Name', 'Fig2_K', 'Position', [100 100 760 520]);
    draw_mean_error_curve(K_vec, mean_R, std_R, schemes, 'Number of users K', 'Effective Spectral Efficiency vs. Number of Users');
    compare_result.K = pack_sweep_result(K_vec, mean_R, std_R, R_all, mean_R_sum, std_R_sum, R_all_sum);
end

if do_N
    params_N = base_params;

    [mean_R, std_R, R_all, mean_R_sum, std_R_sum, R_all_sum] = run_sweep(params_N, schemes, N_vec, 'N', MC, user_pos_pools);
    figure('Name', 'Fig3_N', 'Position', [100 100 760 520]);
    draw_mean_error_curve(N_vec, mean_R, std_R, schemes, 'Number of waveguides N', 'Effective Spectral Efficiency vs. Number of Waveguides');
    compare_result.N = pack_sweep_result(N_vec, mean_R, std_R, R_all, mean_R_sum, std_R_sum, R_all_sum);
end

if do_M
    [mean_R, std_R, R_all, mean_R_sum, std_R_sum, R_all_sum] = run_sweep(base_params, schemes, M_vec, 'M', MC, user_pos_pools);
    figure('Name', 'Fig4_M', 'Position', [100 100 760 520]);
    draw_mean_error_curve(M_vec, mean_R, std_R, schemes, 'Number of PAs per waveguide M', 'Effective Spectral Efficiency vs. Number of PAs');
    compare_result.M = pack_sweep_result(M_vec, mean_R, std_R, R_all, mean_R_sum, std_R_sum, R_all_sum);
end

if do_Dy
    [mean_R, std_R, R_all, mean_R_sum, std_R_sum, R_all_sum] = run_sweep(base_params, schemes, Dy_vec, 'Dy', MC, user_pos_pools);
    figure('Name', 'Fig5_Dy', 'Position', [100 100 760 520]);
    draw_mean_error_curve(Dy_vec, mean_R, std_R, schemes, 'Waveguide length / PA movable range D_y (m)', 'Effective Spectral Efficiency vs. Waveguide Length');
    compare_result.Dy = pack_sweep_result(Dy_vec, mean_R, std_R, R_all, mean_R_sum, std_R_sum, R_all_sum);
end

if do_convergence
    conv_results = run_convergence_cases(base_params, schemes);
    draw_convergence(conv_results, schemes);
    compare_result.convergence = conv_results;
end

if do_cdf
    rate_cells = collect_rate_cdf_data(base_params, schemes, MC, user_pos_pools);
    figure('Name', 'Fig6_CDF', 'Position', [100 100 760 520]);
    draw_rate_cdf(rate_cells, schemes);
    compare_result.cdf.rate_cells = rate_cells;
end


if do_default_geometry
    idx_geo = 1;
    scene_case = base_scene;
    geo_result = run_one_case(base_params, schemes(idx_geo).init_mode, schemes(idx_geo).alg_mode, ...
        base_params.seed, base_params.seed + 100*idx_geo, scene_case);
    figure('Name', 'Fig_Default_Geometry', 'Position', [100 100 760 520]);
    draw_geometry_case(geo_result);
    compare_result.geometry = geo_result;
end

if do_final_bar_ab
    final_bar_ab = run_final_bar_ab_cases(base_params, schemes, MC, user_pos_pools);
    draw_final_bar_ab(final_bar_ab, schemes);
    compare_result.final_bar_ab = final_bar_ab;
end

if do_H2_ab
    H2_ab = draw_H2_ab_cases(base_params);
    compare_result.H2_ab = H2_ab;
end

end

function s = pack_sweep_result(x,mean_R,std_R,R_all,mean_R_sum,std_R_sum,R_all_sum)
s = struct('x',x,'mean_R',mean_R,'std_R',std_R,'R_all',R_all,...
    'mean_R_sum',mean_R_sum,'std_R_sum',std_R_sum,'R_all_sum',R_all_sum);
end

function schemes = build_scheme_list()
schemes = struct('name', {}, 'init_mode', {}, 'alg_mode', {});
schemes(1).name = 'Proposed AO'; schemes(1).init_mode = 'paper'; schemes(1).alg_mode = 'AO';
schemes(2).name = 'Fixed W+S'; schemes(2).init_mode = 'uniform_fixed'; schemes(2).alg_mode = 'fixed_antenna_ws';
schemes(3).name = 'Fixed W+S reW'; schemes(3).init_mode = 'uniform_fixed'; schemes(3).alg_mode = 'fixed_antenna_ws_reW';
schemes(4).name = 'HG-Rsum'; schemes(4).init_mode = 'uniform_neutral'; schemes(4).alg_mode = 'hg_multiuser';
schemes(5).name = 'SA joint'; schemes(5).init_mode = 'uniform_neutral'; schemes(5).alg_mode = 'sa_joint';
% schemes(6).name = 'PSO joint'; schemes(6).init_mode = 'uniform_neutral'; schemes(6).alg_mode = 'pso_joint';
end

function [mean_R, std_R, R_all_eff, mean_R_sum, std_R_sum, R_all_sum] = run_sweep(base_params, schemes, x_vec, sweep_type, MC, user_pos_pools)
ns = numel(schemes); nx = numel(x_vec);
R_all_eff = zeros(nx, ns, MC); R_all_sum = zeros(nx, ns, MC);
for idx_mc = 1:MC
    user_pos_pool = user_pos_pools{idx_mc};
    % 同一个 MC 下，所有 sweep 图、所有横坐标点、所有方案使用同一个 user_pos_pool。
    for idx_x = 1:nx
        params_case = make_params_for_sweep(base_params, sweep_type, x_vec(idx_x));
        scene_case = build_scene_with_fixed_users(params_case, user_pos_pool);
        for idx_scheme = 1:ns
            init_seed_case = base_params.seed + 20000 + idx_mc;
            alg_seed_case  = base_params.seed + 30000 + 100*idx_scheme + idx_mc;
            out_case = run_one_case(params_case, schemes(idx_scheme).init_mode, schemes(idx_scheme).alg_mode,...
                init_seed_case, alg_seed_case, scene_case);
            R_all_eff(idx_x, idx_scheme, idx_mc) = out_case.final_R_eff;
            R_all_sum(idx_x, idx_scheme, idx_mc) = out_case.final_R_sum;
        end
    end
end
mean_R = squeeze(mean(R_all_eff,3)); std_R = squeeze(std(R_all_eff,0,3));
mean_R_sum = squeeze(mean(R_all_sum,3)); std_R_sum = squeeze(std(R_all_sum,0,3));
end

function params_case = make_params_for_sweep(base_params, sweep_type, x_value)
params_case = base_params;
if strcmp(sweep_type, 'snr')
    if isfield(base_params, 'SNR_dB')
        snr_ref_dB = base_params.SNR_dB;
    else
        snr_ref_dB = 20;   % 与上面保持一致
    end
    params_case.sigma2 = base_params.sigma2 * 10.^((snr_ref_dB - x_value)/10);
elseif strcmp(sweep_type, 'K')
    params_case.K = x_value;
elseif strcmp(sweep_type, 'N')
    params_case.N = x_value;
elseif strcmp(sweep_type, 'M')
    params_case.M = x_value;
elseif strcmp(sweep_type, 'Dy')
    if isfield(params_case, 'waveguide_Dy')
        params_case.waveguide_Dy = x_value;
        params_case.Dy = params_case.waveguide_Dy;
    else
        params_case.Dy = x_value;
    end
else
    error('unsupported sweep_type');
end
end

function out_case = run_one_case(params_case, init_mode, alg_mode, init_seed, alg_seed, scene_in)
if nargin >= 6 && ~isempty(scene_in), scene = scene_in; else, rng(init_seed); scene = Channel_model('build_scene', params_case, [], [], []); end
model = Problem_formulation(params_case, scene);
rng(init_seed);
if strcmp(init_mode, 'paper')
    state = Initialization(params_case, scene, model);
elseif strcmp(init_mode, 'uniform_neutral')
    state = Initialization_uniform_neutral(params_case, scene, model);
elseif strcmp(init_mode, 'uniform_fixed')
    state = Initialization_uniform_fixed(params_case, scene, model);
elseif strcmp(init_mode, 'uniform')
    state = Initialization_uniform(params_case, scene, model);
elseif strcmp(init_mode, 'random')
    state = Initialization_ra(params_case, scene, model);
else
    error('unsupported init_mode');
end
rng(alg_seed);
if strcmp(alg_mode, 'AO')
    [state, history] = run_AO_case(params_case, scene, model, state);
elseif strcmp(alg_mode, 'fixed_antenna_ws')
    [state, history] = run_fixed_antenna_ws_case(params_case, scene, model, state);
elseif strcmp(alg_mode, 'fixed_antenna_ws_reW')
    [state, history] = run_fixed_antenna_ws_reW_case(params_case, scene, model, state);
elseif strcmp(alg_mode, 'hg_multiuser')
    [state, history] = HG_multiuser(params_case, scene, model, state);
elseif strcmp(alg_mode, 'sa_joint')
    [state, history] = run_sa_joint_case(params_case, scene, model, state);
elseif strcmp(alg_mode, 'pso_joint')
    [state, history] = run_pso_joint_case(params_case, scene, model, state);
else
    error('unsupported alg_mode');
end
final_R_sum = Signal_model('sum_rate', params_case, scene, state, []);
[final_R_eff, final_detail] = Effective_rate_model(params_case, scene, state, []);
rates_final = Signal_model('individual_rates', params_case, scene, state, []);
out_case = struct('params',params_case,'scene',scene,'model',model,'state',state,'history',history,...
    'final_R',final_R_eff,'final_R_eff',final_R_eff,'final_R_sum',final_R_sum,'final_detail',final_detail,'rates_final',rates_final);
end

function [state, history] = run_AO_case(params, scene, model, state)
if ~isfield(state,'swap_flag'), state.swap_flag = false; end
[Reff0, d0] = Effective_rate_model(params, scene, state, []);
history = init_history_full(params, scene, state, Reff0, d0);
for t=1:params.T_max
    state.t = t;
    state.W = AO_W(params, scene, model, state); [R1e,~] = Effective_rate_model(params, scene, state, []); R1 = Signal_model('sum_rate', params, scene, state, []);
    [state.theta,state.phi] = AO_angle(params, scene, model, state); [R2e,~] = Effective_rate_model(params, scene, state, []); R2 = Signal_model('sum_rate', params, scene, state, []);
    [state.X,dbg] = AO_X(params, scene, model, state); [R3e,~] = Effective_rate_model(params, scene, state, []); R3 = Signal_model('sum_rate', params, scene, state, []);
    [state.S,state.swap_flag] = AO_S(params, scene, model, state); [R4e,d4] = Effective_rate_model(params, scene, state, []); R4 = Signal_model('sum_rate', params, scene, state, []);
    history.R_after_W(end+1,1)=R1; history.R_after_angle(end+1,1)=R2; history.R_after_X(end+1,1)=R3; history.R_after_S(end+1,1)=R4;
    history.R_eff_after_W(end+1,1)=R1e; history.R_eff_after_angle(end+1,1)=R2e; history.R_eff_after_X(end+1,1)=R3e; history.R_eff_after_S(end+1,1)=R4e;
    history.R_sum(end+1,1)=R4; history.R_eff(end+1,1)=R4e; history.T_X(end+1,1)=d4.T_X; history.T_theta(end+1,1)=d4.T_theta; history.T_phi(end+1,1)=d4.T_phi; history.T_rec(end+1,1)=d4.T_rec; history.time_factor(end+1,1)=d4.time_factor;
    history.S_cells{t,1}=state.S; history.X_cells{t,1}=state.X; history.theta_cells{t,1}=state.theta; history.phi_cells{t,1}=state.phi; history.swap_flag(end+1,1)=state.swap_flag; history.DEBUG_X_cells{t,1}=dbg;
    if abs(history.R_eff(end)-history.R_eff(end-1)) < params.eps_outer, break; end
end
end

function [state, history] = run_fixed_antenna_ws_case(params, scene, model, state)
if ~isfield(state,'swap_flag'), state.swap_flag = false; end
X_fixed = state.X; theta_fixed = pi*ones(params.N,params.M); phi_fixed = zeros(params.N,params.M);
state.theta=theta_fixed; state.phi=phi_fixed;
[Reff0, d0] = Effective_rate_model(params, scene, state, []);
history = init_history_full(params, scene, state, Reff0, d0); history.X_update_mode='fixed_antenna_ws';
for t=1:params.T_max
    state.t=t; state.X=X_fixed; state.theta=theta_fixed; state.phi=phi_fixed; state.W=AO_W(params, scene, model, state);
    [state.S,state.swap_flag]=AO_S_fixed(params, scene, model, state); state.X=X_fixed; state.theta=theta_fixed; state.phi=phi_fixed;
    R=Signal_model('sum_rate', params, scene, state, []); [Re,dt]=Effective_rate_model(params, scene, state, []);
    history.R_sum(end+1,1)=R; history.R_eff(end+1,1)=Re; history.T_X(end+1,1)=dt.T_X; history.T_theta(end+1,1)=dt.T_theta; history.T_phi(end+1,1)=dt.T_phi; history.T_rec(end+1,1)=dt.T_rec; history.time_factor(end+1,1)=dt.time_factor;
    history.S_cells{t,1}=state.S; history.X_cells{t,1}=state.X; history.theta_cells{t,1}=state.theta; history.phi_cells{t,1}=state.phi;
    if abs(history.R_eff(end)-history.R_eff(end-1)) < params.eps_outer, break; end
end
end

function [state, history] = run_fixed_antenna_ws_reW_case(params, scene, model, state)
if exist('AO_S_fixed_reW.m','file') ~= 2
    error('AO_S_fixed_reW.m not found. Cannot run Fixed W+S reW scheme.');
end
X_fixed = state.X; theta_fixed = pi*ones(params.N,params.M); phi_fixed = zeros(params.N,params.M);
state.theta=theta_fixed; state.phi=phi_fixed;
[Reff0, d0] = Effective_rate_model(params, scene, state, []);
history = init_history_full(params, scene, state, Reff0, d0); history.X_update_mode='fixed_antenna_ws_reW';
for t=1:params.T_max
    state.t=t; state.X=X_fixed; state.theta=theta_fixed; state.phi=phi_fixed; state.W=AO_W(params, scene, model, state);
    [state.S,state.W,state.swap_flag]=AO_S_fixed_reW(params, scene, model, state);
    state.X=X_fixed; state.theta=theta_fixed; state.phi=phi_fixed;
    R=Signal_model('sum_rate', params, scene, state, []); [Re,dt]=Effective_rate_model(params, scene, state, []);
    history.R_sum(end+1,1)=R; history.R_eff(end+1,1)=Re; history.T_X(end+1,1)=dt.T_X; history.T_theta(end+1,1)=dt.T_theta; history.T_phi(end+1,1)=dt.T_phi; history.T_rec(end+1,1)=dt.T_rec; history.time_factor(end+1,1)=dt.time_factor;
    history.S_cells{t,1} = state.S;
    history.X_cells{t,1} = state.X;
    history.theta_cells{t,1} = state.theta;
    history.phi_cells{t,1} = state.phi;
    history.swap_flag(end+1,1) = state.swap_flag;
    if abs(history.R_eff(end)-history.R_eff(end-1)) < params.eps_outer, break; end
end
end

function [state, history] = run_sa_joint_case(params, scene, model, state)
[state, history] = SA_joint(params, scene, model, state); history = patch_history(history, scene, params, state, 'sa_joint');
end
function [state, history] = run_pso_joint_case(params, scene, model, state)
[state, history] = PSO_joint(params, scene, model, state); history = patch_history(history, scene, params, state, 'pso_joint');
end

function history = patch_history(history, scene, params, state, mode)
fields = {'DEBUG_X_cells','R_after_W','R_after_angle','R_after_X','R_after_S','R_eff_after_W','R_eff_after_angle','R_eff_after_X','R_eff_after_S','R_before_final_W','R_after_final_W'};
for i=1:numel(fields), if ~isfield(history,fields{i}), history.(fields{i}) = []; end, end
if ~isfield(history,'X_update_mode'), history.X_update_mode = mode; end
Rsum = Signal_model('sum_rate', params, scene, state, []); [Reff,dt] = Effective_rate_model(params, scene, state, []);
if ~isfield(history,'R_sum')||isempty(history.R_sum), history.R_sum=Rsum; else, history.R_sum(end,1)=Rsum; end
if ~isfield(history,'R_eff')||isempty(history.R_eff), history.R_eff=Reff; else, history.R_eff(end,1)=Reff; end
if ~isfield(history,'T_X')||isempty(history.T_X), history.T_X=dt.T_X; else, history.T_X(end,1)=dt.T_X; end
if ~isfield(history,'T_theta')||isempty(history.T_theta), history.T_theta=dt.T_theta; else, history.T_theta(end,1)=dt.T_theta; end
if ~isfield(history,'T_phi')||isempty(history.T_phi), history.T_phi=dt.T_phi; else, history.T_phi(end,1)=dt.T_phi; end
if ~isfield(history,'T_rec')||isempty(history.T_rec), history.T_rec=dt.T_rec; else, history.T_rec(end,1)=dt.T_rec; end
if ~isfield(history,'time_factor')||isempty(history.time_factor), history.time_factor=dt.time_factor; else, history.time_factor(end,1)=dt.time_factor; end
end

function h = init_history_full(params, scene, state, Reff0, d0)
r0 = Signal_model('individual_rates', params, scene, state, []);
h = struct(); h.X0=state.X; h.theta0=state.theta; h.phi0=state.phi; h.S0=state.S; h.rates0=r0;
h.R_sum=sum(r0); h.R_eff=Reff0; h.T_X=d0.T_X; h.T_theta=d0.T_theta; h.T_phi=d0.T_phi; h.T_rec=d0.T_rec; h.time_factor=d0.time_factor;
h.R_after_W=[]; h.R_after_angle=[]; h.R_after_X=[]; h.R_after_S=[]; h.R_eff_after_W=[]; h.R_eff_after_angle=[]; h.R_eff_after_X=[]; h.R_eff_after_S=[];
h.S_cells={}; h.X_cells={}; h.theta_cells={}; h.phi_cells={}; h.DEBUG_X_cells={}; h.swap_flag=false;
end

function final_bar_ab = run_final_bar_ab_cases(base_params, schemes, MC, user_pos_pools)
ab_cases = [0.5 0.3; 0.3 0.18]; ns = numel(schemes);
Rsum = zeros(2,ns,MC); Reff = zeros(2,ns,MC);
for i=1:2
    params_ab = base_params; params_ab.a = ab_cases(i,1); params_ab.b = ab_cases(i,2);
    for mc=1:MC
        user_pos_pool = user_pos_pools{mc};
        scene_case = build_scene_with_fixed_users(params_ab,user_pos_pool);
        for s=1:ns
            init_seed_case = params_ab.seed + 20000 + mc;
            alg_seed_case  = params_ab.seed + 30000 + 100*s + mc;
            out = run_one_case(params_ab, schemes(s).init_mode, schemes(s).alg_mode, ...
                init_seed_case, alg_seed_case, scene_case);
            Rsum(i,s,mc)=out.final_R_sum; Reff(i,s,mc)=out.final_R_eff;
        end
    end
end
final_bar_ab.ab_cases = ab_cases;
final_bar_ab.mean_R_sum_ab = squeeze(mean(Rsum,3)); final_bar_ab.std_R_sum_ab = squeeze(std(Rsum,0,3));
final_bar_ab.mean_R_eff_ab = squeeze(mean(Reff,3)); final_bar_ab.std_R_eff_ab = squeeze(std(Reff,0,3));
end

function draw_final_bar_ab(final_bar_ab, schemes)
for i = 1:2
    figure('Name', sprintf('Fig_FinalBar_ab_%d', i), ...
        'Position', [100 100 900 520]);

    Y = [final_bar_ab.mean_R_sum_ab(i,:).', ...
         final_bar_ab.mean_R_eff_ab(i,:).'];

    hb = bar(Y, 'grouped', 'BarWidth', 0.72);
    hold on;

    hb(1).FaceColor = [0.00 0.45 0.74];
    hb(2).FaceColor = [0.85 0.33 0.10];
    hb(1).EdgeColor = 'none';
    hb(2).EdgeColor = 'none';

    xticks(1:numel(schemes));
    xticklabels({schemes.name});
    xtickangle(25);

    ylabel('Rate (bit/s/Hz)');
    title(sprintf('Final performance, a=%.2f, b=%.2f', ...
        final_bar_ab.ab_cases(i,1), final_bar_ab.ab_cases(i,2)));

    legend({'R_{sum}', 'R_{eff}'}, ...
        'Location', 'northoutside', ...
        'Orientation', 'horizontal');

    ymax = max(Y(:));
    ylim([0, 1.15 * ymax]);

    grid on;
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    ax.GridAlpha = 0.18;
    ax.LineWidth = 1.0;
    ax.FontSize = 11;
    box on;

    hold off;
end
end

function H2_ab = draw_H2_ab_cases(base_params)
ab_cases = [0.5 0.3;
            0.3 0.18];

params_h2 = base_params;
params_h2.N = 1;
params_h2.M = 1;
params_h2.K = 1;
params_h2.NRF = 1;
params_h2.K_max = 1;
params_h2.K_serv = 1;

if isfield(params_h2,'waveguide_Dx')
    wg_Dx = params_h2.waveguide_Dx;
else
    wg_Dx = params_h2.Dx;
end
if isfield(params_h2,'waveguide_Dy')
    wg_Dy = params_h2.waveguide_Dy;
else
    wg_Dy = params_h2.Dy;
end
if isfield(params_h2,'area_Dx')
    area_Dx = params_h2.area_Dx;
else
    area_Dx = params_h2.Dx;
end
if isfield(params_h2,'area_Dy')
    area_Dy = params_h2.area_Dy;
else
    area_Dy = params_h2.Dy;
end

x_grid = linspace(0, area_Dx, 81);
y_grid = linspace(0, area_Dy, 81);
z_grid = linspace(0, params_h2.d, 31);

state = struct();
state.X = wg_Dy / 2;
state.theta = pi;
state.phi = 0;

H2_ab = struct();
H2_ab.ab_cases = ab_cases;
H2_ab.x_grid = x_grid;
H2_ab.y_grid = y_grid;
H2_ab.z_grid = z_grid;
H2_ab.state = state;

for ia = 1:size(ab_cases,1)
    params_h2.a = ab_cases(ia,1);
    params_h2.b = ab_cases(ia,2);

    scene = Channel_model('build_scene', params_h2, [], [], []);
    scene.xW = wg_Dx / 2;
    scene.feed_pos = [scene.xW; 0; params_h2.d];
    scene.N = 1;
    scene.M = 1;

    [Yg, Xg] = meshgrid(y_grid, x_grid);
    H3 = zeros(numel(x_grid), numel(y_grid), numel(z_grid));

    for iz = 1:numel(z_grid)
        Zg = z_grid(iz) * ones(size(Xg));
        scene.user_pos = [Xg(:).'; Yg(:).'; Zg(:).'];
        scene.K = numel(Xg);

        extra = struct();
        extra.use_all = true;
        ch_out = Channel_model('all_users', params_h2, scene, state, extra);
        H = ch_out.H;

        pow_map = abs(H).^2;
        H3(:,:,iz) = reshape(pow_map, size(Xg));
    end

    figure('Name', sprintf('Fig_H2_3D_ab_%d', ia), 'Position', [100 100 1100 480]);

    subplot(1,2,1);
    draw_main_lobe_pattern(params_h2, scene, state, area_Dx, area_Dy, H3, x_grid, y_grid, z_grid);

    subplot(1,2,2);
    H2_z0 = H3(:,:,1);
    H2_z0_plot = max(H2_z0, 1e-30);
    imagesc(y_grid, x_grid, H2_z0_plot);
    set(gca, 'YDir', 'normal');
    set(gca, 'ColorScale', 'log');
    hold on;
    plot(state.X, scene.xW, 'w.', 'MarkerSize', 18);
    hold off;
    colorbar;
    xlabel('y (m)');
    ylabel('x (m)');
    title(sprintf('z = 0 plane |H|^2, a=%.2f, b=%.2f', params_h2.a, params_h2.b));

    H2_ab.scene_xW = scene.xW;
    H2_ab.H3{ia} = H3;
    H2_ab.H2_z0{ia} = H2_z0;
end
end


function draw_main_lobe_pattern(params_h2, scene, state, area_Dx, area_Dy, H3, x_grid, y_grid, z_grid)
beam_center = [state.X, scene.xW, 0];
H2_z0 = H3(:,:,1);
P3 = H3 / max(H3(:));
P2 = H2_z0 / max(H3(:));
p_th = 0.10;
color_z_top = 1;
dx = x_grid(2) - x_grid(1);
dy = y_grid(2) - y_grid(1);

r_eq = zeros(numel(z_grid),1);
for iz = 1:numel(z_grid)
    P = H3(:,:,iz);
    Pn = P / max(P(:));
    mask = (Pn >= p_th);
    A = nnz(mask) * dx * dy;
    r_eq(iz) = sqrt(A / pi);
end
r_eq = smoothdata(r_eq, 'movmean', 3);

idx_lower = find(z_grid <= color_z_top);
z_lower = z_grid(idx_lower);
r_profile_lower = r_eq(idx_lower);

[Yg, Xg] = meshgrid(y_grid, x_grid);
R2 = sqrt((Yg - state.X).^2 + (Xg - scene.xW).^2);
r_max = max(r_profile_lower);
r_bins = linspace(0, r_max, 120);
r_bin_center = 0.5 * (r_bins(1:end-1) + r_bins(2:end));
p_radial = zeros(numel(r_bin_center),1);
for ir = 1:numel(r_bin_center)
    ring_mask = (R2 >= r_bins(ir)) & (R2 < r_bins(ir+1));
    if nnz(ring_mask) > 0
        p_radial(ir) = mean(P2(ring_mask));
    else
        p_radial(ir) = NaN;
    end
end
p_radial = fillmissing(p_radial, 'nearest');
p_radial = smoothdata(p_radial, 'movmean', 5);
c_lower = interp1(r_bin_center, p_radial, r_profile_lower, 'linear', 'extrap');

z_upper = 2*color_z_top - z_lower(end-1:-1:1);
r_profile_upper = r_profile_lower(end-1:-1:1);
c_outer = c_lower(end);
c_upper = linspace(c_outer, 0.85*c_outer, numel(z_upper)).';

z_profile = [z_lower(:); z_upper(:)];
r_profile = [r_profile_lower(:); r_profile_upper(:)];
c_profile = [c_lower(:); c_upper(:)];

theta = linspace(0,2*pi,80);
[U, Z] = meshgrid(theta, z_profile);
R = repmat(r_profile, 1, numel(theta));
Y = state.X + R .* cos(U);
X = scene.xW + R .* sin(U);
C = repmat(c_profile, 1, numel(theta));
pa_pos = [state.X, scene.xW, z_profile(end)];

surf(Y, X, Z, C, 'EdgeColor', 'none', 'FaceAlpha', 0.90);
hold on;

foot_r = r_profile_lower(1);
foot_rx = 1.05 * foot_r;
foot_ry = 0.85 * foot_r;
theta_fp = linspace(0,2*pi,180);
Yf = beam_center(1) + foot_ry*cos(theta_fp);
Xf = beam_center(2) + foot_rx*sin(theta_fp);
Zf = zeros(size(theta_fp));
fill3(Yf, Xf, Zf, [0.4 0.8 1.0], 'FaceAlpha', 0.25, 'EdgeColor', [0.2 0.5 0.9], 'LineStyle', '--');

plot3([pa_pos(1), beam_center(1)], [pa_pos(2), beam_center(2)], [pa_pos(3), beam_center(3)], 'k--', 'LineWidth', 1.2);
plot3(pa_pos(1), pa_pos(2), pa_pos(3), 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

xlabel('y (m)');
ylabel('x (m)');
zlabel('z (m)');
title(sprintf('3D main lobe pattern, a=%.2f, b=%.2f', params_h2.a, params_h2.b));
colormap(jet);
grid on;
axis tight;
daspect([1 1 0.6]);
view(45,25);
box on;
hold off;
end

function conv_results = run_convergence_cases(base_params, schemes)
params_conv = base_params;
params_conv.T_max = 30;
params_conv.SA_max_iter = 5000;
ns=numel(schemes); conv_results=struct('name',cell(ns,1),'alg_mode',cell(ns,1),'R_eff',cell(ns,1),'T_max',cell(ns,1),'SA_max_iter',cell(ns,1));
scene_case=build_scene_with_fixed_users(params_conv, build_fixed_user_pool(params_conv,1,'conv',params_conv.seed+50001));
for s=1:ns
out=run_one_case(params_conv,schemes(s).init_mode,schemes(s).alg_mode,params_conv.seed+1,params_conv.seed+100+s,scene_case);
conv_results(s).name = schemes(s).name;
conv_results(s).alg_mode = schemes(s).alg_mode;
conv_results(s).R_eff = out.history.R_eff(:);
conv_results(s).T_max = params_conv.T_max;
conv_results(s).SA_max_iter = params_conv.SA_max_iter;
end
end
function rate_cells = collect_rate_cdf_data(base_params, schemes, MC, user_pos_pools)
ns=numel(schemes); rate_cells=cell(ns,1);
for mc=1:MC
scene_case=build_scene_with_fixed_users(base_params, user_pos_pools{mc});
for s=1:ns, init_seed_case = base_params.seed + 20000 + mc; alg_seed_case  = base_params.seed + 30000 + 100*s + mc; out = run_one_case(base_params, schemes(s).init_mode, schemes(s).alg_mode, init_seed_case, alg_seed_case, scene_case); rate_cells{s}=[rate_cells{s}; out.rates_final(:)]; end
end
end

function user_pos_pool = build_fixed_user_pool(base_params, x_vec, sweep_type, scene_seed)
params_pool=base_params; if strcmp(sweep_type,'K'), params_pool.K=max(x_vec); else, params_pool.K=base_params.K; end
rng(scene_seed); scene_pool=Channel_model('build_scene', params_pool, [], [], []); user_pos_pool=scene_pool.user_pos;
end
function scene_case = build_scene_with_fixed_users(params_case, user_pos_pool)
scene_case = Channel_model('build_scene', params_case, [], [], []); scene_case.user_pos = user_pos_pool(:,1:params_case.K); scene_case.K=params_case.K; scene_case.M=params_case.M; scene_case.N=params_case.N;
end
function draw_mean_error_curve(x_vec, mean_R, std_R, schemes, x_label_text, title_text)
for s=1:numel(schemes), plot(x_vec,mean_R(:,s),'-o','LineWidth',1.4,'MarkerSize',5); hold on; end
xlabel(x_label_text); ylabel('Average effective spectral efficiency (bit/s/Hz)'); title(title_text,'FontSize',11);
legend({schemes.name},'Location','southoutside','NumColumns',2,'FontSize',8); grid on; set(gca,'FontSize',10);
end
function draw_convergence(conv_results, schemes)
break_iter = 30; x_end_real = 5000; x_end_plot = 5000; x_break_plot = x_end_plot/3;
figure('Name','Fig5_Convergence_BrokenAxis','Position',[100 100 1100 560]);
for s=1:numel(conv_results)
    r = conv_results(s).R_eff(:);
    if strcmp(conv_results(s).alg_mode,'SA_joint')
        if numel(r) == conv_results(s).SA_max_iter + 1, x_real = (0:conv_results(s).SA_max_iter).';
        elseif numel(r) == conv_results(s).SA_max_iter, x_real = (1:conv_results(s).SA_max_iter).';
        else, x_real = round(linspace(0, conv_results(s).SA_max_iter, numel(r))).'; end
    else
        x_real = (0:numel(r)-1).';
    end
    [r_best, idx_best] = max(r);
    x_real_plot = [x_real(1:idx_best); x_end_real];
    r_plot = [r(1:idx_best); r_best];
    plot(compress_conv_x(x_real_plot, break_iter, x_break_plot, x_end_real, x_end_plot), r_plot, '-o', 'LineWidth', 1.6, 'MarkerSize', 6, 'MarkerFaceColor', 'none'); hold on;
end
y_all = cell2mat(arrayfun(@(s) s.R_eff(:), conv_results, 'UniformOutput', false));
r_min = min(y_all); r_max = max(y_all); pad = max(1e-6, 0.08*(r_max-r_min)); ylim([r_min-pad, r_max+pad]);
yl = ylim; plot([x_break_plot x_break_plot], yl, 'k--', 'LineWidth', 1.2);
text(x_break_plot + 80, yl(1) + 0.08*(yl(2)-yl(1)), 'x-axis compressed after 20 iterations', 'FontSize', 11);
tick_real = [0 5 10 15 20 25 30 500 1000 1500 2000 2500 3000 3500 4000 4500 5000];
xticks(compress_conv_x(tick_real, break_iter, x_break_plot, x_end_real, x_end_plot)); xticklabels(string(tick_real));
legend({conv_results.name},'Location','northeastoutside','FontSize',9);
xlabel('Iteration index'); ylabel('R_{eff} (bit/s/Hz)'); title('Convergence behavior of different schemes with compressed x-axis');
grid on; box on; set(gca,'FontSize',10);
end
function x_plot = compress_conv_x(x_real, break_iter, x_break_plot, x_end_real, x_end_plot)
x_real = x_real(:); x_plot = zeros(size(x_real)); idx = x_real <= break_iter;
x_plot(idx) = x_real(idx) / break_iter * x_break_plot;
x_plot(~idx) = x_break_plot + (x_real(~idx) - break_iter) / (x_end_real - break_iter) * (x_end_plot - x_break_plot);
end
function draw_rate_cdf(rate_cells, schemes)
for s=1:numel(schemes), r=sort(rate_cells{s}(:)); F=(1:numel(r))/numel(r); plot(r,F,'LineWidth',1.2); hold on; end
xlabel('Per-user rate (bit/s/Hz)'); ylabel('CDF'); title('CDF of Per-user Rate','FontSize',11); legend({schemes.name},'Location','southoutside','NumColumns',2,'FontSize',8); grid on; set(gca,'FontSize',10);
end
function draw_geometry_case(geo_result)
params_case=geo_result.params; scene=geo_result.scene; state=geo_result.state; history=geo_result.history; M=params_case.M; user_pos=scene.user_pos; S=state.S;
if isfield(params_case,'waveguide_Dy'), wg_Dy=params_case.waveguide_Dy; else, wg_Dy=params_case.Dy; end
if isfield(params_case,'area_Dx'), area_Dx=params_case.area_Dx; else, area_Dx=params_case.Dx; end
if isfield(params_case,'area_Dy'), area_Dy=params_case.area_Dy; else, area_Dy=params_case.Dy; end
scatter(user_pos(1,:), user_pos(2,:), 25, 'filled'); hold on; scatter(user_pos(1,S), user_pos(2,S), 70);
for n=1:params_case.N, line([scene.xW(n),scene.xW(n)],[0,wg_Dy]); plot(scene.xW(n)*ones(1,M),history.X0(n,:),'o'); plot(scene.xW(n)*ones(1,M),state.X(n,:),'x'); end
x_pa = repmat(scene.xW(:),1,M); x_pa=x_pa(:); y_pa=state.X(:); u=sin(state.theta(:)).*cos(state.phi(:)); v=sin(state.theta(:)).*sin(state.phi(:)); nm=sqrt(u.^2+v.^2); quiver(x_pa,y_pa,u./(nm+eps),v./(nm+eps),0.6,'LineWidth',0.8);
xlim([0 area_Dx]); ylim([0 area_Dy]); xlabel('x (m)'); ylabel('y (m)'); title('Final PA/User Configuration','FontSize',11);
legend({'All users','Served users','Waveguide','Initial PA','Final PA','PA orientation'},'Location','eastoutside','FontSize',8); grid on; set(gca,'FontSize',10);
end
function sweep_id = get_sweep_id(sweep_type)
if strcmp(sweep_type,'snr'), sweep_id=1; elseif strcmp(sweep_type,'K'), sweep_id=2; elseif strcmp(sweep_type,'N'), sweep_id=3; elseif strcmp(sweep_type,'M'), sweep_id=4; elseif strcmp(sweep_type,'Dy'), sweep_id=8; else, sweep_id=9; end
end
