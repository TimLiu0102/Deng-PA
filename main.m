function main()
% 论文复现主流程：参数设置 -> 建模 -> 初始化 -> AO循环 -> 输出绘图

%% 1) 参数设置（集中在此处）
params.N = 3; params.M = 4; params.K = 24; params.NRF = 12;
params.K_serv = 6; params.K_max = 24;

params.Dx = 8; params.Dy = 20; params.d = 3; params.Delta = 0.8;
params.lambda = 0.01; params.n_eff = 1.6; params.alphaW = 0.06; params.alphaL = 0.98;
params.a = 0.8; params.b = 0.6; params.v = 1.0; params.n_refr = 1.5;
params.B_norm = 1.0; params.eta = 1.0;

params.P_max = 1; params.sigma2 = 1e-3;
params.lambda_mov = 0.05;

params.I_W = 40; params.eps_W = 1e-4;
params.I_theta = 6; params.Delta_theta0 = 0.08; params.Delta_phi0 = 0.08;
params.beta_theta = 0.6; params.beta_phi = 0.6; params.eps_theta = 1e-5;

params.step_fd = 1e-3; params.line_search_alpha0 = 0.5;
params.line_search_beta = 0.5; params.line_search_max_iter = 8;
params.eps_X = 1e-5;

params.T_S = 2; params.L_in = 2; params.L_out = 4; params.eps_S = 1e-5;
params.T_max = 30; params.eps_outer = 1e-4;

params.user_x_rng = [0, params.Dx];
params.user_y_rng = [0, params.Dy];
params.user_z_rng = [0, 2];

rng(7);

%% 2) 场景与信道几何
scene = Channel_model('build_scene', [], [], [], [], [], params);
params.scene = scene;

%% 3) 问题定义（仅定义，不求解）
problem = Problem_formulation(params);

%% 4) 初始化（不初始化W）
init = Initialization(scene, params);
S = init.S; X = init.X; theta = init.theta; phi = init.phi;
C = init.C; Emax = init.Emax;

%% 5) 初始sum rate（真实目标）
W = [];
R0 = Signal_model('sum_rate', S, X, W, theta, phi, scene, params);
history.R = R0; history.swap = false;

%% 6) 外层AO循环：AO_W -> AO_angle -> AO_X -> AO_S
for t = 1:params.T_max
    W = AO_W(S, X, theta, phi, W, scene, params);
    [theta, phi] = AO_angle(S, X, W, theta, phi, scene, params);
    X = AO_X(S, X, W, theta, phi, scene, params);
    [S, swap_flag] = AO_S(t, S, C, Emax, X, W, theta, phi, scene, params);

    R_new = Signal_model('sum_rate', S, X, W, theta, phi, scene, params);
    history.R(end+1,1) = R_new; %#ok<AGROW>
    history.swap(end+1,1) = swap_flag; %#ok<AGROW>

    % 非下降更新 + 外层停止：论文里各子块候选均按非下降准则接受，故目标单调非减
    if abs(history.R(end) - history.R(end-1)) < params.eps_outer
        break;
    end
end

%% 7) 打印与绘图
Print_and_Plot(history, S, X, theta, phi, problem, scene, params);

end
