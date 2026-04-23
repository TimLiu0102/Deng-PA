function test_reproduce_fig2_channel_model()
% test_reproduce_fig2_channel_model：基于现有 Channel_model 入口复现 Figure 2 风格测试

%% 第1部分：清空环境与参数设置
clc; clear; close all;

params = struct();

% 1) 系统规模参数：单波导、单PA、单链路
params.N = 1;
params.M = 1;
params.K = 1;
params.NRF = 1;
params.K_max = 1;
params.K_serv = 1;

% 2) 几何参数
params.Dx = 10;
params.Dy = 10;
params.d = 5;
params.Delta = 0.5;

% 3) 信道参数（按测试要求设置）
params.lambda = 0.01;
params.n_eff = 1.6;
params.alphaW = 0.01;
params.alphaL = 0.96;
params.v = 1.1;
params.n_refr = 1.5;
params.eta = params.lambda^2 / (4*pi);
params.P_max = 1.0;
params.sigma2 = 1e-9;

% 仅用于 build_scene 中随机用户生成，不影响后续手动覆盖网格用户点
params.seed = 7;
rng(params.seed);

% 先调用 build_scene，保持与主工程一致的场景构造入口
scene = Channel_model('build_scene', params, [], [], []);

% 这里手动覆盖单PA测试几何：
% 目的不是改模型，而是把测试场景固定成“论文 Figure 2 风格的单PA可视化设置”
scene.xW = 2;
scene.feed_pos = [2; 0; params.d];
scene.N = 1;
scene.M = 1;

% 状态变量：固定 PA 在波导上的 y 位置与方位角，俯仰角后续循环
state = struct();
state.X = 2;
state.phi = pi/2;

% 三组横截面尺寸与三个俯仰角
ab_list = [5, 3; 10, 6; 20, 12];
theta_list = [pi, 5*pi/6, 2*pi/3];

%% 第2部分：复现 Figure 2(a) 风格——z=0 平面 |H|^2 空间分布
% 说明：把平面网格点直接当成“用户位置”，通过 all_users 一次性获得每个点的 H
% 在 N=M=1 时，H 的尺寸为 1 x Kgrid，因此 abs(H).^2 就是每个网格点的功率
x_grid = linspace(0, 10, 241);
y_grid = linspace(0, 10, 241);
[Yg, Xg] = meshgrid(y_grid, x_grid);
Kgrid = numel(Xg);

figure('Name', 'Fig2a_style_z0_power_map', 'Color', 'w');
for ia = 1:size(ab_list,1)
    params.a = ab_list(ia,1);
    params.b = ab_list(ia,2);

    for it = 1:numel(theta_list)
        state.theta = theta_list(it);

        scene.user_pos = [Xg(:).'; Yg(:).'; zeros(1, Kgrid)];
        scene.K = Kgrid;

        extra = struct();
        extra.use_all = true;
        ch_out = Channel_model('all_users', params, scene, state, extra);
        H = ch_out.H;

        pow_map = abs(H).^2;
        pow_map_2d = reshape(pow_map, size(Xg));

        idx = (ia-1)*numel(theta_list) + it;
        subplot(3, 3, idx);
        imagesc(y_grid, x_grid, pow_map_2d);
        set(gca, 'YDir', 'normal');
        set(gca, 'ColorScale', 'log');
        hold on;
        plot(2, 2, 'w.', 'MarkerSize', 18);
        hold off;
        colorbar;
        xlabel('y');
        ylabel('x');
        title(sprintf('a=%g, b=%g, \\theta=%.4f rad', params.a, params.b, state.theta));
    end
end
sgtitle('基于 Channel\_model 的 z=0 平面 |H|^2 分布');

%% 第3部分：复现 Figure 2(b) 风格——沿波导轴线 |H|^2 曲线
% 固定 x=scene.xW=2, z=0，仅扫描 y，比较不同 theta 下主瓣峰值与宽度变化
y_line = linspace(0, 10, 1201);
x_line = scene.xW * ones(size(y_line));
z_line = zeros(size(y_line));

figure('Name', 'Fig2b_style_axis_line', 'Color', 'w');
for ia = 1:size(ab_list,1)
    params.a = ab_list(ia,1);
    params.b = ab_list(ia,2);

    subplot(3,1,ia);
    hold on;

    for it = 1:numel(theta_list)
        state.theta = theta_list(it);

        scene.user_pos = [x_line; y_line; z_line];
        scene.K = numel(y_line);

        extra = struct();
        extra.use_all = true;
        ch_out = Channel_model('all_users', params, scene, state, extra);
        H = ch_out.H;

        pow_line = abs(H).^2;
        semilogy(y_line, pow_line, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('\\theta = %.4f rad', state.theta));
    end

    hold off;
    grid on;
    xlabel('y');
    ylabel('|H|^2');
    title(sprintf('沿波导轴线功率曲线（a=%g, b=%g）', params.a, params.b));
    legend('Location', 'best');
end
sgtitle('基于 Channel\_model 的沿波导轴线 |H|^2 曲线');

%% 第4部分：前向半空间诊断图（仅排查用）
% 说明：这里不画信道功率，而是直接调用 global_to_local 检查局部坐标 y~ 的符号分布
% 若 y~ > 0 的区域与预期前向半空间不符，通常说明坐标变换或角度定义存在问题
params.a = 5;
params.b = 3;

figure('Name', 'Forward_halfspace_diagnosis', 'Color', 'w');
for it = 1:numel(theta_list)
    state.theta = theta_list(it);

    mask_forward = false(size(Xg));
    for k = 1:Kgrid
        q_k = [Xg(k); Yg(k); 0];
        p_nm = [scene.xW; state.X; params.d];

        extra_loc = struct();
        extra_loc.q_k = q_k;
        extra_loc.p_nm = p_nm;
        extra_loc.theta_nm = state.theta;
        extra_loc.phi_nm = state.phi;

        local_xyz = Channel_model('global_to_local', params, scene, state, extra_loc);
        mask_forward(k) = (local_xyz(2) > 0);
    end

    subplot(1,3,it);
    imagesc(y_grid, x_grid, double(mask_forward));
    set(gca, 'YDir', 'normal');
    colormap(gca, gray);
    colorbar;
    xlabel('y');
    ylabel('x');
    title(sprintf('\\theta=%.4f rad 时 y~ > 0 区域', state.theta));
end
sgtitle('前向半空间诊断（仅排查用）');

end
