function test_reproduce_fig2_channel_model_v2()
% test_reproduce_fig2_channel_model_v2：在 v1 基础上针对测试盲区做增强

%% 第1部分：清空环境与单波导单PA参数设置
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

% 3) 信道参数
params.lambda = 0.01;
params.n_eff = 1.6;
params.alphaW = 0.01;
params.alphaL = 0.96;
params.v = 1.1;
params.n_refr = 1.5;
params.eta = params.lambda^2 / (4*pi);
params.P_max = 1.0;
params.sigma2 = 1e-9;

rng(7);

% 保持与工程一致：先调用 build_scene，再按测试目的手动覆盖单PA位置
scene = Channel_model('build_scene', params, [], [], []);
scene.xW = 2;
scene.feed_pos = [2; 0; params.d];
scene.N = 1;
scene.M = 1;

state = struct();
state.X = 2;

%% 第2部分：统一色标的 Figure 2(a) 风格图
% 这里统一 caxis 的目的：避免 9 个子图各自自动拉伸颜色，导致不能直接比较绝对峰值强弱
ab_list = [5 3; 10 6; 20 12];
theta_list = [pi, 5*pi/6, 2*pi/3];
state.phi = pi/2;

x_grid = linspace(0, 10, 241);
y_grid = linspace(0, 10, 241);
[Yg, Xg] = meshgrid(y_grid, x_grid);
Kgrid = numel(Xg);

pow_maps = cell(size(ab_list,1), numel(theta_list));

% 第一步：先算完 9 组功率图并缓存，不一边算一边画
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
        pow_maps{ia,it} = reshape(pow_map, size(Xg));
    end
end

% 第二步：收集 9 张图中的全部正值，求统一全局色标范围
all_positive = [];
for ia = 1:size(ab_list,1)
    for it = 1:numel(theta_list)
        tmp = pow_maps{ia,it};
        tmp = tmp(tmp > 0);
        all_positive = [all_positive; tmp(:)]; %#ok<AGROW>
    end
end
cmin = min(all_positive);
cmax = max(all_positive);

% 第三步：统一 caxis 画图
figure('Name', 'Unified_caxis_Fig2a_style', 'Color', 'w');
for ia = 1:size(ab_list,1)
    params.a = ab_list(ia,1);
    params.b = ab_list(ia,2);

    for it = 1:numel(theta_list)
        state.theta = theta_list(it);
        pow_map_2d = pow_maps{ia,it};

        subplot(3,3,(ia-1)*numel(theta_list)+it);
        imagesc(y_grid, x_grid, pow_map_2d);
        set(gca,'YDir','normal');
        set(gca,'ColorScale','log');
        caxis([cmin, cmax]);
        colorbar;
        hold on;
        % PA 投影点
        plot(2, 2, 'wo', 'MarkerSize', 6, 'LineWidth', 1.2);

        % 理论主瓣落点：phi=pi/2 时与 z=0 平面交点
        x_hit = 2;
        y_hit = state.X - params.d * tan(state.theta);
        if y_hit >= 0 && y_hit <= 10
            plot(y_hit, x_hit, 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
        hold off;

        xlabel('y');
        ylabel('x');
        title(sprintf('a=%g, b=%g, \\theta=%.4f rad', params.a, params.b, state.theta));
    end
end
sgtitle('统一色标下的 z=0 平面 |H|^2 分布');

%% 第3部分：扩大窗口的前向/背向诊断图
% 之前窗口只看 [0,10]x[0,10]，可能恰好都落在 y~>0 区域，整图全灰不代表模型一定没问题
% 所以这里扩大到 [-10,10]x[-10,10] 直接检查 y~ 正负分区
params.a = 5;
params.b = 3;

x_grid2 = linspace(-10, 10, 301);
y_grid2 = linspace(-10, 10, 301);
[Yg2, Xg2] = meshgrid(y_grid2, x_grid2);
Kgrid2 = numel(Xg2);

theta_case = [5*pi/6, 5*pi/6];
phi_case = [pi/2, 0];

mask_forward_cases = cell(1,2);
mask_backward_cases = cell(1,2);

for ic = 1:2
    state.theta = theta_case(ic);
    state.phi = phi_case(ic);

    mask_forward = false(size(Xg2));
    mask_backward = false(size(Xg2));

    for k = 1:Kgrid2
        q_k = [Xg2(k); Yg2(k); 0];
        p_nm = [scene.xW; state.X; params.d];

        extra_loc = struct();
        extra_loc.q_k = q_k;
        extra_loc.p_nm = p_nm;
        extra_loc.theta_nm = state.theta;
        extra_loc.phi_nm = state.phi;

        local_xyz = Channel_model('global_to_local', params, scene, state, extra_loc);
        mask_forward(k) = (local_xyz(2) > 0);
        mask_backward(k) = (local_xyz(2) <= 0);
    end

    mask_forward_cases{ic} = mask_forward;
    mask_backward_cases{ic} = mask_backward;
end

figure('Name', 'Expanded_forward_backward_diagnosis', 'Color', 'w');
for ic = 1:2
    subplot(2,2,ic);
    imagesc(y_grid2, x_grid2, double(mask_forward_cases{ic}));
    set(gca,'YDir','normal');
    colorbar;
    xlabel('y');
    ylabel('x');
    title(sprintf('y~ > 0, \\theta=%.4f, \\phi=%.4f', theta_case(ic), phi_case(ic)));

    subplot(2,2,2+ic);
    imagesc(y_grid2, x_grid2, double(mask_backward_cases{ic}));
    set(gca,'YDir','normal');
    colorbar;
    xlabel('y');
    ylabel('x');
    title(sprintf('y~ <= 0, \\theta=%.4f, \\phi=%.4f', theta_case(ic), phi_case(ic)));
end
sgtitle('扩大窗口后的前向/背向半空间诊断');

%% 第4部分：前向功率图 / 背向功率图 分离显示
% 不修改 Channel_model，仅通过 global_to_local 先做掩膜，再用 all_users 得到同一组网格的真实信道功率
% 这样就可以直接观察背向区域是否也有明显能量
pow_front_cases = cell(1,2);
pow_back_cases = cell(1,2);

for ic = 1:2
    state.theta = theta_case(ic);
    state.phi = phi_case(ic);

    scene.user_pos = [Xg2(:).'; Yg2(:).'; zeros(1,Kgrid2)];
    scene.K = Kgrid2;

    extra = struct();
    extra.use_all = true;
    ch_out = Channel_model('all_users', params, scene, state, extra);
    H = ch_out.H;
    pow_map2 = reshape(abs(H).^2, size(Xg2));

    mask_forward = mask_forward_cases{ic};
    mask_backward = mask_backward_cases{ic};

    pow_front = pow_map2;
    pow_front(~mask_forward) = NaN;

    pow_back = pow_map2;
    pow_back(~mask_backward) = NaN;

    pow_front_cases{ic} = pow_front;
    pow_back_cases{ic} = pow_back;
end

% 统一色标（前向+背向四张图共享）
all_positive_fb = [];
for ic = 1:2
    t1 = pow_front_cases{ic};
    t1 = t1(isfinite(t1) & t1 > 0);
    all_positive_fb = [all_positive_fb; t1(:)]; %#ok<AGROW>

    t2 = pow_back_cases{ic};
    t2 = t2(isfinite(t2) & t2 > 0);
    all_positive_fb = [all_positive_fb; t2(:)]; %#ok<AGROW>
end
cmin_fb = min(all_positive_fb);
cmax_fb = max(all_positive_fb);

figure('Name', 'Front_Back_Power_Separated', 'Color', 'w');
for ic = 1:2
    subplot(2,2,ic);
    imagesc(y_grid2, x_grid2, pow_front_cases{ic});
    set(gca,'YDir','normal');
    set(gca,'ColorScale','log');
    caxis([cmin_fb, cmax_fb]);
    colorbar;
    hold on;
    plot(2,2,'wo','MarkerSize',6,'LineWidth',1.2);
    hold off;
    xlabel('y');
    ylabel('x');
    title(sprintf('前向功率图, \\theta=%.4f, \\phi=%.4f', theta_case(ic), phi_case(ic)));

    subplot(2,2,2+ic);
    imagesc(y_grid2, x_grid2, pow_back_cases{ic});
    set(gca,'YDir','normal');
    set(gca,'ColorScale','log');
    caxis([cmin_fb, cmax_fb]);
    colorbar;
    hold on;
    plot(2,2,'wo','MarkerSize',6,'LineWidth',1.2);
    hold off;
    xlabel('y');
    ylabel('x');
    title(sprintf('背向功率图, \\theta=%.4f, \\phi=%.4f', theta_case(ic), phi_case(ic)));
end
sgtitle('前向区域与背向区域的 |H|^2 分离显示');

%% 第5部分：最关键的定量测试——波束轴正向/反向 3D 采样对比
% 这是最关键测试：直接沿同一束轴方向做正向/反向 3D 采样并比较峰值
% 不需要改 Channel_model，只要把采样点塞进 scene.user_pos，再调用 all_users 即可
% 若反向轴线持续出现明显功率，就能量化暴露潜在背向泄漏
params.a = 5;
params.b = 3;
state.theta = 5*pi/6;
state.phi = 0;

p_nm = [scene.xW; state.X; params.d];
dir_vec = [sin(state.theta)*cos(state.phi); ...
           sin(state.theta)*sin(state.phi); ...
           cos(state.theta)];

s_line = linspace(0.2, 8, 400);
q_forward = p_nm + dir_vec * s_line;
q_backward = p_nm - dir_vec * s_line;

scene.user_pos = q_forward;
scene.K = numel(s_line);
extra = struct();
extra.use_all = true;
ch_fwd = Channel_model('all_users', params, scene, state, extra);
H_forward = ch_fwd.H;
pow_forward = abs(H_forward).^2;

scene.user_pos = q_backward;
scene.K = numel(s_line);
ch_bwd = Channel_model('all_users', params, scene, state, extra);
H_backward = ch_bwd.H;
pow_backward = abs(H_backward).^2;

figure('Name', '3D_axis_forward_backward_compare', 'Color', 'w');
semilogy(s_line, pow_forward, 'b-', 'LineWidth', 1.3, 'DisplayName', '正向轴线');
hold on;
semilogy(s_line, pow_backward, 'r--', 'LineWidth', 1.3, 'DisplayName', '反向轴线');
hold off;
grid on;
xlabel('s');
ylabel('|H|^2');
title('波束轴正向/反向 3D 采样对比');
legend('Location', 'best');

ratio_back_to_front = max(pow_backward) / max(pow_forward);
fprintf('\n===== 波束轴正向/反向 3D 采样定量结果 =====\n');
fprintf('max(pow_forward) = %.6e\n', max(pow_forward));
fprintf('max(pow_backward) = %.6e\n', max(pow_backward));
fprintf('ratio_back_to_front = %.6e\n', ratio_back_to_front);

figure('Name', 'Backward_to_Forward_ratio_vs_s', 'Color', 'w');
semilogy(s_line, pow_backward ./ pow_forward, 'k-', 'LineWidth', 1.3);
grid on;
xlabel('s');
ylabel('pow\_backward / pow\_forward');
title('波束轴正向/反向 3D 采样对比（比值曲线）');

end
