function Print_and_Plot(params, scene, model, result)
% Print_and_Plot：仅负责结果打印与基础可视化

state = result.state;
history = result.history;

%% ======================== 命令行打印 ========================
fprintf('\n================ 结果汇总 ================\n');
fprintf('系统规模: N=%d, M=%d, K=%d\n', params.N, params.M, params.K);
fprintf('服务用户数 K_serv=%d\n', model.K_serv);

if isfield(history,'R_sum') && ~isempty(history.R_sum)
    R_hist = history.R_sum(:);
    fprintf('初始 sum rate: %.6f\n', R_hist(1));
    fprintf('最终 sum rate: %.6f\n', R_hist(end));
    fprintf('外层迭代次数: %d\n', numel(R_hist)-1);
else
    fprintf('sum rate 历史缺失，无法打印收敛信息。\n');
end

if isfield(state,'S')
    fprintf('最终服务用户集合 S: ');
    fprintf('%d ', state.S);
    fprintf('\n');
end

if isfield(state,'W') && ~isempty(state.W)
    final_power = real(trace(state.W * state.W'));
    fprintf('最终总发射功率 tr(WW^H): %.6f\n', final_power);
    fprintf('是否满足功率约束 (<= P_max): %d\n', final_power <= params.P_max + 1e-10);
else
    fprintf('W 不存在，无法打印最终功率。\n');
end

if isfield(state,'S')
    fprintf('最终服务用户数是否等于 K_serv: %d\n', numel(state.S) == model.K_serv);
end

if isfield(history,'swap_flag') && ~isempty(history.swap_flag)
    idx_swap = find(history.swap_flag(:));
    if isempty(idx_swap)
        fprintf('未发生用户交换。\n');
    else
        fprintf('发生用户交换的迭代轮次: ');
        fprintf('%d ', idx_swap);
        fprintf('\n');
    end
end

%% ======================== 图1：sum rate 收敛曲线 ========================
if isfield(history,'R_sum') && ~isempty(history.R_sum)
    figure;
    plot(0:numel(history.R_sum)-1, history.R_sum, '-o');
    xlabel('外层迭代次数');
    ylabel('sum rate');
    title('外层迭代 sum rate 收敛曲线');
    grid on;
end

%% ======================== 图2：系统几何示意图 ========================
if isfield(scene,'user_pos') && isfield(scene,'xW') && isfield(state,'X')
    figure;
    hold on;

    user_pos = scene.user_pos;
    scatter(user_pos(1,:), user_pos(2,:), 25, 'filled');

    N = size(state.X,1);
    M = size(state.X,2);
    for n = 1:N
        % 波导位置线（x固定，y从0到Dy）
        line([scene.xW(n), scene.xW(n)], [0, params.Dy]);

        % 最终PA位置投影
        scatter(scene.xW(n)*ones(1,M), state.X(n,:), 40);
    end

    if isfield(state,'S') && ~isempty(state.S)
        s_idx = state.S(:).';
        scatter(user_pos(1,s_idx), user_pos(2,s_idx), 60);
    end

    xlabel('x 方向位置');
    ylabel('y 方向位置');
    title('系统几何示意图（用户、波导、PA）');
    legend('全部用户','波导','最终PA位置','最终服务用户','Location','best');
    grid on;
    hold off;
end

%% ======================== 图3：初始与最终 PA 位置对比 ========================
has_init_X = false;
if isfield(history,'X') && ~isempty(history.X)
    X_init = history.X{1};
    has_init_X = true;
elseif isfield(state,'X_init')
    X_init = state.X_init;
    has_init_X = true;
end

if has_init_X && isfield(state,'X')
    figure;
    hold on;
    N = size(state.X,1);
    M = size(state.X,2);
    for n = 1:N
        plot(scene.xW(n)*ones(1,M), X_init(n,:), 'o');
        plot(scene.xW(n)*ones(1,M), state.X(n,:), 'x');
    end
    xlabel('x 方向位置');
    ylabel('y 方向位置');
    title('初始与最终 PA 位置对比');
    legend('初始PA位置','最终PA位置','Location','best');
    grid on;
    hold off;
end

%% ======================== 图4：最终角度分布 ========================
if isfield(state,'theta')
    theta_vec = state.theta(:);
    figure;
    plot(1:numel(theta_vec), theta_vec, '-o');
    xlabel('PA 索引');
    ylabel('\theta (rad)');
    title('最终 \theta 分布');
    grid on;
end

if isfield(state,'phi')
    phi_vec = state.phi(:);
    figure;
    plot(1:numel(phi_vec), phi_vec, '-o');
    xlabel('PA 索引');
    ylabel('\phi (rad)');
    title('最终 \phi 分布');
    grid on;
end

%% ======================== 图5：服务用户选择结果 ========================
if isfield(scene,'user_pos') && isfield(state,'S') && ~isempty(state.S)
    figure;
    hold on;
    user_pos = scene.user_pos;
    scatter(user_pos(1,:), user_pos(2,:), 20, 'filled');

    s_idx = state.S(:).';
    scatter(user_pos(1,s_idx), user_pos(2,s_idx), 70);
    for i = 1:numel(s_idx)
        k = s_idx(i);
        text(user_pos(1,k), user_pos(2,k), sprintf(' %d', k));
    end

    xlabel('x 方向位置');
    ylabel('y 方向位置');
    title('最终服务用户选择结果');
    legend('全部用户','服务用户','Location','best');
    grid on;
    hold off;
end

end
