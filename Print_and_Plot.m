function Print_and_Plot(params, scene, model, result)
% Print_and_Plot：仅负责结果打印与基础可视化（不改变算法流程）

state = result.state;
history = result.history;

% 约定：X(n,m) 表示第 n 条波导上第 m 个 PA 的 y 坐标

%% ======================== 最终用户指标（可选） ========================
has_full_state = isfield(state,'S') && isfield(state,'X') && isfield(state,'theta') ...
    && isfield(state,'phi') && isfield(state,'W') && ~isempty(state.S);

rates_final = [];
sinr_final = [];
if has_full_state
    rates_final = Signal_model('individual_rates', params, scene, state, []);
    s_idx = state.S(:).';
    sinr_final = zeros(numel(s_idx),1);
    for i = 1:numel(s_idx)
        sinr_final(i) = Signal_model('sinr', params, scene, state, struct('k', s_idx(i)));
    end
end

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

if ~isempty(rates_final)
    fprintf('最终各服务用户 rate: ');
    fprintf('%.6f ', rates_final(:));
    fprintf('\n');
end
if ~isempty(sinr_final)
    fprintf('最终各服务用户 SINR: ');
    fprintf('%.6f ', sinr_final(:));
    fprintf('\n');
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
    N = size(state.X,1);
    M = size(state.X,2);

    h_users = scatter(user_pos(1,:), user_pos(2,:), 25, 'filled');

    h_waveguide = gobjects(0);
    h_pa = gobjects(0);
    for n = 1:N
        h_line = line([scene.xW(n), scene.xW(n)], [0, params.Dy]);
        if isempty(h_waveguide), h_waveguide = h_line; end

        h_sc = scatter(scene.xW(n)*ones(1,M), state.X(n,:), 40);
        if isempty(h_pa), h_pa = h_sc; end
    end

    h_serv = gobjects(0);
    if isfield(state,'S') && ~isempty(state.S)
        s_idx = state.S(:).';
        h_serv = scatter(user_pos(1,s_idx), user_pos(2,s_idx), 60);
    end

    if isfield(state,'theta') && isfield(state,'phi') && numel(state.theta)==N*M && numel(state.phi)==N*M
        x_pa = repmat(scene.xW(:), 1, M);
        x_pa = x_pa(:);
        y_pa = state.X(:);
        u = cos(state.theta(:));
        v = sin(state.theta(:));
        quiver(x_pa, y_pa, u, v, 0.4, 'Color', [0.2 0.2 0.2], 'LineWidth', 0.8, 'MaxHeadSize', 1);
    end

    xlabel('x 方向位置');
    ylabel('y 方向位置');
    title('系统几何示意图（用户、波导、PA）');

    h_list = [h_users, h_waveguide, h_pa];
    lbl_list = {'全部用户','波导','最终PA位置'};
    if ~isempty(h_serv)
        h_list = [h_list, h_serv];
        lbl_list{end+1} = '最终服务用户';
    end
    legend(h_list, lbl_list, 'Location', 'best');

    grid on;
    hold off;
end

%% ======================== 图3：初始与最终 PA 位置对比 ========================
has_init_X = false;
if isfield(history,'X0') && ~isempty(history.X0)
    X_init = history.X0;
    has_init_X = true;
elseif isfield(history,'X') && ~isempty(history.X)
    X_init = history.X{1};
    has_init_X = true;
elseif isfield(state,'X_init') && ~isempty(state.X_init)
    X_init = state.X_init;
    has_init_X = true;
end

if has_init_X && isfield(state,'X') && isfield(scene,'xW')
    figure;
    hold on;
    N = size(state.X,1);
    M = size(state.X,2);

    h_init = gobjects(0);
    h_final = gobjects(0);
    for n = 1:N
        h1 = plot(scene.xW(n)*ones(1,M), X_init(n,:), 'o');
        if isempty(h_init), h_init = h1; end
        h2 = plot(scene.xW(n)*ones(1,M), state.X(n,:), 'x');
        if isempty(h_final), h_final = h2; end
    end

    xlabel('x 方向位置');
    ylabel('y 方向位置');
    title('初始与最终 PA 位置对比');
    legend([h_init, h_final], {'初始PA位置','最终PA位置'}, 'Location', 'best');
    grid on;
    hold off;
end

%% ======================== 图3B：每轮 X 更新量 ========================
if isfield(history,'X_cells') && ~isempty(history.X_cells) && isfield(history,'X0')
    T = numel(history.X_cells);
    deltaX = zeros(T,1);
    for t = 1:T
        X_now = history.X_cells{t};
        if t == 1
            X_prev = history.X0;
        else
            X_prev = history.X_cells{t-1};
        end
        deltaX(t) = norm(X_now - X_prev, 'fro');
    end

    figure;
    plot(1:T, deltaX, '-o');
    xlabel('外层迭代轮次');
    ylabel('||X^{(t)} - X^{(t-1)}||_F');
    title('每轮 X 更新量');
    grid on;
end

%% ======================== 图3C：各 PA 位置随迭代变化轨迹 ========================
if isfield(history,'X_cells') && ~isempty(history.X_cells) && isfield(history,'X0')
    T = numel(history.X_cells);
    [N, M] = size(history.X0);

    figure;
    hold on;
    for n = 1:N
        for m = 1:M
            traj = zeros(T+1,1);
            traj(1) = history.X0(n,m);
            for t = 1:T
                traj(t+1) = history.X_cells{t}(n,m);
            end
            plot(0:T, traj, '-o');
        end
    end
    xlabel('外层迭代轮次');
    ylabel('PA y 方向位置');
    title('各 PA 位置随迭代变化轨迹');
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
    h_all = scatter(user_pos(1,:), user_pos(2,:), 20, 'filled');

    s_idx = state.S(:).';
    h_sel = scatter(user_pos(1,s_idx), user_pos(2,s_idx), 70);
    for i = 1:numel(s_idx)
        k = s_idx(i);
        text(user_pos(1,k), user_pos(2,k), sprintf(' %d', k));
    end

    xlabel('x 方向位置');
    ylabel('y 方向位置');
    title('最终服务用户选择结果');
    legend([h_all, h_sel], {'全部用户','服务用户'}, 'Location', 'best');
    grid on;
    hold off;
end

%% ======================== 图6：各块增益图 ========================
if isfield(history,'R_after_W') && isfield(history,'R_after_angle') ...
        && isfield(history,'R_after_X') && isfield(history,'R_after_S') ...
        && isfield(history,'R_sum')
    R_prev = history.R_sum(:);
    Rw = history.R_after_W(:);
    Ra = history.R_after_angle(:);
    Rx = history.R_after_X(:);
    Rs = history.R_after_S(:);

    T = min([numel(R_prev)-1, numel(Rw), numel(Ra), numel(Rx), numel(Rs)]);
    if T > 0
        R_prev = R_prev(1:T);
        Rw = Rw(1:T);
        Ra = Ra(1:T);
        Rx = Rx(1:T);
        Rs = Rs(1:T);

        deltaW = Rw - R_prev;
        deltaAngle = Ra - Rw;
        deltaX = Rx - Ra;
        deltaS = Rs - Rx;

        figure;
        bar(1:T, [deltaW, deltaAngle, deltaX, deltaS], 'stacked');
        xlabel('外层迭代轮次');
        ylabel('sum rate 增益');
        title('各块更新的 sum rate 增益');
        legend('W块','角度块','位置块','用户集块','Location','best');
        grid on;
    end
end

%% ======================== 图7：最终服务用户 rate 条形图 ========================
if ~isempty(rates_final) && isfield(state,'S') && ~isempty(state.S)
    figure;
    bar(state.S(:), rates_final(:));
    xlabel('最终服务用户编号');
    ylabel('rate (bit/s/Hz)');
    title('最终服务用户 rate');
    grid on;
end

%% ======================== 图8：最终服务用户 SINR 图 ========================
if ~isempty(sinr_final) && isfield(state,'S') && ~isempty(state.S)
    figure;
    bar(state.S(:), sinr_final(:));
    xlabel('最终服务用户编号');
    ylabel('SINR');
    title('最终服务用户 SINR');
    grid on;
end

%% ======================== 图9：用户集合演化图 ========================
if isfield(history,'S_cells') && ~isempty(history.S_cells)
    S_cells = history.S_cells;
    T = numel(S_cells);
    S_map = zeros(T, params.K);
    for t = 1:T
        s_t = S_cells{t};
        s_t = s_t(:);
        s_t = s_t(s_t>=1 & s_t<=params.K);
        S_map(t, s_t) = 1;
    end

    figure;
    imagesc(S_map);
    xlabel('用户编号');
    ylabel('外层迭代');
    title('用户集合演化图');
    cb = colorbar;
    ylabel(cb, '是否被服务(0/1)');
end

%% ======================== 图10：X 位置演化图 ========================
if isfield(history,'X_cells') && ~isempty(history.X_cells)
    X_cells = history.X_cells;
    T = numel(X_cells);

    x0 = X_cells{1};
    L = numel(x0);
    X_map = nan(T, L);

    for t = 1:T
        xt = X_cells{t};
        if ~isempty(xt) && numel(xt)==L
            X_map(t,:) = xt(:).';
        end
    end

    figure;
    imagesc(X_map);
    xlabel('PA 索引（展平后）');
    ylabel('外层迭代');
    title('X 位置演化图');
    cb = colorbar;
    ylabel(cb, 'y 坐标值');
end


%% ======================== DEBUG_X START ========================
if isfield(history,'DEBUG_X_cells') && ~isempty(history.DEBUG_X_cells)
    DEBUG_X_cells = history.DEBUG_X_cells;
    T = numel(DEBUG_X_cells);
    L = params.line_search_max_iter;

    DEBUG_X_x_proj_all = [];
    DEBUG_X_x_raw_all = [];
    DEBUG_X_delta_f_all = [];

    for t = 1:T
        DEBUG_X_x_proj_t = nan(0, L);
        DEBUG_X_x_raw_t = nan(0, L);
        DEBUG_X_delta_f_t = nan(L,1);

        DEBUG_X_t = DEBUG_X_cells{t};
        if isstruct(DEBUG_X_t) && isfield(DEBUG_X_t,'DEBUG_X_waveguides') && ~isempty(DEBUG_X_t.DEBUG_X_waveguides)
            DEBUG_X_wg1 = DEBUG_X_t.DEBUG_X_waveguides{1};
            if isstruct(DEBUG_X_wg1) && isfield(DEBUG_X_wg1,'DEBUG_X_iters') && ~isempty(DEBUG_X_wg1.DEBUG_X_iters)
                DEBUG_X_it1 = DEBUG_X_wg1.DEBUG_X_iters{1};
                if isstruct(DEBUG_X_it1)
                    if isfield(DEBUG_X_it1,'DEBUG_X_x_proj') && ~isempty(DEBUG_X_it1.DEBUG_X_x_proj)
                        DEBUG_X_x_proj_t = DEBUG_X_it1.DEBUG_X_x_proj;
                    end
                    if isfield(DEBUG_X_it1,'DEBUG_X_x_raw') && ~isempty(DEBUG_X_it1.DEBUG_X_x_raw)
                        DEBUG_X_x_raw_t = DEBUG_X_it1.DEBUG_X_x_raw;
                    end
                    if isfield(DEBUG_X_it1,'DEBUG_X_delta_f') && ~isempty(DEBUG_X_it1.DEBUG_X_delta_f)
                        DEBUG_X_delta_f_t = DEBUG_X_it1.DEBUG_X_delta_f(:);
                    end
                end
            end
        end

        if isempty(DEBUG_X_x_proj_all)
            DEBUG_X_M = size(DEBUG_X_x_proj_t,1);
            if DEBUG_X_M == 0
                DEBUG_X_M = size(DEBUG_X_x_raw_t,1);
            end
            if DEBUG_X_M == 0
                DEBUG_X_M = params.M;
            end
            DEBUG_X_x_proj_all = nan(T*L, DEBUG_X_M);
            DEBUG_X_x_raw_all = nan(T*L, DEBUG_X_M);
            DEBUG_X_delta_f_all = nan(T*L, 1);
        end

        DEBUG_X_idx = (t-1)*L + (1:L);
        if ~isempty(DEBUG_X_x_proj_t) && size(DEBUG_X_x_proj_t,2) == L
            DEBUG_X_x_proj_all(DEBUG_X_idx,1:size(DEBUG_X_x_proj_t,1)) = DEBUG_X_x_proj_t.';
        end
        if ~isempty(DEBUG_X_x_raw_t) && size(DEBUG_X_x_raw_t,2) == L
            DEBUG_X_x_raw_all(DEBUG_X_idx,1:size(DEBUG_X_x_raw_t,1)) = DEBUG_X_x_raw_t.';
        end
        if ~isempty(DEBUG_X_delta_f_t) && numel(DEBUG_X_delta_f_t) >= L
            DEBUG_X_delta_f_all(DEBUG_X_idx,1) = DEBUG_X_delta_f_t(1:L);
        elseif ~isempty(DEBUG_X_delta_f_t)
            DEBUG_X_delta_f_all(DEBUG_X_idx(1:numel(DEBUG_X_delta_f_t)),1) = DEBUG_X_delta_f_t;
        end
    end

    if ~isempty(DEBUG_X_x_proj_all)
        DEBUG_X_step_all = 1:(T*L);

        figure;
        plot(DEBUG_X_step_all, DEBUG_X_x_proj_all, '-o');
        hold on;
        for t = 1:(T-1)
            xline(t*L + 0.5, '--k');
        end
        hold off;
        xlabel('全局步编号');
        ylabel('投影后候选位置');
        title('DEBUG_X：各外层轮次中 8 次投影后候选位置（第1条波导，第1次位置内迭代）');
        grid on;

        figure;
        plot(DEBUG_X_step_all, DEBUG_X_x_raw_all, '-o');
        hold on;
        for t = 1:(T-1)
            xline(t*L + 0.5, '--k');
        end
        hold off;
        xlabel('全局步编号');
        ylabel('投影前候选位置');
        title('DEBUG_X：各外层轮次中 8 次投影前候选位置（第1条波导，第1次位置内迭代）');
        grid on;

        figure;
        stem(DEBUG_X_step_all, DEBUG_X_delta_f_all, 'filled');
        hold on;
        yline(0, '--k');
        for t = 1:(T-1)
            xline(t*L + 0.5, '--k');
        end
        hold off;
        xlabel('全局步编号');
        ylabel('候选点目标增量');
        title('DEBUG_X：各外层轮次中 8 次候选点目标增量（第1条波导，第1次位置内迭代）');
        grid on;
    end
end
%% ======================== DEBUG_X END ==========================

end
