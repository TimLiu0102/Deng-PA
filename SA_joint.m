function [state_best, history_sa] = SA_joint(params, scene, model, state0)
% SA_joint：纯启发式联合优化（联合扰动 S/X/theta/phi/W）
% 说明：不再调用 AO_W，W 也作为搜索变量一起做随机扰动

% 兼容：若主脚本未显式提供 SA 参数，则使用默认值
if ~isfield(params, 'SA_max_iter') || isempty(params.SA_max_iter)
    params.SA_max_iter = 500;
end
if ~isfield(params, 'SA_T0') || isempty(params.SA_T0)
    params.SA_T0 = 1.0;
end
if ~isfield(params, 'SA_alpha') || isempty(params.SA_alpha)
    params.SA_alpha = 0.985;
end
if ~isfield(params, 'SA_step_X') || isempty(params.SA_step_X)
    params.SA_step_X = 0.5;
end
if ~isfield(params, 'SA_step_theta') || isempty(params.SA_step_theta)
    params.SA_step_theta = 0.08;
end
if ~isfield(params, 'SA_step_phi') || isempty(params.SA_step_phi)
    params.SA_step_phi = 0.08;
end
if ~isfield(params, 'SA_step_W') || isempty(params.SA_step_W)
    params.SA_step_W = 0.05;
end

% 纯启发式：初始点直接使用 state0，不做 AO_W 精修
state_cur = state0;
R_cur = Signal_model('sum_rate', params, scene, state_cur, []);

state_best = state_cur;
R_best = R_cur;

history_sa = struct();
history_sa.R_sum = R_best;
history_sa.R_current = R_cur;
history_sa.R_best = R_best;
history_sa.S0 = state0.S;
history_sa.X0 = state0.X;
history_sa.theta0 = state0.theta;
history_sa.phi0 = state0.phi;
history_sa.rates0 = Signal_model('individual_rates', params, scene, state_cur, []);
history_sa.R_after_W = [];
history_sa.R_after_angle = [];
history_sa.R_after_X = [];
history_sa.R_after_S = [];
history_sa.swap_flag = false;

history_sa.accept_flag = false(params.SA_max_iter, 1);
history_sa.move_type = cell(params.SA_max_iter, 1);
history_sa.S_cells = cell(params.SA_max_iter, 1);
history_sa.X_cells = cell(params.SA_max_iter, 1);
history_sa.theta_cells = cell(params.SA_max_iter, 1);
history_sa.phi_cells = cell(params.SA_max_iter, 1);
history_sa.W_cells = cell(params.SA_max_iter, 1);

for iter = 1:params.SA_max_iter
    T = max(params.SA_T0 * params.SA_alpha^(iter - 1), eps);
    state_try = state_cur;

    % 五种扰动：S / X / theta / phi / W
    r_move = randi(5);
    if r_move == 1
        move_type = 'S';

        S_row = state_try.S(:).';
        if isfield(state_try, 'C') && ~isempty(state_try.C)
            pool = state_try.C(:).';
        else
            pool = 1:scene.K;
        end

        pos = randi(numel(S_row));
        cand = setdiff(pool, S_row);
        if ~isempty(cand)
            user_new = cand(randi(numel(cand)));
            S_row(pos) = user_new;
        end
        state_try.S = S_row(:).';

    elseif r_move == 2
        move_type = 'X';

        [N, M] = size(state_try.X);
        n = randi(N);
        m = randi(M);
        state_try.X(n, m) = state_try.X(n, m) + params.SA_step_X * randn;

        row_proj = Constraint_Checker('project_position', params, state_try.X(n, :).');
        state_try.X(n, :) = row_proj(:).';

    elseif r_move == 3
        move_type = 'theta';

        [N, M] = size(state_try.theta);
        n = randi(N);
        m = randi(M);
        state_try.theta(n, m) = state_try.theta(n, m) + params.SA_step_theta * randn;
        state_try.theta(n, m) = min(max(state_try.theta(n, m), pi/2), pi);

    elseif r_move == 4
        move_type = 'phi';

        [N, M] = size(state_try.phi);
        n = randi(N);
        m = randi(M);
        state_try.phi(n, m) = state_try.phi(n, m) + params.SA_step_phi * randn;
        state_try.phi(n, m) = mod(state_try.phi(n, m) + pi, 2*pi) - pi;
        if state_try.phi(n, m) <= -pi
            state_try.phi(n, m) = pi;
        end

    else
        move_type = 'W';

        W_try = state_try.W;
        W_try = W_try + params.SA_step_W * (randn(size(W_try)) + 1j*randn(size(W_try)));

        p_try = real(trace(W_try * W_try'));
        if p_try > 0
            W_try = W_try * sqrt(params.P_max / p_try);
        end

        state_try.W = W_try;
    end

    % 纯启发式：试探点评价不再调用 AO_W
    R_try = Signal_model('sum_rate', params, scene, state_try, []);

    delta = R_try - R_cur;
    accepted = false;
    if delta >= 0
        accepted = true;
    else
        if rand < exp(delta / T)
            accepted = true;
        end
    end

    if accepted
        state_cur = state_try;
        R_cur = R_try;
    end

    if R_try > R_best
        state_best = state_try;
        R_best = R_try;
    end

    history_sa.R_sum(iter+1,1) = R_best;
    history_sa.R_current(iter+1,1) = R_cur;
    history_sa.R_best(iter+1,1) = R_best;
    history_sa.accept_flag(iter,1) = accepted;
    history_sa.move_type{iter,1} = move_type;
    history_sa.S_cells{iter,1} = state_cur.S;
    history_sa.X_cells{iter,1} = state_cur.X;
    history_sa.theta_cells{iter,1} = state_cur.theta;
    history_sa.phi_cells{iter,1} = state_cur.phi;
    history_sa.W_cells{iter,1} = state_cur.W;
end

end
