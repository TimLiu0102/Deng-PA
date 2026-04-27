function [X_new, DEBUG_X_info] = AO_X_grid_refine(params, scene, model, state)
% AO_X_grid_refine：固定 (S,W,theta,phi) 后，按波导逐条更新位置 X（粗到细网格细化版）

if nargin < 3
    model = struct(); %#ok<NASGU>
end

% 默认参数（保持简洁，缺省时直接使用）
if ~isfield(params, 'I_X_refine'), params.I_X_refine = 3; end
if ~isfield(params, 'X_refine_grid_num'), params.X_refine_grid_num = 41; end
if ~isfield(params, 'X_refine_radius_list'), params.X_refine_radius_list = [inf, 0.5, 0.1, 0.02, 0.005]; end
if ~isfield(params, 'eps_X_refine'), params.eps_X_refine = 0; end
if ~isfield(params, 'X_refine_use_ystar'), params.X_refine_use_ystar = true; end

X_new = state.X;
N = size(X_new,1);

DEBUG_X_info = struct();
DEBUG_X_info.DEBUG_X_waveguides = cell(N,1);

% 逐条波导更新，不并行
for n = 1:N
    [X_new, DEBUG_X_wg_n] = update_single_waveguide_grid_refine(n, X_new, params, scene, state);
    DEBUG_X_info.DEBUG_X_waveguides{n} = DEBUG_X_wg_n;
end

end

%% ======================== 内部子函数 ========================
function [X_cur, DEBUG_X_wg] = update_single_waveguide_grid_refine(n, X_cur, params, scene, state)
% 第n条波导的位置子问题更新（粗到细 + 逐个PA交替）
M = size(X_cur,2);

DEBUG_X_wg = struct();
DEBUG_X_wg.DEBUG_X_iters = {};

for it = 1:params.I_X_refine
    x_old_round = X_cur(n,:);

    % 细化半径：超过给定列表长度时，重复使用最后一个
    radius_idx = min(it, numel(params.X_refine_radius_list));
    radius = params.X_refine_radius_list(radius_idx);

    % 仅对第1条波导第1轮保留完整调试轨迹，避免调试信息过大
    need_detail = (n == 1 && it == 1);
    if need_detail
        DEBUG_X_waveguide_states = zeros(M, M+1);
        DEBUG_X_waveguide_states(:,1) = X_cur(n,:).';
        DEBUG_X_pa_cells = cell(M,1);
    end

    for m = 1:M
        % 由相邻PA约束得到该PA的可行区间
        if m == 1
            lb = 0;
        else
            lb = X_cur(n,m-1) + params.Delta;
        end

        if m == M
            ub = params.Dy;
        else
            ub = X_cur(n,m+1) - params.Delta;
        end

        x_now = X_cur(n,m);

        % 当前点评价
        state_cur = state;
        state_cur.X = X_cur;
        R_cur = Signal_model('sum_rate', params, scene, state_cur, []);

        % 无可行空间：跳过
        if lb > ub
            if need_detail
                DEBUG_X_pa_m = struct();
                DEBUG_X_pa_m.DEBUG_X_x_cur = x_now;
                DEBUG_X_pa_m.DEBUG_X_grid_xm_proj = x_now;
                DEBUG_X_pa_m.DEBUG_X_grid_delta_f = 0;
                DEBUG_X_pa_m.DEBUG_X_best_grid_index = 1;
                DEBUG_X_pa_m.DEBUG_X_best_x = x_now;
                DEBUG_X_pa_m.DEBUG_X_best_R = R_cur;
                DEBUG_X_pa_m.DEBUG_X_R_cur = R_cur;
                DEBUG_X_pa_m.DEBUG_X_accept = false;
                DEBUG_X_pa_m.DEBUG_X_lb = lb;
                DEBUG_X_pa_m.DEBUG_X_ub = ub;
                DEBUG_X_pa_m.DEBUG_X_local_lb = x_now;
                DEBUG_X_pa_m.DEBUG_X_local_ub = x_now;
                DEBUG_X_pa_cells{m} = DEBUG_X_pa_m;
                DEBUG_X_waveguide_states(:,m+1) = X_cur(n,:).';
            end
            continue;
        end

        % 粗到细局部搜索区间
        if isinf(radius)
            local_lb = lb;
            local_ub = ub;
        else
            local_lb = max(lb, x_now - radius);
            local_ub = min(ub, x_now + radius);
        end

        % 候选点：基础网格 + 当前点 + y_star引导点
        if local_lb == local_ub
            grid = local_lb;
        else
            grid = linspace(local_lb, local_ub, params.X_refine_grid_num);
        end

        candidates = grid(:);
        candidates(end+1,1) = x_now;

        if params.X_refine_use_ystar && isfield(state, 'y_star') && ~isempty(state.y_star)
            if isfield(state, 'S') && ~isempty(state.S)
                for ii = 1:numel(state.S)
                    k = state.S(ii);
                    if k >= 1 && k <= size(state.y_star,1)
                        candidates(end+1,1) = state.y_star(k,n,m);
                    end
                end
            end
        end

        candidates = min(max(candidates, lb), ub);
        candidates = unique(candidates);

        % 候选点评价
        R_candidates = zeros(numel(candidates),1);
        best_R = -inf;
        best_x = x_now;
        best_idx = 1;

        for c = 1:numel(candidates)
            x_try = candidates(c);
            X_candidate = X_cur;
            X_candidate(n,m) = x_try;

            state_tmp = state;
            state_tmp.X = X_candidate;
            R_try = Signal_model('sum_rate', params, scene, state_tmp, []);

            R_candidates(c) = R_try;
            if R_try > best_R
                best_R = R_try;
                best_x = x_try;
                best_idx = c;
            end
        end

        % 非下降接受
        accept_flag = (best_R >= R_cur + params.eps_X_refine);
        if accept_flag
            X_cur(n,m) = best_x;
        else
            X_cur(n,m) = x_now;
        end

        if need_detail
            DEBUG_X_pa_m = struct();
            DEBUG_X_pa_m.DEBUG_X_x_cur = x_now;
            DEBUG_X_pa_m.DEBUG_X_grid_xm_proj = candidates(:);
            DEBUG_X_pa_m.DEBUG_X_grid_delta_f = R_candidates(:) - R_cur;
            DEBUG_X_pa_m.DEBUG_X_best_grid_index = best_idx;
            DEBUG_X_pa_m.DEBUG_X_best_x = best_x;
            DEBUG_X_pa_m.DEBUG_X_best_R = best_R;
            DEBUG_X_pa_m.DEBUG_X_R_cur = R_cur;
            DEBUG_X_pa_m.DEBUG_X_accept = accept_flag;
            DEBUG_X_pa_m.DEBUG_X_lb = lb;
            DEBUG_X_pa_m.DEBUG_X_ub = ub;
            DEBUG_X_pa_m.DEBUG_X_local_lb = local_lb;
            DEBUG_X_pa_m.DEBUG_X_local_ub = local_ub;
            DEBUG_X_pa_cells{m} = DEBUG_X_pa_m;

            DEBUG_X_waveguide_states(:,m+1) = X_cur(n,:).';
        end
    end

    % 每轮末尾做一次投影（保险）
    x_proj = Constraint_Checker('project_position', params, X_cur(n,:).');
    X_cur(n,:) = x_proj(:).';

    DEBUG_X_it = struct();
    DEBUG_X_it.DEBUG_X_radius = radius;
    DEBUG_X_it.DEBUG_X_x_old_round = x_old_round;
    DEBUG_X_it.DEBUG_X_x_new_round = X_cur(n,:);
    DEBUG_X_it.DEBUG_X_round_changed = ~isequal(X_cur(n,:), x_old_round);
    if need_detail
        DEBUG_X_it.DEBUG_X_pa_cells = DEBUG_X_pa_cells;
        DEBUG_X_it.DEBUG_X_waveguide_states = DEBUG_X_waveguide_states;
    else
        DEBUG_X_it.DEBUG_X_pa_cells = {};
        DEBUG_X_it.DEBUG_X_waveguide_states = [];
    end
    DEBUG_X_wg.DEBUG_X_iters{it} = DEBUG_X_it;

    if isequal(X_cur(n,:), x_old_round)
        break;
    end
end
end
