function [X_new, DEBUG_X_info] = AO_X_ex(params, scene, model, state)
% AO_X_ex：固定 (S,W,theta,phi) 后，按波导逐条更新位置 X（离散穷举版）

if nargin < 3
    model = struct(); %#ok<NASGU>
end

X_new = state.X;
N = size(X_new,1);

DEBUG_X_info = struct();
DEBUG_X_info.DEBUG_X_waveguides = cell(N,1);

% 逐条波导更新，不并行
for n = 1:N
    [X_new, DEBUG_X_wg_n] = update_single_waveguide_ex(n, X_new, params, scene, state);
    DEBUG_X_info.DEBUG_X_waveguides{n} = DEBUG_X_wg_n;
end

end

%% ======================== 内部子函数 ========================
function [X_cur, DEBUG_X_wg] = update_single_waveguide_ex(n, X_cur, params, scene, state)
% 第n条波导的位置子问题更新（离散穷举 + 逐个PA交替）
G = params.Dy;
grid = linspace(0, G, 128);
max_round = 20;
M = size(X_cur,2);

DEBUG_X_wg = struct();
DEBUG_X_wg.DEBUG_X_iters = {};

for it = 1:max_round
    x_old_round = X_cur(n,:);

    % 仅记录“第1条波导、第1次位置内迭代”的逐PA轨迹
    if n == 1 && it == 1
        DEBUG_X_waveguide_states = zeros(M, M+1);
        DEBUG_X_waveguide_states(:,1) = X_cur(n,:).';
        DEBUG_X_pa_cells = cell(M,1);
    end

    for m = 1:M
        % 先显式评估当前点，作为候选集合基准
        x_now = X_cur(n,m);
        [R_cur, ~] = evaluate_position_candidate_ex(n, m, x_now, X_cur, params, scene, state);

        best_x = x_now;
        best_R = R_cur;
        best_grid_index = [];

        DEBUG_X_grid_delta_f = zeros(numel(grid),1);

        % 再遍历离散网格候选点，统一按真实 sum rate 比较
        for g = 1:numel(grid)
            x_try = grid(g);
            [R_try, X_try] = evaluate_position_candidate_ex(n, m, x_try, X_cur, params, scene, state);
            DEBUG_X_grid_delta_f(g) = R_try - R_cur;
            if R_try > best_R
                best_R = R_try;
                best_x = X_try(n,m);
                best_grid_index = g;
            end
        end

        current_grid_index = find(abs(grid - x_now) < 1e-12, 1, 'first');
        if isempty(best_grid_index)
            if ~isempty(current_grid_index)
                best_grid_index = current_grid_index;
            else
                [~, best_grid_index] = min(abs(grid - best_x));
            end
        end

        % 非下降接受：只有不劣于当前点评价时才更新
        if best_R >= R_cur
            X_cur(n,m) = best_x;
        else
            X_cur(n,m) = x_now;
        end

        % 仅记录“第1条波导、第1次位置内迭代”的穷搜调试信息
        if n == 1 && it == 1
            DEBUG_X_pa_m = struct();
            DEBUG_X_pa_m.DEBUG_X_x_cur = x_now;
            DEBUG_X_pa_m.DEBUG_X_grid_xm_proj = grid(:);
            DEBUG_X_pa_m.DEBUG_X_grid_delta_f = DEBUG_X_grid_delta_f;
            DEBUG_X_pa_m.DEBUG_X_best_grid_index = best_grid_index;
            DEBUG_X_pa_cells{m} = DEBUG_X_pa_m;

            DEBUG_X_waveguide_states(:,m+1) = X_cur(n,:).';
        end
    end

    if n == 1 && it == 1
        DEBUG_X_it = struct();
        DEBUG_X_it.DEBUG_X_pa_cells = DEBUG_X_pa_cells;
        DEBUG_X_it.DEBUG_X_waveguide_states = DEBUG_X_waveguide_states;
        DEBUG_X_wg.DEBUG_X_iters{it} = DEBUG_X_it;
    end

    if isequal(X_cur(n,:), x_old_round)
        break;
    end
end
end

function [R, X_candidate] = evaluate_position_candidate_ex(n, m, x_try, X_ref, params, scene, state)
% 固定其他变量，仅修改单个 X(n,m) 候选点并计算真实 sum rate
X_candidate = X_ref;
X_candidate(n,m) = x_try;

state_tmp = state;
state_tmp.X = X_candidate;
R = Signal_model('sum_rate', params, scene, state_tmp, struct());
end
