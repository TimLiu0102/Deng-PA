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
    [X_new, DEBUG_X_waveguide] = update_single_waveguide_ex(n, X_new, params, scene, state);
    DEBUG_X_info.DEBUG_X_waveguides{n,1} = DEBUG_X_waveguide;
end

end

%% ======================== 内部子函数 ========================
function [X_cur, DEBUG_X_waveguide] = update_single_waveguide_ex(n, X_cur, params, scene, state)
% 第n条波导的位置子问题更新（离散穷举 + 逐个PA交替）
G = params.Dy;
grid = linspace(0, G, 128);
max_round = 20;
L = params.line_search_max_iter;

DEBUG_X_waveguide = struct();
DEBUG_X_waveguide.DEBUG_X_iters = cell(params.I_X,1);
DEBUG_X_waveguide.DEBUG_X_g_vec = cell(params.I_X,1);
DEBUG_X_waveguide.DEBUG_X_d_vec = cell(params.I_X,1);
DEBUG_X_waveguide.DEBUG_X_R_old = nan(params.I_X,1);
DEBUG_X_waveguide.DEBUG_X_R_trial = nan(params.I_X,1);

M = size(X_cur,2);

for it = 1:max_round
    x_old_round = X_cur(n,:);

    state_round = state;
    state_round.X = X_cur;
    R_round_old = Signal_model('sum_rate', params, scene, state_round, struct());

    DEBUG_X_iter = struct();
    DEBUG_X_iter.DEBUG_X_x_raw = nan(M, L);
    DEBUG_X_iter.DEBUG_X_x_proj = nan(M, L);
    DEBUG_X_iter.DEBUG_X_delta_f = nan(L, 1);
    debug_col = 0;

    for m = 1:M
        x_cur_m = X_cur(n,m);
        [R_cur, ~] = evaluate_position_candidate_ex(n, m, x_cur_m, X_cur, params, scene, state);

        % 当前点作为候选基准，若离散网格中不存在更优点，则保持当前位置不变
        x_best = x_cur_m;
        R_best = R_cur;

        for g = 1:numel(grid)
            x_try = grid(g);

            X_raw = X_cur;
            X_raw(n,m) = x_try;
            x_raw_n = X_raw(n,:).';

            [R_try, X_try] = evaluate_position_candidate_ex(n, m, x_try, X_cur, params, scene, state);

            if debug_col < L
                debug_col = debug_col + 1;
                DEBUG_X_iter.DEBUG_X_x_raw(:,debug_col) = x_raw_n;
                DEBUG_X_iter.DEBUG_X_x_proj(:,debug_col) = X_try(n,:).';
                DEBUG_X_iter.DEBUG_X_delta_f(debug_col,1) = R_try - R_cur;
            end

            if R_try > R_best
                R_best = R_try;
                x_best = X_try(n,m);
            end
        end

        X_cur(n,m) = x_best;
    end

    state_round_new = state;
    state_round_new.X = X_cur;
    R_round_new = Signal_model('sum_rate', params, scene, state_round_new, struct());

    if it <= params.I_X
        DEBUG_X_waveguide.DEBUG_X_iters{it,1} = DEBUG_X_iter;
        DEBUG_X_waveguide.DEBUG_X_g_vec{it,1} = [];
        DEBUG_X_waveguide.DEBUG_X_d_vec{it,1} = [];
        DEBUG_X_waveguide.DEBUG_X_R_old(it,1) = R_round_old;
        DEBUG_X_waveguide.DEBUG_X_R_trial(it,1) = R_round_new;
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

xbar_n = X_candidate(n,:).';
xproj_n = Constraint_Checker('project_position', params, xbar_n);
X_candidate(n,:) = xproj_n(:).';

state_tmp = state;
state_tmp.X = X_candidate;
R = Signal_model('sum_rate', params, scene, state_tmp, struct());
end
