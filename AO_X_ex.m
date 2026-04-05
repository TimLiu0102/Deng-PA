function X_new = AO_X_ex(params, scene, model, state)
% AO_X_ex：固定 (S,W,theta,phi) 后，按波导逐条更新位置 X（离散穷举版）

if nargin < 3
    model = struct(); %#ok<NASGU>
end

X_new = state.X;
N = size(X_new,1);

% 逐条波导更新，不并行
for n = 1:N
    X_new = update_single_waveguide_ex(n, X_new, params, scene, state);
end

end

%% ======================== 内部子函数 ========================
function X_cur = update_single_waveguide_ex(n, X_cur, params, scene, state)
% 第n条波导的位置子问题更新（离散穷举 + 逐个PA交替）
G = params.Dy;
grid = linspace(0, G, 128);
max_round = 20;

for it = 1:max_round
    x_old_round = X_cur(n,:);

    M = size(X_cur,2);
    for m = 1:M
        x_best = X_cur(n,m);
        R_best = -inf;

        for g = 1:numel(grid)
            x_try = grid(g);
            [R_try, X_try] = evaluate_position_candidate_ex(n, m, x_try, X_cur, params, scene, state);
            if R_try > R_best
                R_best = R_try;
                x_best = X_try(n,m);
            end
        end

        X_cur(n,m) = x_best;
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
