function X_new = AO_X_ex(params, scene, model, state)
% AO_X_ex：固定 (S,W,theta,phi) 后，按波导逐条联合枚举更新位置 X

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
% 第n条波导的位置子问题更新（联合组合枚举）
X_candidates = generate_waveguide_position_combinations_ex(params, size(X_cur,2));

% 把当前位置也作为候选
x_cur = X_cur(n,:);
X_candidates = [X_candidates; x_cur];

R_best = -inf;
x_best = x_cur;
for i = 1:size(X_candidates,1)
    x_try = X_candidates(i,:);
    R_try = evaluate_position_candidate_ex(n, x_try, X_cur, params, scene, state);
    if R_try > R_best
        R_best = R_try;
        x_best = x_try;
    end
end

X_cur(n,:) = x_best;
end

function X_candidates = generate_waveguide_position_combinations_ex(params, M)
% 生成单条波导位置组合候选：从网格中选M个点
G = params.Dy;
grid = linspace(0, G, 128);

idx_comb = nchoosek(1:numel(grid), M);
X_candidates = grid(idx_comb);
X_candidates = filter_feasible_position_combinations_ex(X_candidates, params.Delta);
end

function X_ok = filter_feasible_position_combinations_ex(X_in, Delta)
% 只保留满足最小间距约束的组合
if isempty(X_in)
    X_ok = X_in;
    return;
end

diff_x = diff(X_in, 1, 2);
mask = all(diff_x >= Delta, 2);
X_ok = X_in(mask,:);
end

function R = evaluate_position_candidate_ex(n, x_try, X_ref, params, scene, state)
% 固定其他变量，仅替换第n条波导位置并计算真实 sum rate
X_candidate = X_ref;
X_candidate(n,:) = x_try;

state_tmp = state;
state_tmp.X = X_candidate;
R = Signal_model('sum_rate', params, scene, state_tmp, struct());
end
