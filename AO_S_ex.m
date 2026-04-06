function [S_new, swap_flag] = AO_S_ex(params, scene, model, state)
% AO_S_ex：固定 (X,W,theta,phi) 下，按用户组合与排列穷举更新用户集 S

if nargin < 3
    model = struct(); %#ok<NASGU>
end

S_cur = state.S(:).';
K_serv = params.K_serv;

pool = get_user_pool_ex(scene, state);
pool = pool(:).';

all_sets = nchoosek(pool, K_serv);

best_R = -inf;
S_best = S_cur;
for i = 1:size(all_sets,1)
    S_set = all_sets(i,:);
    [S_set_best, R_set_best] = evaluate_user_set_permutations_ex(params, scene, state, S_set);
    if R_set_best > best_R
        best_R = R_set_best;
        S_best = S_set_best;
    end
end

S_new = S_best(:).';
swap_flag = ~isequal(S_new, S_cur);

end

%% ======================== 内部子函数 ========================
function pool = get_user_pool_ex(scene, state)
% 用户搜索池：优先用 state.C，否则用全部用户
if isfield(state,'C') && ~isempty(state.C)
    pool = state.C;
else
    pool = 1:scene.K;
end
end

function [S_best, R_best] = evaluate_user_set_permutations_ex(params, scene, state, S_set)
% 对一个组合的全部排列穷举，返回评分最高的排列
all_perm = perms(S_set);

R_best = -inf;
S_best = S_set(:).';
for p = 1:size(all_perm,1)
    S_candidate = all_perm(p,:);
    R_candidate = evaluate_user_set_ex(params, scene, state, S_candidate);
    if R_candidate > R_best
        R_best = R_candidate;
        S_best = S_candidate;
    end
end
end

function R = evaluate_user_set_ex(params, scene, state, S_candidate)
% 固定当前 (X,W,theta,phi)，仅替换 S 并计算真实 sum rate
state_tmp = state;
state_tmp.S = S_candidate(:).';
R = Signal_model('sum_rate', params, scene, state_tmp, struct());
end
