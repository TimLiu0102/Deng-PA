function [S_new, swap_flag] = AO_S(params, scene, model, state)
% AO_S：固定 (X,W,theta,phi) 下，按 restricted single-swap 更新用户集 S

if nargin < 3
    model = struct(); %#ok<NASGU>
end

S_new = state.S(:).';
swap_flag = false;

%% 1) 周期性触发规则：仅在 mod(t, T_S)==0 时更新
if ~isfield(state,'t')
    t = 0;
else
    t = state.t;
end
if mod(t, params.T_S) ~= 0
    return;
end

max_swaps = params.max_swaps;
if isempty(max_swaps) || max_swaps < 1
    max_swaps = 1;
end

% 可执行多次best single-swap（每次都固定当前W）
for iter_swap = 1:max_swaps
    state_now = state;
    state_now.S = S_new;

    % 2) 当前服务用户个体速率（顺序与 S_new 一致）
    rates = get_current_individual_rates(params, scene, state_now);

    % 3) 内部弱用户集合 I_weak（返回S中的位置索引和用户编号）
    [idx_weak, users_weak] = build_internal_weak_set(S_new, rates, params.L_in);

    % 4) 外部强用户集合 J_strong（按Emax粗筛）
    J_strong = build_external_strong_set(S_new, state.C, state.Emax, params.L_out);

    % 边界情况：无外部候选则不更新
    if isempty(J_strong) || isempty(users_weak)
        break;
    end

    % 当前基准R_eff
    [R_base_eff, ~] = Effective_rate_model(params, scene, state_now, []);

    % 5)-6) 枚举 single-swap 并计算交换增益 Delta_eff
    best_delta = -inf;
    best_pos = [];
    best_user_in = [];
    best_user_out = [];

    for a = 1:numel(idx_weak)
        pos_in_S = idx_weak(a);
        user_in = users_weak(a);

        for b = 1:numel(J_strong)
            user_out = J_strong(b);

            [S_candidate, delta_val] = evaluate_single_swap( ...
                params, scene, state_now, S_new, pos_in_S, user_in, user_out, R_base_eff);

            if delta_val > best_delta
                best_delta = delta_val;
                best_pos = pos_in_S;
                best_user_in = user_in;
                best_user_out = user_out;
                S_best = S_candidate;
            end
        end
    end

    % 7) best-improvement 接受准则（按R_eff增益）
    if best_delta >= params.eps_S
        S_new = S_best;
        swap_flag = true;

        % 保证不重复且长度不变
        if numel(unique(S_new)) ~= numel(S_new)
            S_new = state_now.S;
            break;
        end
    else
        break;
    end
end

end

%% ======================== 内部子函数 ========================
function rates = get_current_individual_rates(params, scene, state_now)
% 当前服务用户个体速率，调用真实模块
rates = Signal_model('individual_rates', params, scene, state_now, struct());
rates = rates(:);
end

function [idx_weak, users_weak] = build_internal_weak_set(S, rates, L_in)
% 速率从小到大选前L_in个弱用户
Kc = numel(S);
L = min(L_in, Kc);
[~, idx_sorted] = sort(rates, 'ascend');
idx_weak = idx_sorted(1:L);
users_weak = S(idx_weak);
end

function J_strong = build_external_strong_set(S, C, Emax, L_out)
% 候选池 C\S 中按 Emax 降序取前L_out个外部强用户
outside = setdiff(C(:).', S(:).', 'stable');
if isempty(outside)
    J_strong = [];
    return;
end

[~, idx] = sort(Emax(outside), 'descend');
L = min(L_out, numel(outside));
J_strong = outside(idx(1:L));
end

function [S_candidate, delta_val] = evaluate_single_swap(params, scene, state_now, S_now, pos_in_S, user_in, user_out, R_base_eff)
S_candidate = S_now;
if S_candidate(pos_in_S) ~= user_in
    % 保证位置与用户对应一致
    pos_in_S = find(S_candidate == user_in, 1);
end
S_candidate(pos_in_S) = user_out;

state_candidate = state_now;
state_candidate.S = S_candidate;

[R_candidate_eff, detail_candidate] = Effective_rate_model(params, scene, state_candidate, []);
if detail_candidate.time_feasible
    delta_val = R_candidate_eff - R_base_eff;
else
    delta_val = -inf;
end
end
