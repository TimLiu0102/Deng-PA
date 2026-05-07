function [S_new, W_new, swap_flag] = AO_S_fixed_reW(params, scene, model, state)
% AO_S_fixed_reW：固定天线构型下，按 restricted single-swap 更新用户集 S；
% 每个候选 swap 都重新调用 AO_W，并用重优化后的 R_eff 评估增益。

if nargin < 3
    model = struct(); %#ok<NASGU>
end

S_new = state.S(:).';
W_new = state.W;
swap_flag = false;

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

for iter_swap = 1:max_swaps
    state_now = state;
    state_now.S = S_new;
    state_now.W = W_new;

    state_base = state_now;
    state_base.W = AO_W(params, scene, model, state_base);
    [R_base_eff, ~] = Effective_rate_model(params, scene, state_base, []);
    W_base = state_base.W;

    rates = Signal_model('individual_rates', params, scene, state_base, struct());
    rates = rates(:);

    [idx_weak, users_weak] = build_internal_weak_set(S_new, rates, params.L_in);
    J_strong = build_external_strong_set_fixed(params, scene, state_base, S_new);

    if isempty(J_strong) || isempty(users_weak)
        W_new = W_base;
        break;
    end

    best_delta = -inf;
    S_best = S_new;
    W_best = W_base;

    for a = 1:numel(idx_weak)
        pos_in_S = idx_weak(a);
        user_in = users_weak(a);

        for b = 1:numel(J_strong)
            user_out = J_strong(b);

            S_candidate = S_new;
            pos_tmp = pos_in_S;
            if S_candidate(pos_tmp) ~= user_in
                pos_tmp = find(S_candidate == user_in, 1);
            end
            S_candidate(pos_tmp) = user_out;

            state_candidate = state_base;
            state_candidate.S = S_candidate;
            state_candidate.W = [];
            state_candidate.W = AO_W(params, scene, model, state_candidate);
            [R_candidate_eff, detail_candidate] = Effective_rate_model(params, scene, state_candidate, []);

            if detail_candidate.time_feasible
                delta_val = R_candidate_eff - R_base_eff;
            else
                delta_val = -inf;
            end

            if delta_val > best_delta
                best_delta = delta_val;
                S_best = S_candidate;
                W_best = state_candidate.W;
            end
        end
    end

    if best_delta >= params.eps_S
        S_new = S_best;
        W_new = W_best;
        swap_flag = true;

        if numel(unique(S_new)) ~= numel(S_new)
            S_new = state_now.S;
            W_new = W_base;
            break;
        end
    else
        W_new = W_base;
        break;
    end
end

end

function [idx_weak, users_weak] = build_internal_weak_set(S, rates, L_in)
Kc = numel(S);
L = min(L_in, Kc);
[~, idx_sorted] = sort(rates, 'ascend');
idx_weak = idx_sorted(1:L);
users_weak = S(idx_weak);
end

function J_strong = build_external_strong_set_fixed(params, scene, state_now, S_now)
extra_ch = struct();
extra_ch.use_all = true;
ch_out = Channel_model('all_users', params, scene, state_now, extra_ch);
H_all = ch_out.H;

fixed_gain = zeros(1, scene.K);
for k = 1:scene.K
    fixed_gain(k) = norm(H_all(:,k));
end

outside = setdiff(1:scene.K, S_now, 'stable');
if isempty(outside)
    J_strong = [];
    return;
end

[~, idx] = sort(fixed_gain(outside), 'descend');
L = min(params.L_out, numel(outside));
J_strong = outside(idx(1:L));
end
