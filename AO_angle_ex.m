function [theta_new, phi_new] = AO_angle_ex(params, scene, model, state)
% AO_angle_ex：固定 (S,X,W) 更新角度 (theta,phi)（离散穷举版）

if nargin < 3
    model = struct(); %#ok<NASGU>
end

theta_new = double(state.theta);
phi_new = double(state.phi);

theta_grid = linspace(pi/2, pi, 10);
phi_grid_full = linspace(-pi, pi, 11);
phi_grid = phi_grid_full(2:end);
max_round = 20;

[N, M] = size(theta_new);

% 按PA顺序逐个更新，不并行
for r = 1:max_round
    theta_old_round = theta_new;
    phi_old_round = phi_new;

    for n = 1:N
        for m = 1:M
            [theta_nm, phi_nm] = update_single_pa_angle_ex(n, m, params, scene, state, theta_new, phi_new, theta_grid, phi_grid);
            theta_new(n,m) = theta_nm;
            phi_new(n,m) = phi_nm;
        end
    end

    if isequal(theta_new, theta_old_round) && isequal(phi_new, phi_old_round)
        break;
    end
end

end

%% ======================== 内部子函数 ========================
function [theta_best, phi_best] = update_single_pa_angle_ex(n, m, params, scene, state, theta_cur_all, phi_cur_all, theta_grid, phi_grid)
% 单个PA角度子问题更新（离散穷举100个候选组合）
theta_best = theta_cur_all(n,m);
phi_best = phi_cur_all(n,m);
R_best = -inf;

for it = 1:numel(theta_grid)
    theta_try = theta_grid(it);
    for ip = 1:numel(phi_grid)
        phi_try = phi_grid(ip);
        R_try = evaluate_angle_candidate_ex(n, m, theta_try, phi_try, params, scene, state, theta_cur_all, phi_cur_all);
        if R_try > R_best
            R_best = R_try;
            theta_best = theta_try;
            phi_best = phi_try;
        end
    end
end
end

function R = evaluate_angle_candidate_ex(n, m, theta_nm, phi_nm, params, scene, state, theta_all, phi_all)
% 仅替换单个PA角度，其他角度固定，调用真实sum rate评价
theta_candidate_all = theta_all;
phi_candidate_all = phi_all;

theta_candidate_all(n,m) = theta_nm;
phi_candidate_all(n,m) = phi_nm;

state_tmp = state;
state_tmp.theta = theta_candidate_all;
state_tmp.phi = phi_candidate_all;
R = Signal_model('sum_rate', params, scene, state_tmp, struct());
end
