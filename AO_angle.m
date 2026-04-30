function [theta_new, phi_new] = AO_angle(params, scene, model, state)
% AO_angle：固定 (S,X,W) 更新角度 (theta,phi)

if nargin < 3
    model = struct(); %#ok<NASGU>
end

theta_new = double(state.theta);
phi_new = double(state.phi);

% 论文给定的8个polling directions
D = [ 1  0;
     -1  0;
      0  1;
      0 -1;
      1  1;
      1 -1;
     -1  1;
     -1 -1];

[N, M] = size(theta_new);

% 按PA顺序逐个更新，不并行
for n = 1:N
    for m = 1:M
        [theta_nm, phi_nm] = update_single_pa_angle(n, m, params, scene, state, theta_new, phi_new, D);
        theta_new(n,m) = theta_nm;
        phi_new(n,m) = phi_nm;
    end
end

end

%% ======================== 内部子函数 ========================
function [theta_best, phi_best] = update_single_pa_angle(n, m, params, scene, state, theta_cur_all, phi_cur_all, D)
% 单个PA的局部pattern search + re-anchoring

theta_c = theta_cur_all(n,m);
phi_c = phi_cur_all(n,m);

dtheta = params.Delta_theta0;
dphi = params.Delta_phi0;

theta_best = theta_c;
phi_best = phi_c;

re_anchor_allowed = false;
if isfield(state,'swap_flag') && state.swap_flag
    re_anchor_allowed = true;
end
no_improve_count = 0;

for r = 1:params.I_theta
    % 当前中心真实R_eff评价
    [R_center, feasible_center] = evaluate_angle_candidate(n,m,theta_c,phi_c,params,scene,state,theta_cur_all,phi_cur_all);
    if ~feasible_center
        no_improve_count = no_improve_count + 1;
        dtheta = params.beta_theta * dtheta;
        dphi   = params.beta_phi   * dphi;
        if re_anchor_allowed || no_improve_count >= 2
            [theta_c, phi_c, theta_cur_all, phi_cur_all] = do_reanchoring( ...
                n,m,params,scene,state,theta_cur_all,phi_cur_all,theta_c,phi_c);
            no_improve_count = 0;
        end
        if dtheta < 1e-4 && dphi < 1e-4
            break;
        end
        continue;
    end

    % 候选生成 + 真实R_eff评价
    R_best_round = -inf;
    theta_hat = theta_c;
    phi_hat = phi_c;
    for i = 1:size(D,1)
        theta_try = theta_c + D(i,1)*dtheta;
        phi_try = phi_c + D(i,2)*dphi;
        [theta_try, phi_try] = project_angle_pair(theta_try, phi_try);

        [R_try, feasible_try] = evaluate_angle_candidate(n,m,theta_try,phi_try,params,scene,state,theta_cur_all,phi_cur_all);
        if feasible_try && (R_try > R_best_round)
            R_best_round = R_try;
            theta_hat = theta_try;
            phi_hat = phi_try;
        end
    end

    % 非下降接受准则
    if R_best_round >= R_center + params.eps_theta
        theta_c = theta_hat;
        phi_c = phi_hat;
        theta_cur_all(n,m) = theta_c;
        phi_cur_all(n,m) = phi_c;
        no_improve_count = 0;
    else
        no_improve_count = no_improve_count + 1;
        dtheta = params.beta_theta * dtheta;
        dphi   = params.beta_phi   * dphi;

        % re-anchoring：当前角度 + 若干用户方向 + (pi,0)
        if re_anchor_allowed || no_improve_count >= 2
            [theta_c, phi_c, theta_cur_all, phi_cur_all] = do_reanchoring( ...
                n,m,params,scene,state,theta_cur_all,phi_cur_all,theta_c,phi_c);
            no_improve_count = 0;
        end
    end

    % 局部搜索终止：步长足够小
    if dtheta < 1e-4 && dphi < 1e-4
        break;
    end
end

theta_best = theta_c;
phi_best = phi_c;
end

function [theta_c, phi_c, theta_all, phi_all] = do_reanchoring(n,m,params,scene,state,theta_all,phi_all,theta_c,phi_c)
% 论文re-anchoring机制
anchors = build_anchor_set(n,m,params,scene,state,theta_all,phi_all,theta_c,phi_c);

bestR = -inf;
bestTheta = theta_c;
bestPhi = phi_c;
for i = 1:size(anchors,1)
    [th, ph] = project_angle_pair(anchors(i,1), anchors(i,2));
    [R, feasible] = evaluate_angle_candidate(n,m,th,ph,params,scene,state,theta_all,phi_all);
    if feasible && (R > bestR)
        bestR = R;
        bestTheta = th;
        bestPhi = ph;
    end
end

theta_c = bestTheta;
phi_c = bestPhi;
theta_all(n,m) = theta_c;
phi_all(n,m) = phi_c;
end

function anchors = build_anchor_set(n,m,params,scene,state,theta_all,phi_all,theta_c,phi_c)
% 锚点集合：当前角度 + 最近代表用户方向 + 默认方向(pi,0)
anchors = [theta_c, phi_c];

rep_users = select_anchor_users(n,m,params,scene,state,3);
p_nm = [scene.xW(n); state.X(n,m); params.d];
for i = 1:numel(rep_users)
    qk = scene.user_pos(:,rep_users(i));
    v = qk - p_nm;
    [th, ph] = angle_from_vector(v);
    [th, ph] = project_angle_pair(th, ph);
    anchors = [anchors; th, ph]; %#ok<AGROW>
end

anchors = [anchors; pi, 0];

% 去重（按数值近似）
anchors = unique(round(anchors,12),'rows');

% 确保仍在可行域
for i = 1:size(anchors,1)
    [anchors(i,1), anchors(i,2)] = project_angle_pair(anchors(i,1), anchors(i,2));
end
end

function users_sel = select_anchor_users(n,m,params,scene,state,Krep)
% 从当前服务用户中选与该PA几何距离最近的若干用户
if isempty(state.S)
    users_sel = [];
    return;
end

p_nm = [scene.xW(n); state.X(n,m); params.d];
S = state.S(:).';
d = zeros(1,numel(S));
for i = 1:numel(S)
    qk = scene.user_pos(:,S(i));
    d(i) = norm(qk - p_nm);
end
[~, idx] = sort(d, 'ascend');
users_sel = S(idx(1:min(Krep,numel(S))));
end

function [R, feasible] = evaluate_angle_candidate(n,m,theta_nm,phi_nm,params,scene,state,theta_all,phi_all)
% 仅替换单个PA角度，其他角度固定，调用真实R_eff评价
state_tmp = state;
state_tmp.theta = theta_all;
state_tmp.phi = phi_all;
state_tmp.theta(n,m) = theta_nm;
state_tmp.phi(n,m) = phi_nm;
[R_eff, detail] = Effective_rate_model(params, scene, state_tmp, []);
feasible = detail.time_feasible;
if feasible
    R = R_eff;
else
    R = -inf;
end
end

function [theta, phi] = project_angle_pair(theta, phi)
% 投影到角度可行域 Omega: pi/2 <= theta <= pi, -pi < phi <= pi
if theta < pi/2
    theta = pi/2;
elseif theta > pi
    theta = pi;
end
phi = atan2(sin(phi), cos(phi));
if phi <= -pi
    phi = phi + 2*pi;
end
end

function [theta, phi] = angle_from_vector(v)
% 几何方向角，和Channel_model中的坐标约定一致
r = norm(v);
if r < 1e-12
    theta = pi;
    phi = 0;
    return;
end

theta = acos(v(3)/r);
phi = atan2(v(2), v(1));
end
