function model = Problem_formulation(params, scene)
% Problem Formulation：仅定义目标、变量、约束与可行域，不做求解

if nargin < 2
    scene = struct(); %#ok<NASGU>
end

model = struct();

%% 总目标（对应论文总优化问题）
model.objective_name = 'Joint maximization of effective rate within one frame';
model.Rsum_definition = 'R_sum(S,X,W,theta,phi) = sum_{k in S} log2(1 + gamma_k)';
model.Trec_definition = 'T_rec(theta,phi,X) = max{T_X(X), T_theta(theta), T_phi(phi)}';
model.Reff_definition = 'R_eff(S,X,W,theta,phi) = [1 - T_rec(theta,phi,X)/T_f]^+ R_sum(S,X,W,theta,phi)';
model.objective = 'max_{S,X,W,theta,phi} R_eff(S,X,W,theta,phi)';

%% 优化变量
model.variables = {'S','X','W','theta','phi'};

%% 关键参数
% 服务用户数按论文定义取 K_serv = min(NRF, K_max)
% 这里在代码中实现为 min(params.NRF, params.K_max)
model.K_serv = min(params.NRF, params.K_max);
model.N = params.N;
model.M = params.M;
model.K = params.K;
model.T_f = params.T_f;
model.v_PA = params.v_PA;
model.omega_theta = params.omega_theta;
model.omega_phi = params.omega_phi;

%% 约束定义（对应论文约束式）
% 1) 用户集合约束
model.user_set_constraint = struct(...
    'set_relation', 'S subseteq C', ...
    'cardinality', '|S| = K_serv');

% 2) 位置相关约束
model.position_bounds = struct(...
    'Dy', params.Dy, ...
    'Delta', params.Delta, ...
    'ordering', '0 <= y_{n,1} <= ... <= y_{n,M} <= Dy', ...
    'spacing', 'y_{n,m} - y_{n,m-1} >= Delta, m>=2');

% 3) 角度约束
model.angle_bounds = struct(...
    'theta_min', pi/2, ...
    'theta_max', pi, ...
    'phi_range', '(-pi, pi]');

% 4) 功率约束
model.power_constraint = struct(...
    'expression', 'trace(W*W^H) <= P_max', ...
    'P_max', params.P_max);

% 5) 重构时间约束
model.reconfiguration_time_constraint = struct(...
    'expression', 'T_rec(theta,phi,X) <= T_f', ...
    'T_f', params.T_f);

%% 波导位置可行域（对应论文可行域定义）
model.position_feasible_set = ...
    'X_n = {x_n in R^M | 0 <= y_{n,1} <= ... <= y_{n,M} <= Dy, y_{n,m}-y_{n,m-1} >= Delta}';

%% 外层交替优化块顺序（供 main/AO 模块统一调用）
model.update_order = {'W','angle','X','S'};

end
