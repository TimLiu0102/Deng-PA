function [R_eff, detail] = Effective_rate_model(params, scene, state, extra)
% Effective Rate Model：统一计算有效速率

if nargin < 4 || isempty(extra), extra = struct(); end

R_sum = Signal_model('sum_rate', params, scene, state, []);

if ~isfield(params,'T_f'), params.T_f = 1; end
if ~isfield(params,'v_PA'), params.v_PA = 1; end
if ~isfield(params,'omega_theta'), params.omega_theta = 1; end
if ~isfield(params,'omega_phi'), params.omega_phi = 1; end

[X_ref, theta_ref, phi_ref] = build_reference_configuration(params, extra);

T_X = max(abs(state.X(:) - X_ref(:))) / params.v_PA;
T_theta = max(abs(state.theta(:) - theta_ref(:))) / params.omega_theta;
T_phi = max(angle_distance_pi(state.phi(:), phi_ref(:))) / params.omega_phi;

T_rec = max([T_X, T_theta, T_phi]);

time_factor = max(1 - T_rec / params.T_f, 0);
R_eff = time_factor * R_sum;

detail = struct();
detail.R_sum = R_sum;
detail.R_eff = R_eff;
detail.T_X = T_X;
detail.T_theta = T_theta;
detail.T_phi = T_phi;
detail.T_rec = T_rec;
detail.time_factor = time_factor;
detail.time_feasible = (T_rec <= params.T_f + 1e-12);

end

%% ======================== 内部子函数 ========================
function [X_ref, theta_ref, phi_ref] = build_reference_configuration(params, extra)

if isfield(extra,'X_ref') && ~isempty(extra.X_ref)
    X_ref = extra.X_ref;
else
    y_ref = zeros(params.N, params.M);
    if params.M == 1
        y_ref(:,1) = 0;
    else
        for m = 1:params.M
            y_ref(:,m) = (m-1) * params.Delta + ((m-1)/(params.M-1)) * (params.Dy - (params.M-1)*params.Delta);
        end
    end
    X_ref = y_ref;
end

if isfield(extra,'theta_ref') && ~isempty(extra.theta_ref)
    theta_ref = extra.theta_ref;
else
    theta_ref = pi * ones(params.N, params.M);
end

if isfield(extra,'phi_ref') && ~isempty(extra.phi_ref)
    phi_ref = extra.phi_ref;
else
    phi_ref = zeros(params.N, params.M);
end

end

function d = angle_distance_pi(a, b)

d = abs(mod(a - b + pi, 2*pi) - pi);

end
