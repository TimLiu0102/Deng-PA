function [state, history] = AO_fixedX(params, scene, model, state)
% AO_fixedX：固定X的交替优化（W -> angle -> fixed X -> S）

rates0 = Signal_model('individual_rates', params, scene, state, []);
R_sum0 = sum(rates0);
[R_eff0, detail0] = Effective_rate_model(params, scene, state, []);
R_old = R_eff0;

history = struct();
history.X0 = state.X;
history.theta0 = state.theta;
history.phi0 = state.phi;
history.S0 = state.S;
history.rates0 = rates0;

history.R_sum = R_sum0;
history.R_eff = R_eff0;
history.T_X = detail0.T_X;
history.T_theta = detail0.T_theta;
history.T_phi = detail0.T_phi;
history.T_rec = detail0.T_rec;
history.time_factor = detail0.time_factor;

history.R_after_W = [];
history.R_after_angle = [];
history.R_after_X = [];
history.R_after_S = [];
history.R_eff_after_W = [];
history.R_eff_after_angle = [];
history.R_eff_after_X = [];
history.R_eff_after_S = [];

history.S_cells = {};
history.X_cells = {};
history.theta_cells = {};
history.phi_cells = {};
history.swap_flag = false;
history.X_update_mode = 'fixedX';

for t = 1:params.T_max
    state.t = t;

    state.W = AO_W(params, scene, model, state);
    R_after_W = Signal_model('sum_rate', params, scene, state, []);
    [R_eff_after_W, ~] = Effective_rate_model(params, scene, state, []);

    [state.theta, state.phi] = AO_angle(params, scene, model, state);
    R_after_angle = Signal_model('sum_rate', params, scene, state, []);
    [R_eff_after_angle, ~] = Effective_rate_model(params, scene, state, []);

    R_after_X = R_after_angle;
    R_eff_after_X = R_eff_after_angle;

    [state.S, state.swap_flag] = AO_S(params, scene, model, state);
    R_after_S = Signal_model('sum_rate', params, scene, state, []);
    [R_eff_after_S, detail_S] = Effective_rate_model(params, scene, state, []);

    history.R_after_W(end+1,1) = R_after_W;
    history.R_after_angle(end+1,1) = R_after_angle;
    history.R_after_X(end+1,1) = R_after_X;
    history.R_after_S(end+1,1) = R_after_S;
    history.R_eff_after_W(end+1,1) = R_eff_after_W;
    history.R_eff_after_angle(end+1,1) = R_eff_after_angle;
    history.R_eff_after_X(end+1,1) = R_eff_after_X;
    history.R_eff_after_S(end+1,1) = R_eff_after_S;

    R_new = R_eff_after_S;
    history.R_eff(end+1,1) = R_new;
    history.R_sum(end+1,1) = R_after_S;
    history.T_X(end+1,1) = detail_S.T_X;
    history.T_theta(end+1,1) = detail_S.T_theta;
    history.T_phi(end+1,1) = detail_S.T_phi;
    history.T_rec(end+1,1) = detail_S.T_rec;
    history.time_factor(end+1,1) = detail_S.time_factor;

    history.S_cells{t,1} = state.S;
    history.X_cells{t,1} = state.X;
    history.theta_cells{t,1} = state.theta;
    history.phi_cells{t,1} = state.phi;
    history.swap_flag(end+1,1) = state.swap_flag;

    if abs(R_new - R_old) < params.eps_outer
        break;
    end

    R_old = R_new;
end
end
