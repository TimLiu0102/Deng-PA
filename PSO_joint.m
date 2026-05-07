function [state_best, history_pso] = PSO_joint(params, scene, model, state0)
% PSO_joint：纯启发式联合优化（联合搜索 S/X/theta/phi/W）
% 说明：不调用 AO_W，直接用 Effective_rate_model 评价

if ~isfield(params, 'PSO_num_particles') || isempty(params.PSO_num_particles)
    params.PSO_num_particles = 20;
end
if ~isfield(params, 'PSO_max_iter') || isempty(params.PSO_max_iter)
    params.PSO_max_iter = 250;
end
if ~isfield(params, 'PSO_w') || isempty(params.PSO_w)
    params.PSO_w = 0.7;
end
if ~isfield(params, 'PSO_c1') || isempty(params.PSO_c1)
    params.PSO_c1 = 1.5;
end
if ~isfield(params, 'PSO_c2') || isempty(params.PSO_c2)
    params.PSO_c2 = 1.5;
end
if ~isfield(params, 'PSO_vmax_X') || isempty(params.PSO_vmax_X)
    params.PSO_vmax_X = 1.0;
end
if ~isfield(params, 'PSO_vmax_theta') || isempty(params.PSO_vmax_theta)
    params.PSO_vmax_theta = 0.15;
end
if ~isfield(params, 'PSO_vmax_phi') || isempty(params.PSO_vmax_phi)
    params.PSO_vmax_phi = 0.15;
end
if ~isfield(params, 'PSO_prob_S') || isempty(params.PSO_prob_S)
    params.PSO_prob_S = 0.4;
end
if ~isfield(params, 'PSO_beta_W') || isempty(params.PSO_beta_W)
    params.PSO_beta_W = 0.3;
end
if ~isfield(params, 'PSO_step_W') || isempty(params.PSO_step_W)
    params.PSO_step_W = 0.02;
end
if ~isfield(params, 'PSO_init_X_jitter') || isempty(params.PSO_init_X_jitter)
    params.PSO_init_X_jitter = 2.0;
end

particles = struct();
particles(1).state = state0;
particles(1).vX = zeros(params.N, params.M);
particles(1).vtheta = zeros(params.N, params.M);
particles(1).vphi = zeros(params.N, params.M);

for p = 2:params.PSO_num_particles
    state_p = state0;
    state_p.S = randperm(scene.K, params.K_serv);
    state_p.X = state0.X + params.PSO_init_X_jitter * randn(size(state0.X));
    state_p.X = project_all_waveguides(state_p.X, params);
    state_p.theta = pi/2 + (pi/2) * rand(params.N, params.M);
    state_p.phi = -pi + 2*pi * rand(params.N, params.M);

    W_try = state0.W + params.PSO_step_W * (randn(size(state0.W)) + 1j*randn(size(state0.W)));
    state_p.W = normalize_power(W_try, params.P_max);

    particles(p).state = state_p;
    particles(p).vX = zeros(params.N, params.M);
    particles(p).vtheta = zeros(params.N, params.M);
    particles(p).vphi = zeros(params.N, params.M);
end

gbest_score = -inf;
gbest_state = particles(1).state;
gbest_R_sum = -inf;
gbest_detail = struct('T_X',0,'T_theta',0,'T_phi',0,'T_rec',0,'time_factor',0);

for p = 1:params.PSO_num_particles
    [score_p, R_sum_p, detail_p] = evaluate_pso_state(params, scene, particles(p).state);
    particles(p).best_state = particles(p).state;
    particles(p).best_score = score_p;
    particles(p).best_R_sum = R_sum_p;
    particles(p).best_detail = detail_p;

    if score_p > gbest_score
        gbest_score = score_p;
        gbest_state = particles(p).state;
        gbest_R_sum = R_sum_p;
        gbest_detail = detail_p;
    end
end

history_pso = struct();
history_pso.R_eff = gbest_score;
history_pso.R_eff_current = gbest_score;
history_pso.R_eff_best = gbest_score;
history_pso.R_sum = gbest_R_sum;
history_pso.R_sum_current = gbest_R_sum;
history_pso.R_sum_best = gbest_R_sum;
history_pso.R_current = gbest_score;
history_pso.R_best = gbest_score;
history_pso.T_X = gbest_detail.T_X;
history_pso.T_theta = gbest_detail.T_theta;
history_pso.T_phi = gbest_detail.T_phi;
history_pso.T_rec = gbest_detail.T_rec;
history_pso.time_factor = gbest_detail.time_factor;
history_pso.S0 = state0.S;
history_pso.X0 = state0.X;
history_pso.theta0 = state0.theta;
history_pso.phi0 = state0.phi;
history_pso.rates0 = Signal_model('individual_rates', params, scene, state0, []);
history_pso.R_after_W = [];
history_pso.R_after_angle = [];
history_pso.R_after_X = [];
history_pso.R_after_S = [];
history_pso.X_update_mode = 'pso_joint';
history_pso.DEBUG_X_cells = {};
history_pso.swap_flag = false;
history_pso.S_cells = cell(params.PSO_max_iter, 1);
history_pso.X_cells = cell(params.PSO_max_iter, 1);
history_pso.theta_cells = cell(params.PSO_max_iter, 1);
history_pso.phi_cells = cell(params.PSO_max_iter, 1);
history_pso.W_cells = cell(params.PSO_max_iter, 1);

for iter = 1:params.PSO_max_iter
    for p = 1:params.PSO_num_particles
        state_particle = particles(p).state;
        pbest = particles(p).best_state;

        r1 = rand(params.N, params.M);
        r2 = rand(params.N, params.M);
        particles(p).vX = params.PSO_w * particles(p).vX ...
                        + params.PSO_c1 * r1 .* (pbest.X - state_particle.X) ...
                        + params.PSO_c2 * r2 .* (gbest_state.X - state_particle.X);
        particles(p).vX = max(min(particles(p).vX, params.PSO_vmax_X), -params.PSO_vmax_X);
        state_particle.X = state_particle.X + particles(p).vX;
        state_particle.X = project_all_waveguides(state_particle.X, params);

        r1 = rand(params.N, params.M);
        r2 = rand(params.N, params.M);
        particles(p).vtheta = params.PSO_w * particles(p).vtheta ...
                            + params.PSO_c1 * r1 .* (pbest.theta - state_particle.theta) ...
                            + params.PSO_c2 * r2 .* (gbest_state.theta - state_particle.theta);
        particles(p).vtheta = max(min(particles(p).vtheta, params.PSO_vmax_theta), -params.PSO_vmax_theta);
        state_particle.theta = state_particle.theta + particles(p).vtheta;
        state_particle.theta = min(max(state_particle.theta, pi/2), pi);

        r1 = rand(params.N, params.M);
        r2 = rand(params.N, params.M);
        dphi_p = angle_diff(pbest.phi, state_particle.phi);
        dphi_g = angle_diff(gbest_state.phi, state_particle.phi);
        particles(p).vphi = params.PSO_w * particles(p).vphi ...
                          + params.PSO_c1 * r1 .* dphi_p ...
                          + params.PSO_c2 * r2 .* dphi_g;
        particles(p).vphi = max(min(particles(p).vphi, params.PSO_vmax_phi), -params.PSO_vmax_phi);
        state_particle.phi = wrap_phi(state_particle.phi + particles(p).vphi);

        if rand < params.PSO_prob_S
            state_particle.S = update_user_set_from_gbest(state_particle.S, gbest_state.S, scene, params);
        end

        extra_ch = struct();
        extra_ch.use_all = false;
        ch_out = Channel_model('all_users', params, scene, state_particle, extra_ch);
        H = ch_out.H;
        W_try = state_particle.W;
        kcol = randi(size(W_try, 2));
        hk = H(:, kcol);

        if norm(hk) > 0
            hk_dir = hk / norm(hk);
            W_try(:, kcol) = (1 - params.PSO_beta_W) * W_try(:, kcol) ...
                           + params.PSO_beta_W * hk_dir ...
                           + params.PSO_step_W * (randn(size(W_try(:, kcol))) + 1j*randn(size(W_try(:, kcol))));
        else
            W_try(:, kcol) = W_try(:, kcol) ...
                           + params.PSO_step_W * (randn(size(W_try(:, kcol))) + 1j*randn(size(W_try(:, kcol))));
        end
        state_particle.W = normalize_power(W_try, params.P_max);

        [score_p, R_sum_p, detail_p] = evaluate_pso_state(params, scene, state_particle);

        particles(p).state = state_particle;

        if score_p > particles(p).best_score
            particles(p).best_state = state_particle;
            particles(p).best_score = score_p;
            particles(p).best_R_sum = R_sum_p;
            particles(p).best_detail = detail_p;
        end

        if score_p > gbest_score
            gbest_state = state_particle;
            gbest_score = score_p;
            gbest_R_sum = R_sum_p;
            gbest_detail = detail_p;
        end
    end

    history_pso.R_eff(iter+1,1) = gbest_score;
    history_pso.R_eff_current(iter+1,1) = gbest_score;
    history_pso.R_eff_best(iter+1,1) = gbest_score;
    history_pso.R_sum(iter+1,1) = gbest_R_sum;
    history_pso.R_sum_current(iter+1,1) = gbest_R_sum;
    history_pso.R_sum_best(iter+1,1) = gbest_R_sum;
    history_pso.R_current(iter+1,1) = gbest_score;
    history_pso.R_best(iter+1,1) = gbest_score;
    history_pso.T_X(iter+1,1) = gbest_detail.T_X;
    history_pso.T_theta(iter+1,1) = gbest_detail.T_theta;
    history_pso.T_phi(iter+1,1) = gbest_detail.T_phi;
    history_pso.T_rec(iter+1,1) = gbest_detail.T_rec;
    history_pso.time_factor(iter+1,1) = gbest_detail.time_factor;
    history_pso.S_cells{iter,1} = gbest_state.S;
    history_pso.X_cells{iter,1} = gbest_state.X;
    history_pso.theta_cells{iter,1} = gbest_state.theta;
    history_pso.phi_cells{iter,1} = gbest_state.phi;
    history_pso.W_cells{iter,1} = gbest_state.W;
end

state_best = gbest_state;

end

function [score, R_sum, detail] = evaluate_pso_state(params, scene, state)
[R_eff, detail] = Effective_rate_model(params, scene, state, []);
R_sum = detail.R_sum;
if detail.time_feasible
    score = R_eff;
else
    score = -inf;
end
end

function X = project_all_waveguides(X, params)
for n = 1:size(X, 1)
    row_proj = Constraint_Checker('project_position', params, X(n, :).');
    X(n, :) = row_proj(:).';
end
end

function W = normalize_power(W, P_max)
p = real(trace(W * W'));
if p > 0
    W = W * sqrt(P_max / p);
end
end

function phi = wrap_phi(phi)
phi = atan2(sin(phi), cos(phi));
phi(phi <= -pi) = phi(phi <= -pi) + 2*pi;
end

function d = angle_diff(a, b)
d = atan2(sin(a - b), cos(a - b));
end

function S_new = update_user_set_from_gbest(S, S_best, scene, params)
S_new = S(:).';
pos = randi(params.K_serv);
cand_from_best = setdiff(S_best(:).', S_new);
if ~isempty(cand_from_best)
    S_new(pos) = cand_from_best(randi(numel(cand_from_best)));
else
    cand_all = setdiff(1:scene.K, S_new);
    if ~isempty(cand_all)
        S_new(pos) = cand_all(randi(numel(cand_all)));
    end
end

S_new = unique(S_new, 'stable');
if numel(S_new) < params.K_serv
    cand_all = setdiff(1:scene.K, S_new);
    need = params.K_serv - numel(S_new);
    if need > 0 && ~isempty(cand_all)
        add_num = min(need, numel(cand_all));
        add_idx = randperm(numel(cand_all), add_num);
        S_new = [S_new, cand_all(add_idx)]; %#ok<AGROW>
    end
end
S_new = S_new(1:params.K_serv);
end
