function out = Signal_model(mode, params, scene, state, extra)
% Signal Model and Performance Metric：仅负责信号与速率计算

if nargin < 5 || isempty(extra), extra = struct(); end

switch lower(mode)
    case 'tx_signal'
        out = tx_signal_model(state, extra);

    case 'rx_signal'
        out = rx_signal_model(params, scene, state, extra);

    case 'sinr'
        out = compute_sinr_single(params, scene, state, extra.k);

    case 'user_rate'
        gamma_k = compute_sinr_single(params, scene, state, extra.k);
        out = log2(1 + gamma_k);

    case 'sum_rate'
        rates = compute_individual_rates(params, scene, state);
        out = sum(rates);

    case 'individual_rates'
        out = compute_individual_rates(params, scene, state);

    otherwise
        error('Signal_model: unsupported mode');
end

end

%% ======================== 内部子函数 ========================
function x_rad = tx_signal_model(state, extra)
% 发射信号：x_rad = sum_k w_k s_k
W = state.W;
Kserv = size(W,2);

if isfield(extra,'s')
    s = extra.s(:);
else
    s = ones(Kserv,1);
end

x_rad = W * s; % NM x 1
end

function y_k = rx_signal_model(params, scene, state, extra)
% 接收信号：y_k = desired + interference + noise
k = extra.k;
if isfield(extra,'s')
    s = extra.s(:);
else
    s = ones(numel(state.S),1);
end
if isfield(extra,'noise')
    n_k = extra.noise;
else
    n_k = 0;
end

[H, user_idx] = get_service_channel_matrix(params, scene, state);
col_k = find(user_idx == k, 1);
if isempty(col_k)
    error('Signal_model: rx_signal requires k in state.S');
end

hk = H(:, col_k);
W = state.W;

desired = (hk' * W(:,col_k)) * s(col_k);
interference = 0;
for j = 1:numel(user_idx)
    if j ~= col_k
        interference = interference + (hk' * W(:,j)) * s(j);
    end
end

y_k = desired + interference + n_k;
end

function gamma_k = compute_sinr_single(params, scene, state, k)
% 单用户SINR：gamma_k = |h_k^H w_k|^2 / (sum_{j!=k}|h_k^H w_j|^2 + sigma2)
[H, user_idx] = get_service_channel_matrix(params, scene, state);
col_k = find(user_idx == k, 1);
if isempty(col_k)
    error('Signal_model: k must be in current service set state.S');
end

hk = H(:, col_k);
W = state.W;

signal_pow = abs(hk' * W(:,col_k))^2;
interf_pow = 0;
for j = 1:numel(user_idx)
    if j ~= col_k
        interf_pow = interf_pow + abs(hk' * W(:,j))^2;
    end
end

gamma_k = signal_pow / (interf_pow + params.sigma2);
end

function rates = compute_individual_rates(params, scene, state)
% 当前服务用户集合的速率向量，顺序与 state.S 一致
[H, user_idx] = get_service_channel_matrix(params, scene, state);
W = state.W;
Kserv = numel(user_idx);

rates = zeros(Kserv,1);
for kcol = 1:Kserv
    hk = H(:,kcol);
    signal_pow = abs(hk' * W(:,kcol))^2;

    interf_pow = 0;
    for j = 1:Kserv
        if j ~= kcol
            interf_pow = interf_pow + abs(hk' * W(:,j))^2;
        end
    end

    gamma = signal_pow / (interf_pow + params.sigma2);
    rates(kcol) = log2(1 + gamma);
end
end

function [H, user_idx] = get_service_channel_matrix(params, scene, state)
% 通过 Channel_model 获取当前服务用户复合信道
extra_ch = struct();
extra_ch.use_all = false;
ch_out = Channel_model('all_users', params, scene, state, extra_ch);
H = ch_out.H;
user_idx = ch_out.user_idx(:).';
end
