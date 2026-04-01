function out = Signal_model(mode, S, X, W, theta, phi, scene, params)
% 仅负责信号模型/SINR/速率（论文信号模型部分）
params.scene = scene;
H = Channel_model('composite_channels', S, X, theta, phi, [], params);

if isempty(W)
    W = simple_w_init(H, params.P_max);
end

switch mode
    case 'sum_rate'
        out = compute_sum_rate(H, W, params.sigma2);
    case 'individual_rates'
        out = compute_individual_rates(H, W, params.sigma2);
    case 'sinr'
        out = compute_sinr(H, W, params.sigma2);
    otherwise
        error('Signal_model: unknown mode');
end

end

function gamma = compute_sinr(H, W, sigma2)
K = size(H,2); gamma = zeros(K,1);
for k = 1:K
    hk = H(:,k);
    sig = abs(hk' * W(:,k))^2;
    interf = 0;
    for j = 1:K
        if j ~= k
            interf = interf + abs(hk' * W(:,j))^2;
        end
    end
    gamma(k) = sig / (interf + sigma2);
end
end

function rk = compute_individual_rates(H, W, sigma2)
gamma = compute_sinr(H, W, sigma2);
rk = log2(1 + gamma);
end

function R = compute_sum_rate(H, W, sigma2)
rk = compute_individual_rates(H, W, sigma2);
R = sum(rk);
end

function W = simple_w_init(H, Pmax)
K = size(H,2);
W = H;
if norm(W,'fro') < 1e-12
    W = randn(size(H))+1j*randn(size(H));
end
W = W / norm(W,'fro') * sqrt(Pmax);
if K == 0
    W = zeros(size(H,1),0);
end
end
