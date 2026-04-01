function [S_new, swap_flag] = AO_S(t, S, C, Emax, X, W, theta, phi, scene, params)
% 论文用户子问题：restricted single-swap best-improvement
S_new = S; swap_flag = false;
if mod(t, params.T_S) ~= 0
    return;
end

Rk = Signal_model('individual_rates', S, X, W, theta, phi, scene, params);
[~,idWeak] = sort(Rk,'ascend');
I_weak = S(idWeak(1:min(params.L_in,numel(S))));

outs = setdiff(C,S,'stable');
[~,idOut] = sort(Emax(outs),'descend');
J_strong = outs(idOut(1:min(params.L_out,numel(outs))));

R_base = Signal_model('sum_rate', S, X, W, theta, phi, scene, params);
bestDelta = -inf; bestPair = [];
for i = 1:numel(I_weak)
    for j = 1:numel(J_strong)
        Sp = S;
        Sp(Sp==I_weak(i)) = J_strong(j);
        Rp = Signal_model('sum_rate', Sp, X, W, theta, phi, scene, params);
        Delta = Rp - R_base;
        if Delta > bestDelta
            bestDelta = Delta; bestPair = [I_weak(i), J_strong(j)];
        end
    end
end

if ~isempty(bestPair) && bestDelta >= params.eps_S
    S_new = S;
    S_new(S_new==bestPair(1)) = bestPair(2);
    swap_flag = true;
end

end
