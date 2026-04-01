function Print_and_Plot(history, S, X, theta, phi, problem, scene, params)
% 仅负责结果打印与绘图
fprintf('=== MU-MISO Pinching Antenna AO复现 ===\n');
fprintf('最终sum rate: %.6f bit/s/Hz\n', history.R(end));
fprintf('外层迭代次数: %d\n', numel(history.R)-1);
fprintf('最终服务用户集合 S: ');
fprintf('%d ', S); fprintf('\n');

figure; plot(0:numel(history.R)-1, history.R, '-o','LineWidth',1.2);
xlabel('Outer Iteration'); ylabel('True Sum Rate (bit/s/Hz)'); grid on;
title('AO过程真实总和速率');

figure; hold on;
for n=1:params.N
    plot(scene.xW(n)*ones(1,params.M), X(n,:), 's', 'MarkerSize',7);
end
scatter(scene.users(:,1), scene.users(:,2), 20, 'k', 'filled');
xlabel('x'); ylabel('y'); grid on; title('波导PA位置与用户投影');
lgd = arrayfun(@(n)sprintf('WG%d',n),1:params.N,'UniformOutput',false);
legend([lgd, {'Users'}]);

disp(problem.objective);
end
