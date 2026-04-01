function problem = Problem_formulation(params)
% 仅定义论文总优化问题与可行域，不做求解

problem.objective = 'max_{S,X,W,theta,phi} R_sum(S,X,W,theta,phi)';
problem.variables = {'S','X','W','theta','phi'};

problem.constraints = {
    'S \subseteq C, |S| = K_serv';
    'y_{n,m} - y_{n,m-1} >= Delta, m>=2';
    'x_n in X_n';
    'pi/2 <= theta_{n,m} <= pi';
    '-pi < phi_{n,m} <= pi';
    'trace(W*W^H) <= P_max'
    };

problem.feasible_set_Xn = 'X_n = {x_n in R^M | 0<=y_{n,1}<=...<=y_{n,M}<=Dy, y_{n,m}-y_{n,m-1}>=Delta}';
problem.params_snapshot = params;

end
