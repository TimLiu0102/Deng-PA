function out = Constraint_Checker(mode, params, x_in)
% Constraint Checker：仅负责单条波导PA位置的可行性检查与投影

switch lower(mode)
    case 'check_position'
        out = check_position(params, x_in);

    case 'project_position'
        out = project_position(params, x_in);

    otherwise
        error('Constraint_Checker: unsupported mode');
end

end

%% ======================== 内部子函数 ========================
function is_feasible = check_position(params, x_in)
% PA位置物理可行性检查：
% 0 <= y1 <= ... <= yM <= Dy, 且 y_m - y_{m-1} >= Delta
x = x_in(:);
M = numel(x);
Dy = params.Dy;
Delta = params.Delta;

if M == 0
    is_feasible = false;
    return;
end

cond1 = all(x >= 0) && all(x <= Dy);
cond2 = all(diff(x) >= 0);
if M == 1
    cond3 = true;
else
    cond3 = all(diff(x) >= Delta);
end

is_feasible = cond1 && cond2 && cond3;
end

function x_proj = project_position(params, x_in)
% 将名义位置向量投影到论文可行域 chi_n（即 Pi_{chi_n}(x_in)）
% chi_n: 0<=y1<=...<=yM<=Dy, 且 y_m-y_{m-1}>=Delta
x = x_in(:);
M = numel(x);
Dy = params.Dy;
Delta = params.Delta;

if (M-1)*Delta > Dy
    error('Constraint_Checker: infeasible params, (M-1)*Delta > Dy');
end

% 变量替换：z_m = y_m - (m-1)Delta
% 则最小间距约束 y_m-y_{m-1}>=Delta 等价为 z_1<=z_2<=...<=z_M
% 同时由 0<=y_1, y_M<=Dy 可得 z 位于 [0, Dy-(M-1)Delta]
offset = ((0:M-1)' * Delta);
z_in = x - offset;

z_lb = 0;
z_ub = Dy - (M-1)*Delta;

% 对 z_in 做“单调非降 + 区间边界”的欧氏投影（bounded isotonic regression, PAVA）
z_proj = bounded_isotonic_projection(z_in, z_lb, z_ub);

% 映射回原变量 y_m = z_m + (m-1)Delta
x_proj = z_proj + offset;
end

function z_proj = bounded_isotonic_projection(z_in, z_lb, z_ub)
% PAVA：求解 min ||z-z_in||^2 s.t. z_1<=...<=z_M, z_lb<=z<=z_ub
% 每个块维护均值，若出现单调性违背则合并相邻块
v = z_in(:);
M = numel(v);

blk_start = zeros(M,1);
blk_end = zeros(M,1);
blk_sum = zeros(M,1);
blk_cnt = zeros(M,1);
blk_val = zeros(M,1);
B = 0;

for i = 1:M
    B = B + 1;
    blk_start(B) = i;
    blk_end(B) = i;
    blk_sum(B) = v(i);
    blk_cnt(B) = 1;
    blk_val(B) = min(max(blk_sum(B)/blk_cnt(B), z_lb), z_ub);

    while B >= 2 && blk_val(B-1) > blk_val(B)
        blk_end(B-1) = blk_end(B);
        blk_sum(B-1) = blk_sum(B-1) + blk_sum(B);
        blk_cnt(B-1) = blk_cnt(B-1) + blk_cnt(B);
        blk_val(B-1) = min(max(blk_sum(B-1)/blk_cnt(B-1), z_lb), z_ub);
        B = B - 1;
    end
end

z_proj = zeros(M,1);
for b = 1:B
    z_proj(blk_start(b):blk_end(b)) = blk_val(b);
end
end
