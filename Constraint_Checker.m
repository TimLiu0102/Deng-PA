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
% 将原始位置向量投影到可行集 chi_n
% chi_n: 0<=y1<=...<=yM<=Dy, 且相邻间距>=Delta
x = x_in(:);
M = numel(x);
Dy = params.Dy;
Delta = params.Delta;

if (M-1)*Delta > Dy
    error('Constraint_Checker: infeasible params, (M-1)*Delta > Dy');
end

if M == 1
    x_proj = min(max(x,0),Dy);
    return;
end

% 1) 基本区间裁剪
x = min(max(x, 0), Dy);

% 2) 前向修正：保证 y_m >= y_{m-1} + Delta
for m = 2:M
    x(m) = max(x(m), x(m-1) + Delta);
end

% 3) 末端裁剪并反向修正：保证 y_M <= Dy 且保持最小间距
x(M) = min(x(M), Dy);
for m = M-1:-1:1
    x(m) = min(x(m), x(m+1) - Delta);
end

% 4) 再做一次区间与前向修正，确保最终满足chi_n
x = min(max(x, 0), Dy);
for m = 2:M
    x(m) = max(x(m), x(m-1) + Delta);
end

% 若末端超界，则整体平移回界内（保持间距结构）
if x(M) > Dy
    shift = x(M) - Dy;
    x = x - shift;
    x = min(max(x,0),Dy);
end

x_proj = x;
end
