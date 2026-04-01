function varargout = Constraint_Checker(mode, x, params)
% 仅负责PA位置物理可行域检查与投影
% chi_n: 0<=y1<=...<=yM<=Dy, 且相邻间隔>=Delta

switch mode
    case 'project'
        varargout{1} = project_to_chi(x(:), params.Dy, params.Delta);
    case 'is_feasible'
        varargout{1} = is_feasible(x(:), params.Dy, params.Delta);
    otherwise
        error('Constraint_Checker: unknown mode');
end

end

function xproj = project_to_chi(x, Dy, Delta)
M = numel(x);
z = x - (0:M-1)'*Delta;
ub = Dy - (M-1)*Delta;
z = min(max(z,0),ub);
z = isotonic_nondec(z);
xproj = z + (0:M-1)'*Delta;
xproj = min(max(xproj,0),Dy);
end

function tf = is_feasible(x, Dy, Delta)
tf = all(x>=-1e-12) && all(x<=Dy+1e-12) && all(diff(x)>=Delta-1e-12);
end

function z = isotonic_nondec(y)
% 简洁PAV实现
z = y; n = numel(y); w = ones(n,1); i = 1;
while i < n
    if z(i) <= z(i+1)
        i = i + 1;
    else
        newv = (w(i)*z(i)+w(i+1)*z(i+1))/(w(i)+w(i+1));
        z(i) = newv; z(i+1) = newv;
        w(i) = w(i)+w(i+1); w(i+1)=w(i);
        j = i;
        while j>1 && z(j-1)>z(j)
            newv = (w(j-1)*z(j-1)+w(j)*z(j))/(w(j-1)+w(j));
            z(j-1)=newv; z(j)=newv;
            w(j-1)=w(j-1)+w(j); w(j)=w(j-1);
            j=j-1;
        end
        i = i + 1;
    end
end
for i = 2:n
    z(i) = max(z(i), z(i-1));
end
end
