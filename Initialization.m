function init = Initialization(scene, params)
% 按论文顺序初始化：C -> S^(0) -> X^(0) -> (theta,phi)^(0)

N=params.N; M=params.M; K=params.K;
Gpot = zeros(K,N,M); ystar = zeros(K,N,M); dstar = zeros(K,N,M);

% 1) Potential gain matrix 与候选池
for k = 1:K
    qk = scene.users(k,:);
    for n = 1:N
        Akn = (qk(1)-scene.xW(n))^2 + (qk(3)-params.d)^2;
        den = (2*log(params.alphaL))^2 - params.alphaW^2;
        gamma_star = sqrt(max(Akn*params.alphaW^2/max(den,1e-12), 0));
        for m = 1:M
            ystar(k,n,m) = qk(2) - gamma_star;
            p = [scene.xW(n), ystar(k,n,m), params.d];
            dstar(k,n,m) = norm(qk - p);
            Gpot(k,n,m) = sqrt(1/M)*exp(-(params.alphaW/2)*ystar(k,n,m)) ...
                * (params.alphaL^dstar(k,n,m)) ...
                * (params.lambda*params.n_refr*params.v*sqrt(2*params.a*params.b))/(2*max(dstar(k,n,m),1e-9));
        end
    end
end
Emax = squeeze(max(max(Gpot,[],3),[],2));
[~,idxE] = sort(Emax,'descend');
C = idxE(1:min(2*params.K_serv,K));

% 2) Global matching（收益最大匹配）
y_ref = zeros(N,M);
for n=1:N
    for m=1:M
        y_ref(n,m) = (m-1)*params.Delta + ((m-1)/max(M-1,1))*(params.Dy-(M-1)*params.Delta);
    end
end

U = -inf(numel(C), N*M);
for ic = 1:numel(C)
    k = C(ic);
    for n=1:N
        for m=1:M
            dmov = abs(y_ref(n,m)-ystar(k,n,m));
            U(ic,(n-1)*M+m) = Gpot(k,n,m) - params.lambda_mov*dmov;
        end
    end
end
match = max_weight_match(U, params.K_serv);

S = unique(C(match(:,1))); S = S(:).';
if numel(S) < params.K_serv
    extra = setdiff(C,S,'stable');
    S = [S, extra(1:min(params.K_serv-numel(S),numel(extra)))];
end
S = S(1:params.K_serv);

% 3) 初始X^(0)
X = y_ref;
owner = zeros(N,M);
for r = 1:size(match,1)
    k = C(match(r,1)); pa = match(r,2);
    n = ceil(pa/M); m = pa - (n-1)*M;
    X(n,m) = ystar(k,n,m); owner(n,m)=k;
end
for n=1:N
    X(n,:) = Constraint_Checker('project', X(n,:).', params).';
end

% 4) 初始角度
theta = pi*ones(N,M); phi = zeros(N,M);
for n=1:N
    for m=1:M
        if owner(n,m)~=0
            k = owner(n,m);
            p = [scene.xW(n), X(n,m), params.d];
            v = scene.users(k,:) - p;
            [th,ph] = vec_to_ang(v);
            theta(n,m)=th; phi(n,m)=ph;
        end
    end
end

init.S=S; init.X=X; init.theta=theta; init.phi=phi;
init.C=C; init.Emax=Emax; init.Gpot=Gpot;
init.matching=match; init.y_ref=y_ref; init.owner=owner;
end

function match = max_weight_match(U, Kserv)
% 行:候选用户, 列:PA槽位。输出[行,列]
[nr,nc]=size(U); pairs=[]; Uwork=U;
for t=1:min([Kserv,nr,nc])
    [mx,id]=max(Uwork(:));
    if ~isfinite(mx), break; end
    [r,c]=ind2sub([nr,nc],id);
    pairs=[pairs; r,c]; %#ok<AGROW>
    Uwork(r,:)=-inf; Uwork(:,c)=-inf;
end
match=pairs;
end

function [theta,phi] = vec_to_ang(v)
r = norm(v);
if r<1e-12
    theta=pi; phi=0; return;
end
theta = acos(max(min(v(3)/r,1),-1));
theta = min(max(theta,pi/2),pi);
phi = atan2(v(2),v(1));
if phi<=-pi, phi=phi+2*pi; end
if phi>pi, phi=phi-2*pi; end
end
