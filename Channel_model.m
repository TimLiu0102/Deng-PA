function out = Channel_model(mode, S, X, theta, phi, W, params)
% 仅负责系统几何与信道模型（论文系统模型部分）

switch mode
    case 'build_scene'
        out = build_scene(params);
    case 'composite_channels'
        out = composite_channels(S, X, theta, phi, params);
    otherwise
        error('Channel_model: unknown mode');
end

end

function scene = build_scene(params)
N = params.N; K = params.K;
xW = ((2*(1:N)-1)/(2*N))*params.Dx; % 波导x坐标
feed = [xW(:), zeros(N,1), params.d*ones(N,1)];

qx = params.user_x_rng(1) + (params.user_x_rng(2)-params.user_x_rng(1))*rand(K,1);
qy = params.user_y_rng(1) + (params.user_y_rng(2)-params.user_y_rng(1))*rand(K,1);
qz = params.user_z_rng(1) + (params.user_z_rng(2)-params.user_z_rng(1))*rand(K,1);
users = [qx, qy, qz];

scene.xW = xW(:); scene.feed = feed; scene.users = users;
end

function H = composite_channels(S, X, theta, phi, params)
% 输出H: (NM) x |S|，每列对应用户k的复合信道h_k(X)
N = params.N; M = params.M; NM = N*M;
scene = params.scene;
Kserv = numel(S);
H = zeros(NM, Kserv);

for idx = 1:Kserv
    k = S(idx);
    qk = scene.users(k,:).';
    hk = zeros(NM,1);
    for n = 1:N
        xn = X(n,:).';
        gn = guided_channel(xn, params);
        hfree = free_space_channel(qk, n, xn, theta(n,:).', phi(n,:).', scene, params);
        hk((n-1)*M+1:n*M) = gn .* hfree;
    end
    H(:,idx) = hk;
end
end

function gn = guided_channel(xn, params)
lambda_g = params.lambda/params.n_eff;
M = params.M;
gn = sqrt(1/M) * exp(-(params.alphaW/2)*xn) .* exp(-1j*2*pi/lambda_g*xn);
end

function htilde = free_space_channel(qk, n, xn, thetan, phin, scene, params)
M = params.M;
htilde = zeros(M,1);
for m = 1:M
    pnm = [scene.xW(n); xn(m); params.d];
    dkm = norm(qk - pnm);
    xyz_local = global_to_local(qk - pnm, thetan(m), phin(m));
    ups = local_pattern(xyz_local(1), xyz_local(2), xyz_local(3), params);
    htilde(m) = sqrt(params.eta * (params.alphaL^dkm)) * ups;
end
end

function xyz_tilde = global_to_local(v, theta, phi)
Rz = [cos(pi/2-phi), -sin(pi/2-phi), 0; sin(pi/2-phi), cos(pi/2-phi), 0; 0,0,1];
Rx = [1,0,0; 0,cos(theta-pi/2),-sin(theta-pi/2); 0,sin(theta-pi/2),cos(theta-pi/2)];
xyz_tilde = Rx * Rz * v;
end

function val = local_pattern(x, y, z, params)
% 局部方向图 Upsilon(x,y,z)
y_eff = max(y, 1e-6);
k0 = 2*pi/params.lambda;
w1 = params.v*params.a*params.lambda;
w2 = params.v*params.b*params.lambda;
W1 = params.lambda*y_eff/(pi*params.n_refr*w1);
W2 = params.lambda*y_eff/(pi*params.n_refr*w2);
R1 = y_eff; R2 = y_eff;
Th1 = atan(params.lambda*y_eff/(pi*params.n_refr*w1));
Th2 = atan(params.lambda*y_eff/(pi*params.n_refr*w2));

amp = sqrt((w1*w2)/(W1*W2))*params.B_norm;
phase = -(x^2/W1^2 + z^2/W2^2) ...
    -1j*k0*params.n_refr*(x^2/(2*R1) + z^2/(2*R2) + y_eff) ...
    +1j*(Th1+Th2)/2;
val = amp*exp(phase);
end
