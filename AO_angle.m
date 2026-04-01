function [theta, phi] = AO_angle(S, X, W, theta, phi, scene, params)
% 论文角度子问题：pattern search + re-anchoring
D = [1 0; -1 0; 0 1; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
N=params.N; M=params.M;

for n=1:N
    for m=1:M
        th_c = theta(n,m); ph_c = phi(n,m);
        dth = params.Delta_theta0; dph = params.Delta_phi0;

        for r = 1:params.I_theta
            R_center = Signal_model('sum_rate', S, X, W, theta, phi, scene, params);
            bestR = -inf; best = [th_c, ph_c];

            for ii = 1:size(D,1)
                cand = project_angle([th_c + D(ii,1)*dth, ph_c + D(ii,2)*dph]);
                theta_t = theta; phi_t = phi;
                theta_t(n,m)=cand(1); phi_t(n,m)=cand(2);
                R_cand = Signal_model('sum_rate', S, X, W, theta_t, phi_t, scene, params);
                if R_cand > bestR
                    bestR = R_cand; best = cand;
                end
            end

            if bestR >= R_center + params.eps_theta
                th_c = best(1); ph_c = best(2);
                theta(n,m)=th_c; phi(n,m)=ph_c;
            else
                anchors = make_anchors(n,m,X,scene,params);
                [th_a, ph_a] = choose_best_anchor(anchors,n,m,S,X,W,theta,phi,scene,params);
                th_c = th_a; ph_c = ph_a; theta(n,m)=th_c; phi(n,m)=ph_c;
                dth = params.beta_theta*dth; dph = params.beta_phi*dph;
            end
        end
    end
end

end

function a = project_angle(a)
a(1) = min(max(a(1),pi/2),pi);
a(2) = atan2(sin(a(2)), cos(a(2)));
if a(2)<=-pi, a(2)=a(2)+2*pi; end
end

function anchors = make_anchors(n,m,X,scene,params)
anchors = [pi,0];
p = [scene.xW(n), X(n,m), params.d];
rep = 1:min(3,size(scene.users,1));
for k = rep
    v = scene.users(k,:) - p;
    r = norm(v);
    if r<1e-12, th=pi; ph=0; else
        th = acos(max(min(v(3)/r,1),-1));
        th = min(max(th,pi/2),pi);
        ph = atan2(v(2),v(1));
    end
    anchors = [anchors; project_angle([th,ph])]; %#ok<AGROW>
end
end

function [th,ph] = choose_best_anchor(anchors,n,m,S,X,W,theta,phi,scene,params)
bestR=-inf; th=theta(n,m); ph=phi(n,m);
for i=1:size(anchors,1)
    theta_t=theta; phi_t=phi;
    theta_t(n,m)=anchors(i,1); phi_t(n,m)=anchors(i,2);
    R=Signal_model('sum_rate',S,X,W,theta_t,phi_t,scene,params);
    if R>bestR, bestR=R; th=anchors(i,1); ph=anchors(i,2); end
end
end
