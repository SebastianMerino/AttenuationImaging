
% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function [u,G] = IRLS_ANIS_TV_weighted(b,A,mu,M,N,tol,mask,minimask,W)

[u,~] = cgs(A'*A,A'*b,1e-6,20);
%figure(109); imagesc(8.686*reshape(u,M,N)); colormap pink; caxis([0 1.2])

G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_anisotropic(u,M,N,W);

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dx = kron(speye(N),D);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dy = kron(D,speye(M));

D = [Dx' Dy']';

ite_irls = 0;
error = 1;

while error > tol && ite_irls < 20
    
    X = reshape(u,M,N);
    ite_irls = ite_irls + 1;
    Dh = [diff(X,[],1);zeros(1,N)];
    Dv = [diff(X,[],2),zeros(M,1)];

    % WEIGHTING DIFFERENCES
%     figure,imagesc(Dh)
%     figure,imagesc(Dv)
    Dh = Dh.*W;
    Dv = Dv.*W;
%     figure,imagesc(Dh)
%     figure,imagesc(Dv)

    Px = abs(Dh + eps);
    Px = 1./Px;
    Px = Px(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(Px,0,omega);
    Wx = kron(speye(1),omega);
    
    Py = abs(Dv + eps);
    Py = 1./Py;
    Py = Py(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(Py,0,omega);
    Wy = kron(speye(1),omega);
    
    AtA = A'*A;
    %mu=5000;
    %[u] = cgs(AtA + mu*D'*W*D, A'*b,1e-6,200);
    [u,~] = cgs( AtA + mu*Dx'*Wx*Dx + mu*Dy'*Wy*Dy , A'*b, 1e-6 , 20);
    
    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_anisotropic(u,M,N,W);
    error = abs(G(ite_irls+1) - G(ite_irls));

%     figure,imagesc(reshape(u*8.686,M,N),[0.5 1.2])

end

%figure(909); plot(1:length(G),G);

end


