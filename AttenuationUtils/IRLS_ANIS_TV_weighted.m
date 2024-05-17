
% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function [u,G] = IRLS_ANIS_TV_weighted(b,A,mu,M,N,tol,mask,minimask,W)

[u,~] = cgs(A'*A,A'*b,1e-6,20);
%figure(109); imagesc(8.686*reshape(u,M,N)); colormap pink; caxis([0 1.2])

G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_anisotropic(u,M,N,W);
Wdiag = spdiags(W(:),0,M*N,M*N);

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dy = kron(speye(N),D);
Dy = Wdiag*Dy;

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dx = kron(D,speye(M));
Dx = Wdiag*Dx;

ite_irls = 0;
error = 1;

while error > tol && ite_irls < 20
    
    ite_irls = ite_irls + 1;
    
    Dh = Dx*u;
    Dv = Dy*u;
    
    eps = 0.1;
    
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
    [u,~] = pcg( AtA + mu*Dx'*Wx*Dx + mu*Dy'*Wy*Dy , A'*b, 1e-6 , 20, [], [], u);
    
    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_anisotropic(u,M,N,W);
    error = abs(G(ite_irls+1) - G(ite_irls));
    % figure, imagesc(reshape(u,M,N)); % NUNCA DESCOMENTAR O CRASHEARA
end

end


