function [u,G] = IRLS_TV_weighted(b,A,mu,M,N,tol,mask,weights)
% Optimizes the following cost function of weighted isotropic TV:
%   0.5*||A*u(:)-b||_2^2 + mu*SWTV(u)
% Inputs: 
%       b               vector containing measurements
%       A               matrix for linear system of eq
%       mu              regularization parameter
%       M,N             image size of u
%       tol             tolerance for error
%       mask            binary mask of analyzed region
%       weights         weights as a MxN matrix
%  
% Outputs:
%       u               vector of image samples, size MN
%       G               vector containing cost function for each iteration
%
% Author: A. Coila
% Weights and docs added by Sebastian Merino

AtA = A'*A;
Atb = A'*b;

Wdiag = spdiags(weights(:),0,M*N,M*N);

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

D = [Dx' Dy']';

ite_irls = 0;
error = 1;

%[u,~] = cgs(AtA+mu*(D')*D,Atb);
[u,~] = cgs(AtA,Atb);
G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(u,M,N,weights);

while error > tol && ite_irls < 200
    
    X = reshape(u,M,N);
    ite_irls = ite_irls + 1;
    Dh = Dx*u;
    Dv = Dy*u;
    
    vksquare = Dh.^2 + Dv.^2;
    vksquare = vksquare(:);
    
    eps = 0.3;
    P = sqrt(vksquare + eps^2);
    P = 1./P;
    
    %P = sqrt(P.^2 + eps^2);
    P = P(:).*mask;
    omega = speye(M*N);   % sparse diagonal of just ones
    omega = spdiags(P,0,omega);   % sparse diagonal of P values instead of ones.
    W = kron(speye(2),omega);
    
    [u,~] = cgs(AtA + mu*D'*W*D, Atb, 1e-6, 200);
    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(u,M,N,weights);
    error = abs(G(ite_irls+1) - G(ite_irls));
    %error = (G(ite_irls+1) - G(ite_irls)).^2/G(ite_irls).^2;
    
end

%figure(909); plot(1:length(G),G);

end
