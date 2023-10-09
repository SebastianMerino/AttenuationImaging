function [B,C] = AlterOpti_ADMM(A1,A2,b,mu1,mu2,m,n,tol,mask)

p = length(mask)/(m*n);
minimask = reshape(mask,[m n p]);
minimask = minimask(:,:,1);
minimask = minimask(:);
b = mask.*b;
b(isnan(b)) = 0;
A = [A1 A2];
[u,~] = cgs2(A'*A,A'*b,1e-6,20);
B = reshape(u(1:end/2),m,n);
C = reshape(u(end/2+1:end),m,n);

B = B(:);
C = C(:);
D = 0;
v = 0;

F(1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc(B,m,n,minimask) + mu2*TVcalc(C,m,n,minimask);

ite  = 0;
error = 1;

while abs(error) > tol && ite < 20
    ite = ite + 1;

    rho = 1;
    % First part of ADMM algorithm: B
    B = IRLS_TV(b-A2*C-D-v,A1,mu1/rho,m,n,tol,minimask);

    % Second part of ADMM algorithm: C
    C = IRLS_TV(b-A1*B-D-v,A2,mu2/rho,m,n,tol,minimask);

    % Third part of ADMM algorithm: D
    % Least squares: 1/2*||D||_2^2 + rho/2*||D-w||_2^2
    w = b - A1*B - A2*C - v;
    D = (rho/(rho+1))*w;

    % Fourth part of ADMM algorithm: v
    v = v + A1*B + A2*C + D - b;
    F(ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + ...
        mu1*TVcalc(B,m,n,minimask) + mu2*TVcalc(C,m,n,minimask);
end

end


function [TV] = TVcalc(B,M,N,mask)
% Returns the Total Variation
% Inputs: 
%       B               vector containing an image
%       M,N             image size
%       mask            binary mask of analyzed region
%  
% Outputs:
%       TV              number
%
% Author: Andres Leonel Coila

mask(isnan(mask)) = 0;
mask = mask(:);

X = reshape(B,M,N);
Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

P = Dh.^2 + Dv.^2;
P = sqrt(P);
TV = norm(P(:).*mask,1);

end

function [u,C] = IRLS_TV(b,A,mu,M,N,tol,mask)
% Returns the solution to a linear system regularized with Total Variation
%   COST = 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
% Inputs: 
%       b               vector containing measurements
%       A               matrix for linear system of eq
%       mu              regularization parameter
%       M,N             image size of u
%       tol             tolerance for relative error
%       mask            binary mask of analyzed region. Vector from 2D
%  
% Outputs:
%       u               vector of image samples, size MN
%       C               vector containing cost function for each iteration
%
% Author: Andres Leonel Coila
% Modified by Sebastian Merino

AtA = A'*A; % This takes A LOT OF TIME
Atb = A'*b;
    
[u,~] = cgs2(AtA,Atb,1e-6,20);
C(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask);

[Dx,Dy] = diffOperator(M,N,true);
D = [Dx' Dy']';

ite_irls = 0;
error = 1;
eps = 0.01; % 0.01 for isotropic is ok
while error > tol        
    ite_irls = ite_irls + 1;

    Dh = Dx*u;
    Dv = Dy*u;
    P = sqrt(Dh.^2 + Dv.^2 + eps^2);
    P = P.^(-0.5);
    P = P(:).*mask;
    
    omega = speye(M*N);
    omega = spdiags(P,0,omega);
    
    W = kron(speye(2),omega);
    [u,~] = cgs2( AtA + mu*D'*W*D, Atb, tol , 20, u );
    
    C(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,mask);
    error = abs(C(ite_irls+1) - C(ite_irls));
end
end

function [Dx,Dy] = diffOperator(M,N,colMajor)
% Returns differential operators
% Inputs: 
%       M,N             image size
%       colMajor        true if the matrix is stored col-major, false if
%                       stored row-major
%  
% Outputs:
%       Dx,Dy           Difference operators
% Author: Sebastian Merino
if ~colMajor
    foo = M;
    M = N;
    N = foo;
end

G = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
G(:,end) = [];
G(M,M) = 0;
Dx = kron(speye(N),G);

G = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
G(:,end) = [];
G(N,N) = 0;
Dy = kron(G,speye(M));

end

function [x,ite_cgs] = cgs2(A,b,tol,maxit,varargin)
% Solves A*x = b using the conjugate squared gradient method

if length(varargin) == 1
    x = varargin{1};
else
    x= zeros(size(A,2),1);
end
r = A*x - b;
p = -r;
ite_cgs = 0;

while norm(r,2) > tol && ite_cgs < maxit

    alpha = (r'*r)/(p'*(A*p));
    x = x + alpha*p;
    rn = r + alpha*(A*p);
    beta = (rn'*rn)/(r'*r);
    p = -rn + beta*p;
    r = rn;
    ite_cgs = ite_cgs + 1;
end

end
