function [B,C,ite] = optimSwift(A1,A2,b,mu1,mu2,m,n,tol,W)
% Solver for SLD with Weighted Isotropic Total Variation regularization for ACS 
% and weighted Tikhonov regularization for BSC

p = length(b)/(m*n);
b(isnan(b)) = 0;
A = [A1 A2];
AtA = A'*A;
Atb = A'*b;
[u,~] = cgs(AtA,Atb);
B = reshape(u(1:end/2),m,n);
C = reshape(u(end/2+1:end),m,n);

B = B(:);
C = C(:);
D = 0;
v = 0;

Wdiag = spdiags(W(:),0,m*n,m*n);
F(1) = 1/2*(norm( b - A1*B - A2*C ))^2 + ...
    mu1*TVcalc_isotropic(B,m,n,W) + mu2*sum(abs(Wdiag*C(:)),'all');

ite  = 0;
error = 1;
rho = 1;
Bprev = B;
while abs(error) > tol && ite < 100
    ite = ite + 1;
    
    % First part of ADMM algorithm: B
    B = IRLS_TV_weighted(b-A2*C-D-v,A1,mu1/rho,m,n,tol,W, B);

    % Second part of ADMM algorithm: C
    Params.alpha2 = mu2/rho; Params.tolerance = tol;
    Params.beta = 0.1; Params.k = 1;
    Params.operator = 'L'; Params.L = Wdiag;
    C = optimTikhonovReg(A2,b-A1*B-D-v,Params,C);
    % C = fista(A2, b-A1*B-D-v, Wdiag, mu2/rho, p, 200, tol, C);
    % C = ista(A2, b-A1*B-D-v, Wdiag, mu2/rho, p, 200, tol, C);
    % C = optimAdmmTikhonov(A2, b-A1*B-D-v, Wdiag, mu2/rho, 10, 200, tol, C);

    % Third part of ADMM algorithm: D
    % Least squares: 1/2*||D||_2^2 + rho/2*||D-w||_2^2
    w = b - A1*B - A2*C - v;
    D = (rho/(rho+1))*w;
    
    % Fourth part of ADMM algorithm: v
    v = v + A1*B + A2*C + D - b;
    F(ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + ...
    mu1*TVcalc_isotropic(B,m,n,W) + mu2*sum(abs(Wdiag*C(:)),'all');
   
    error = norm(B - Bprev)/norm(Bprev);
    Bprev = B;
end
% disp(F(ite+1))

end


function [u,G] = IRLS_TV_weighted(b,A,mu,M,N,tol,weights, u0)
% Optimizes the following cost function of weighted isotropic TV:
%   0.5*||A*u(:)-b||_2^2 + mu*SWTV(u)
% Inputs: 
%       b               vector containing measurements
%       A               matrix for linear system of eq
%       mu              regularization parameter
%       M,N             image size of u
%       tol             tolerance for error
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
relativeChange = NaN;

u = u0;
uPrev = u0;
G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(u,M,N,weights);

% setup = struct('type','ilutp','droptol',1e-6);
% [L,U] = ilu(AtA + mu*(D'*D),setup);
% M_func = @(x) x ./ diag(A);

% fprintf("CGS Iterations: ")
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
    P = P(:);
    omega = speye(M*N);   % sparse diagonal of just ones
    omega = spdiags(P,0,omega);   % sparse diagonal of P values instead of ones.
    W = kron(speye(2),omega);
    
    [u,~,~,innerIter] = cgs(AtA + mu*D'*W*D, Atb, tol, 200, [], [],u);
    % [u, info] = pcg_solver(AtA + mu*D'*W*D, Atb, M_func, u, tol, 200);
    % [u,~,~,innerIter] = cgs(AtA + mu*D'*W*D, Atb, tol, 200,L,U,u);
    % fprintf("%d,", innerIter)

    G(ite_irls+1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc_isotropic(u,M,N,weights);
    relativeChange(ite_irls+1) = norm(u - uPrev)/norm(uPrev);
    
    % error = abs(G(ite_irls+1) - G(ite_irls));    
    error = norm(u - uPrev)/norm(uPrev);
    uPrev = u;

end
% fprintf("\nIterations of IRLS: %d\n",ite_irls)
% figure; tiledlayout(2,1),
% nexttile,
% plot(0:ite_irls,G);
% title('Loss function')
% xlim([1,ite_irls])
% nexttile,
% plot(0:ite_irls,relativeChange);
% title('Relative change in solution')

end


function [x,ite] = optimTikhonovReg(A,B,Params,x0)
% Optimizes the system Ax = B using the generalized Tikhonov regularization
%   ||Ax-b||^2 + alpha2 * ||Lx||_k
switch (Params.operator)
    case 'I'        % Penalizing greater values
        L=speye(size(A,2));
    case 'G'        % Penalizing greater gradients
        L=speye(size(A,2))-circshift(speye(size(A,2)),[0 1]);
        L(end,:)=0;
    case 'L'        % Penalizing customized matrix
        L = Params.L;
end

% [x0,~] = cgs(A'*A,A'*B,Params.tolerance,200);
err = 1;
ite = 0;
% fprintf("CGS Iterations: ")
while err > Params.tolerance && ite < 200
    Lx = L*x0;
    W = spdiags( Params.k/2*( abs(Lx.^2+Params.beta).^(Params.k/2 - 1) ),...
        0, length(Lx), length(Lx));
    % x = ((A'*A+Params.alpha2 *L'*W*L)\A') *B;
    [x,~,~,innerIter] = cgs( A'*A+Params.alpha2 *L'*W*L , A'*B, 1e-6 , 200,[],[],x0);
    % [x,innerIter] = cgs2( A'*A+Params.alpha2 *L'*W*L , A'*B, 1e-6 , 200,x0);
    % fprintf("%d,", innerIter)
    err = norm(x-x0)/norm(x); 
    x0 = x;
    ite = ite + 1;
end

% fprintf("\nIterations of Tik: %d\n",ite)


end


function [x, info] = pcg_solver(A, b, M, x0, tol, max_iter)
% PCG_SOLVER Solve Ax = b using the Preconditioned Conjugate Gradient (PCG) method.
%
% Inputs:
%   A        - Function handle or matrix representing the linear system.
%   b        - Right-hand side vector.
%   M        - Function handle or matrix representing the preconditioner (optional).
%   x0       - Initial guess for the solution (optional, defaults to zero vector).
%   tol      - Convergence tolerance (optional, defaults to 1e-6).
%   max_iter - Maximum number of iterations (optional, defaults to 1000).
%
% Outputs:
%   x        - Approximate solution vector.
%   info     - Struct containing convergence information.
%              Fields: converged, iterations, residual_norm.

    % Set default values for optional inputs
    if nargin < 3 || isempty(M)
        M = @(x) x; % Identity preconditioner
    elseif ~isa(M, 'function_handle')
        M = @(x) M \ x; % Convert matrix preconditioner to function handle
    end
    
    if nargin < 4 || isempty(x0)
        x0 = zeros(size(b));
    end
    
    if nargin < 5 || isempty(tol)
        tol = 1e-6;
    end
    
    if nargin < 6 || isempty(max_iter)
        max_iter = 1000;
    end
    
    % Ensure A is a function handle
    if ~isa(A, 'function_handle')
        A = @(x) A * x; % Convert matrix A to a function handle
    end
    
    % Initialize variables
    x = x0;
    r = b - A(x);           % Initial residual
    z = M(r);               % Apply preconditioner
    p = z;                  % Initial search direction
    rz_old = r' * z;        % Dot product of r and z
    
    for k = 1:max_iter
        Ap = A(p);
        alpha = rz_old / (p' * Ap);
        x = x + alpha * p;   % Update solution
        r = r - alpha * Ap;  % Update residual
        
        % Check for convergence
        residual_norm = norm(r);
        if residual_norm < tol
            info.converged = true;
            info.iterations = k;
            info.residual_norm = residual_norm;
            return;
        end
        
        z = M(r);            % Apply preconditioner
        rz_new = r' * z;     % Dot product of r and z (new)
        beta = rz_new / rz_old;
        p = z + beta * p;    % Update search direction
        rz_old = rz_new;
    end
    
    % If max iterations reached
    info.converged = false;
    info.iterations = max_iter;
    info.residual_norm = residual_norm;
end





