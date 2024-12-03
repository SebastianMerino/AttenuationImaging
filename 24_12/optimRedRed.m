function [B,C,ite] = optimRedRed(A1,A2,b,muB,muC,m,n,tol,W)
% Regularized solver for SLD, RED for ACS and L1 norm for BSC

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
    muB*TVcalc_isotropic(B,m,n,W) + muC*sum(abs(Wdiag*C(:)),'all');

ite  = 0;
error = 1;
rho = 1;
median = 5;
Bprev = B;
while abs(error) > tol && ite < 100
    ite = ite + 1;
    
    % First part of ADMM algorithm: B
    B = optimRedLinear(A1,b-A2*C-D-v,muB*2/rho,tol,100,median,m,n,muB*2/rho,B);

    % Second part of ADMM algorithm: C
    C = optimRedLinear(A2,b-A1*B-D-v,muC*2/rho,tol,100,median,m,n,muC*2/rho,C);

    % Third part of ADMM algorithm: D
    % Least squares: 1/2*||D||_2^2 + rho/2*||D-w||_2^2
    w = b - A1*B - A2*C - v;
    D = (rho/(rho+1))*w;
    
    % Fourth part of ADMM algorithm: v
    v = v + A1*B + A2*C + D - b;
    F(ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + ...
    muB*TVcalc_isotropic(B,m,n,W) + muC*sum(abs(Wdiag*C(:)),'all');
   
    error = norm(B - Bprev)/norm(Bprev);
    Bprev = B;
end
% disp(F(ite+1))

end


function out = optimRedLinear(A,y,lambda,tol,max_iter,median,m,n,rho,x0)
x_est = x0;
v_est = x_est;
u_est = zeros(size(v_est));

AtArho = A'*A+rho*speye(length(x0));
Aty = A'*y;
xPrev = x_est;
for k = 1:max_iter
    % Part1 of the ADMM, approximates the solution of:
    % x = argmin_z 1/(2sigma^2)||Ax-y||_2^2 + 0.5*beta||x - v + u||_2^2
    [x_est,~] = cgs(AtArho, Aty + rho*(v_est-u_est), tol, 5,[],[],x_est);
    
    % Part2 of the ADMM, approximates the solution of
    % v = argmin_z lambda*z'*(z-denoiser(z)) +  0.5*beta||z - x - u||_2^2
    for j = 1:1
        v_reshap = reshape(v_est,m,n);
        f_v_est = medfilt2(v_reshap, [median median],'symmetric');
        f_v_est = f_v_est(:);
        v_est = (rho*(x_est + u_est) + lambda*f_v_est)/(lambda + rho);
    end

    % Part3 of the ADMM, update the dual variable
    u_est = u_est + x_est - v_est;

    if norm(xPrev - x_est)/norm(xPrev) < tol && k>2
        break
    end
    xPrev = x_est;

end
out = x_est;
% disp(k)
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




