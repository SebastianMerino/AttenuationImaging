function [B,C] = optimAdmmWeightedTvTikhonov(A1,A2,b,mu1,mu2,m,n,tol,mask,W)
% Solver for SLD with Weighted Isotropic Total Variation regularization for ACS 
% and weighted Tikhonov regularization for BSC

p = length(mask)/(m*n);
minimask = reshape(mask,[m n p]);
minimask = minimask(:,:,1);
minimask = minimask(:);
b = mask.*b;
b(isnan(b)) = 0;
A = [A1 A2];
AtA = A'*A;
Atb = A'*b;
[u,~] = cgs(AtA,Atb);
B = reshape(u(1:end/2),m,n);
C = reshape(u(end/2+1:end),m,n);
%figure(109); imagesc(8.686*reshape(B,m,n)); colormap pink; caxis([0 1.2])
B = B(:);
C = C(:);
D = 0;
v = 0;

Wdiag = spdiags(W(:),0,m*n,m*n);
F(1) = 1/2*(norm( b - A1*B - A2*C ))^2 + ...
    mu1*TVcalc_isotropic(B,m,n,minimask) + mu2*sum(abs(Wdiag*C(:)),'all');

ite  = 0;
error = 1;

while abs(error) > tol && ite < 20
    ite = ite + 1;
    
    rho = 1;
    % First part of ADMM algorithm: B
    B = IRLS_TV_weighted(b-A2*C-D-v,A1,mu1/rho,m,n,tol,minimask,W);
    % B = IRLS_TV(b-A2*C-D-v,A1,mu1/rho,m,n,tol,mask,minimask);

    % Second part of ADMM algorithm: C
    Params.alpha2 = mu2/rho; Params.tolerance = tol;
    Params.beta = 0.1; Params.k = 1;
    % Params.operator = 'I';
    Params.operator = 'L'; Params.L = Wdiag;
    C = optimTikhonovReg(A2,b-A1*B-D-v,Params);
    
    % Third part of ADMM algorithm: D
    % Least squares: 1/2*||D||_2^2 + rho/2*||D-w||_2^2
    w = b - A1*B - A2*C - v;
    D = (rho/(rho+1))*w;
    
    % Fourth part of ADMM algorithm: v
    v = v + A1*B + A2*C + D - b;
    F(ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + ...
    mu1*TVcalc_isotropic(B,m,n,minimask) + mu2*sum(abs(Wdiag*C(:)),'all');
    
end

end


