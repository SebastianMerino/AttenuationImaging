function [res_admm, out] = admmRedMixedDenoiser(A,y,lambda,tol,N,max_iter,median,m,ni,rho)

[x_est,~] = cgs(A'*A, A'*y, tol, 100);
v_est = x_est;
u_est=zeros([N 1]);
% x_est= ones([N 1]);
% v_est=ones([N 1]);
% u_est=ones([N 1]);


AtArho = A'*A+rho*speye(N);
Aty = A'*y;

xPrev = x_est;

for k = 1:1:max_iter

    % Part1 of the ADMM, approximates the solution of:
    % x = argmin_z 1/(2sigma^2)||Ax-y||_2^2 + 0.5*beta||x - v + u||_2^2
    [x_est,~] = cgs(AtArho, Aty + rho*(v_est-u_est), tol, 5,[],[],x_est);
    
    % Part2 of the ADMM, approximates the solution of
    % v = argmin_z lambda*z'*(z-denoiser(z)) +  0.5*beta||z - x - u||_2^2
    % using gradient descent
    for j = 1:5
        v_reshap = reshape(v_est,m,2*ni);
        v_est1=v_reshap((1:m),(1:ni));
        v_est2=v_reshap((1:m),(ni+1:2*ni));
        f_v_est1 = medfilt2(v_est1, [median median],'symmetric');
        f_v_est2 = wthresh(v_est2, 's', 0.001);

        f_v_est = [f_v_est1 f_v_est2];
        f_v_est = f_v_est(:);

        v_est = (rho*(x_est + u_est) + lambda*f_v_est)/(lambda + rho);
    end

    % Part3 of the ADMM, update the dual variable
    u_est = u_est + x_est - v_est;

    r(k) = norm(A*x_est-y)^2/norm(y)^2;

    if norm(xPrev - x_est)/norm(xPrev) < tol && k>2
        break
    end
    xPrev = x_est;

end
out = x_est;
% n_out=k;
res_admm = r;
disp(k)
end