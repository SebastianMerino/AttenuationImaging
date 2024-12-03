function [res_admm,n_out , out ] = admm_red_median(A,y,lambda,error,N,max_iter,inner_iters,inner_iters2,median,m,ni,beta)

x_est= ones([N 1]);
v_est=ones([N 1]);
u_est=ones([N 1]);

n_out = [x_est v_est u_est];
alpha=1.5;
for k = 1:1:max_iter

    % Part1 of the ADMM, approximates the solution of:
    % x = argmin_z 1/(2sigma^2)||Ax-y||_2^2 + 0.5*beta||x - v + u||_2^2
    for j = 1:1:inner_iters
        b = (A*y) + beta*(v_est - u_est);
        A_x_est = (A'*A*x_est) + beta*x_est;
        res = b - A_x_est;
        a_res = (A'*A*res) + beta*res;
        mu_opt = mean(res(:).*res(:))/mean(res(:).*a_res(:));
        x_est = x_est + mu_opt*res;
    end
    
    % relaxation
    x_hat = alpha*x_est + (1-alpha)*v_est;

    % Part2 of the ADMM, approximates the solution of
    % v = argmin_z lambda*z'*(z-denoiser(z)) +  0.5*beta||z - x - u||_2^2
    % using gradient descent
    for j = 1:1:inner_iters2
        v_est = reshape(x_est,m,2*ni);
        v_est1=v_est((1:m),(1:ni));
        v_est2=v_est((1:m),(ni+1:2*ni));
        f_v_est1 = medfilt2(v_est1, [median median],'symmetric');
        f_v_est2 = medfilt2(v_est2, [median median],'symmetric');

        f_v_est = [f_v_est1 f_v_est2];
        f_v_est = f_v_est(:);

        v_est = (beta*(x_hat + u_est) + lambda*f_v_est)/(lambda + beta);
    end

    % Part3 of the ADMM, update the dual variable
    u_est = u_est + x_hat - v_est;

    r(k) = norm(A*x_hat-y)^2/norm(y)^2;

    if r(k) < error
        break
    end


end
out = x_hat;
% n_out=k;
res_admm = r;
end