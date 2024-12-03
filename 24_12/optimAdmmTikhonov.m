function x = optimAdmmTikhonov(A, b, L, lambda, rho, max_iters, tol, x0)
    % Precompute matrix inversions for efficiency
    AtA = A' * A;
    Ab = A' * b;
    P = AtA + rho * (L'*L);
    % M = ichol(P);

    % Initialize variables
    x = x0;
    z = L*x;
    u = zeros(size(z));

    for k = 1:max_iters
        % x-update using Cholesky factorization
        q = Ab + rho * L' * (z - u);
        [x,~] = cgs(P,q,tol,200, [], [],x);

        % z-update with soft-thresholding
        Lx = L * x;
        z_old = z;
        z = soft_threshold(Lx + u, lambda / rho);

        % u-update
        u = u + (Lx - z);

        % Convergence check
        if norm(z - z_old)/norm(z_old) < tol
            break;
        end
    end
    fprintf("Iterations of ADMM-L1: %d\n",k)
end

% Soft-thresholding function
function z = soft_threshold(x, threshold)
    z = sign(x) .* max(abs(x) - threshold, 0);
end
