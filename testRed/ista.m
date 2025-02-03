function x = ista(A, b, L, lambda, L_f, max_iter, tol, x0)
    % ISTA for solving: min_x (1/2)||Ax - b||_2^2 + lambda * ||Lx||_1
    % Inputs:
    %   A        : Data matrix
    %   b        : Observed data
    %   L        : Diagonal matrix (or vector of diagonal elements)
    %   lambda   : Regularization parameter
    %   max_iter : Maximum number of iterations
    %   tol      : Convergence tolerance for stopping
    % Output:
    %   x        : Solution vector
    
    % Dimensions
    [m, n] = size(A);
    
    % Ensure L is treated as a vector if diagonal
    if ismatrix(L) && all(size(L) == [n, n])
        L = diag(L);
    end
    
    % Initialize variables
    x = x0;
    L_inv = 1 ./ L;       % Precompute inverse of diagonal entries
    
    % Precompute Lipschitz constant
    % L_f = norm(A' * A, 2); % Spectral norm of A'A
    eta = 1 / L_f;         % Step size (1/Lipschitz constant)
    
    % Iterative process
    for k = 1:max_iter
        % Store previous solution
        x_prev = x;
        
        % Gradient of the smooth part: f(x) = 1/2 ||Ax - b||_2^2
        grad_f = A' * (A * x - b);
        
        % Gradient descent step
        z = x - eta * grad_f;
        
        % Proximal step for g(x) = lambda * ||Lx||_1
        x = soft_threshold(z, eta * lambda * L_inv);
        
        % Check for convergence
        if norm(x - x_prev, 2) / max(1, norm(x_prev, 2)) < tol
            break;
        end
    end
    fprintf("Iterations of ISTA: %d\n",k)
end

% Soft-thresholding function
function x = soft_threshold(z, thresh)
    x = sign(z) .* max(abs(z) - thresh, 0);
end
