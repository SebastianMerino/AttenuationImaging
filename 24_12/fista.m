function x = fista(A, b, L, lambda, L_f, max_iter, tol, x0)
    % FISTA for solving: min_x (1/2)||Ax - b||_2^2 + lambda * ||Lx||_1
    % Inputs:
    %   A        : Data matrix
    %   b        : Observed data
    %   L        : Diagonal matrix (can be a vector of diagonal elements)
    %   lambda   : Regularization parameter
    %   L_f      : Lipschitz constant (equal to frequency points)
    %   max_iter : Maximum number of iterations
    %   tol      : Convergence tolerance for stopping
    
    % Dimensions
    [m, n] = size(A);
    
    % Ensure L is treated as a vector if diagonal
    if ismatrix(L) && all(size(L) == [n, n])
        L = diag(L);
    end
    
    % Initialize variables
    x_k = x0;  % Initial solution
    y_k = x_k;          % Extrapolated variable
    t_k = 1;            % Momentum parameter
    L_inv = 1 ./ L;     % Precompute inverse of diagonal entries
    
    %  Lipschitz constant: Spectral norm of A'A 
    % L_f = norm(full(A' * A), 2); % EQUAL TO FREQUENCY POINTS
    eta = 1 / L_f; % L_f;         % Step size (1/Lipschitz constant)
    
    % Iterative process
    for k = 1:max_iter
        % Compute gradient of the smooth part (f(x) = 1/2||Ax - b||_2^2)
        grad_f = A' * (A * y_k - b);
        
        % Gradient descent step
        z_k = y_k - eta * grad_f;
        
        % Proximal step for g(x) = lambda * ||Lx||_1
        x_next = soft_threshold(z_k, eta * lambda * L_inv);
        
        % Update momentum parameter
        t_next = 0.5 * (1 + sqrt(1 + 4 * t_k^2));
        
        % Extrapolate
        y_next = x_next + ((t_k - 1) / t_next) * (x_next - x_k);
        
        % Check for convergence
        if norm(x_next - x_k, 2) / max(1, norm(x_k, 2)) < tol
            break;
        end
        
        % Update variables for next iteration
        x_k = x_next;
        y_k = y_next;
        t_k = t_next;
    end
    
    % Return the solution
    x = x_k;

    % fprintf("Iterations of FISTA: %d\n",k)
end

% Soft-thresholding function
function x = soft_threshold(z, thresh)
    x = sign(z) .* max(abs(z) - thresh, 0);
end
