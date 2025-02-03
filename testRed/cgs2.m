function [x, iter, res] = cgs2(A, b, tol, maxIter, x0)
    % Conjugate Gradient Squared (CGS) Method
    % Inputs:
    %   A       - Coefficient matrix (sparse, non-symmetric)
    %   b       - Right-hand side vector
    %   tol     - Tolerance for stopping criterion (default: 1e-6)
    %   maxIter - Maximum number of iterations (default: 1000)
    % Outputs:
    %   x       - Solution vector
    %   iter    - Number of iterations performed
    %   res     - Residual norm ||b - Ax|| at each iteration
    
    % Check inputs
    if nargin < 3, tol = 1e-6; end
    if nargin < 4, maxIter = 1000; end
    
    % Initialization
    n = length(b);
    x = x0;                  % Initial guess
    r = b - A * x;           % Initial residual
    r_tld = r;               % Choose an arbitrary r_tld
    p = zeros(n, 1);         % Initial search direction
    u = zeros(n, 1);
    q = zeros(n, 1);
    v = zeros(n, 1);
    
    rho = 1;                 % Initial rho
    alpha = 1;               % Initial alpha
    omega = 1;               % Initial omega
    
    res = zeros(maxIter, 1); % Residual norms
    iter = 0;                % Iteration counter
    error = 1;

    xPrev = x0;
    % Iterate
    while error > tol && iter < maxIter
        iter = iter + 1;
        
        % Compute rho
        rho_old = rho;
        rho = r_tld' * r;
        if rho == 0
            error('Breakdown: rho = 0');
        end
        
        % Update direction vectors
        if iter == 1
            p = r;
        else
            beta = (rho / rho_old) * (alpha / omega);
            p = r + beta * (p - omega * v);
        end
        
        % Matrix-vector product
        v = A * p;
        
        % Compute alpha
        alpha = rho / (r_tld' * v);
        u = r - alpha * v;
        
        % Matrix-vector product for correction
        q = A * u;
        
        % Compute omega
        omega = (q' * u) / (q' * q);
        if omega == 0
            error('Breakdown: omega = 0');
        end
        
        % Update solution
        x = x + alpha * p + omega * u;
        
        % Update residual
        r = u - omega * q;
        
        % Store residual norm
        res(iter) = norm(r);

        error = norm(x - xPrev)/norm(xPrev);
        % error = norm(r)/norm(b);
        xPrev = x;
    end
    
    % Trim unused entries in residual array
    res = res(1:iter);
    
    if iter == maxIter
        warning('Maximum iterations reached before convergence.');
    end
end
