function [L2] = L2norm(f, x)
    % Computes L2 norm of function using matrix-based integration.
    % f: Value of function at nodes
    % x: Node locations in range [-1, 1] (Gauss, Gauss Lobato...)
    % L2: L2 norm

    N = length(x);
    V = zeros(N, N); 
    for i=1:N
        V(:,i) = JacobiP(x(i), 0, 0, N);
    end

    M = (V'*V)^(-1); % Mass matrix
    L2 = sqrt(f'*M*f); 
end