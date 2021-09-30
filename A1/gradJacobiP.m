function [dPdx] = gradJacobiP(x, alpha, beta, N)
    
    dPdx = zeros(N, 1);
    P = JacobiP(x, alpha+1, beta+1, N-1);
    for i=1:N
        n = i-1;
        if i == 1
            dPdx(i,1) = 0;
        else
            dPdx(i,1) = sqrt(n*(n+alpha+beta+1))*P(i-1,1);
            dPdx(i,1) = 0.5*(alpha+beta+1+n)*P(i-1,1);
        end
    end
    
end