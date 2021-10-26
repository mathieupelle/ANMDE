function [P] = JacobiP(x, alpha, beta, N)

    a = alpha; b = beta;
    P = zeros(N,1);
    P(1,1) = 1;
    if N>1
        P(2,1) = 0.5*(a-b+x*(a+b+2));
    end
    for i=2:N-1
        n = i-1;
        a_l = 2*(n+a)*(n+b)/((2*n+a+b+1)*(2*n+a+b));
        a_c = (a^2-b^2)/((2*n+a+b+2)*(2*n+a+b));
        a_r = 2*(n+1)*(n+a+b+1)/((2*n+a+b+2)*(2*n+a+b+1));
        P(i+1,1) = ((a_c + x)*P(i) - a_l*P(i-1))/a_r;
    end
end