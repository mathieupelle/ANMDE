function P_ortho = OrthoJacobiP(P, a, b)
    P_ortho = zeros(size(P));
    for i=1:length(P)
        n = i-1;
        gam = 2^(a+b+1)*gamma(n+a+1)*gamma(n+b+1)/(factorial(n)*(2*n+a+b+1)*gamma(n+a+b+1));
        P_ortho(i,1) = P(i,1)/sqrt(gam);
    end
end