function Tn = cheby_scale(pn,n)
    Tn = gamma(n+1)*gamma(0.5)/gamma(n+0.5) * pn;
end