function results = voiddistribution(T, p, coeffs, N)
    if(nargin < 4)
        N = 100;
    end
    ls = coeffs(1)*coeffs(2);
    xx = linspace(ls,10,N).';
    pots = potential(ls, coeffs);
    yy = exp(-(pots+p*xx)/T);
    results = [xx ,yy];
end