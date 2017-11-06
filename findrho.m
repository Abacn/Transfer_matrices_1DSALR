% Find density, given temperature and pressure

function rho = findrho(p, beta, coeffs, N)
    % start from hard sphere
    if(nargin<4)
        N=200;
    end
    rho = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, N))),p);
end

