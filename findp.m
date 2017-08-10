% Find pressure, given temperature and density

function p = findp(rho, T, coeffs, tol)
    % start from hard sphere
    if(nargin < 4)
        tol  = 1e-6;
    end
    sigma = coeffs(1);
    beta = 1/T;
    p0 = rho*T/(1-sigma*rho);
    rho0 = beta/deriv(@(bp)(-log(Pfunc_isobaric_NN(bp, beta, coeffs, 200))),p0);
    % First search a border
    if(rho0 < rho)
        while(rho0 < rho)
            p0 = 2*p0;
            rho0 = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 200))),p0);
        end
        plow = p0/2;
        phigh = p0;
    elseif(rho0 > rho)
        while(rho0 > rho)
            p0 = p0/2;
            rho0 = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 200))),p0);
        end
        plow = p0;
        phigh = p0*2;
    else
        p = p0;  % although coincidence is imposible if it is not hard sphere.
        return;
    end
    pcen = (phigh + plow)/2;
    delt = phigh - plow;
    err = delt/pcen;
    % binary search
    while(err > tol)
        rho0 = beta/deriv(@(bp)(-log(Pfunc_isobaric_3NN(bp, beta, coeffs, 200))),pcen);
        if(rho0 < rho)
            plow = pcen;
        else
            phigh = pcen;
        end
        pcen = (phigh + plow)/2;
        delt = phigh - plow;
        err = delt/pcen;
    end
    p = pcen;
end

